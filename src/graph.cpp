#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "graph.h"
#include "postings.h"

#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>
#include <cstring>
#include <thread>

constexpr int NSHARD = 64;

static void bump(cov_t& c) { c += (c != COV_MAX); }

static void addsat(cov_t& d, cov_t v) {
    const uint32_t t = (uint32_t)d + v;
    d = t > COV_MAX ? COV_MAX : (cov_t)t;
}

static std::vector<std::string> collect_shards(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) die("graph", "no input %s", path.c_str());
    if (is_gz(path.c_str())) die("graph", "%s is compressed; snape reads uncompressed FASTQ only -- decompress it first", path.c_str());
    if (!S_ISDIR(st.st_mode)) return {path};
    std::vector<std::string> v;
    size_t skipped_gz = 0;
    DIR* d = opendir(path.c_str());
    if (!d) die("graph", "opendir %s", path.c_str());
    for (struct dirent* e = readdir(d); e; e = readdir(d)) {
        const std::string n = e->d_name;
        if (n.find(".fastq") == std::string::npos && n.find(".fq") == std::string::npos && n.find(".svb") == std::string::npos) continue;
        if (is_gz(n.c_str())) { skipped_gz++; continue; }
        v.push_back(path + "/" + n);
    }
    closedir(d);
    if (skipped_gz) say("graph", "%zu compressed shard(s) in %s ignored (uncompressed FASTQ only)", skipped_gz, path.c_str());
    if (v.empty()) die("graph", "no uncompressed .fastq/.fq/.svb shards in %s", path.c_str());
    std::sort(v.begin(), v.end());
    return v;
}

template <class T>
struct SortJob {
    std::vector<T>* v;
    size_t lo, hi;
    bool (*cmp)(const T&, const T&);
};

template <class T>
static void sort_job(SortJob<T>* j) {
    std::sort(j->v->begin() + j->lo, j->v->begin() + j->hi, j->cmp);
}

template <class T>
struct MergeJob {
    std::vector<T>* v;
    size_t a, b, e;
    bool (*cmp)(const T&, const T&);
};

template <class T>
static void merge_job(MergeJob<T>* j) {
    std::inplace_merge(j->v->begin() + j->a, j->v->begin() + j->b, j->v->begin() + j->e, j->cmp);
}

template <class T>
static void psort(std::vector<T>& v, int threads, bool (*cmp)(const T&, const T&)) {
    const size_t n = v.size();
    int nt = 1;
    while (nt * 2 <= threads && (size_t)nt * 2 <= n / 65536) nt *= 2;
    if (nt < 2) { std::sort(v.begin(), v.end(), cmp); return; }
    std::vector<size_t> bound(nt + 1);
    for (int i = 0; i <= nt; i++) bound[i] = n * i / nt;
    {
        std::vector<SortJob<T>> jobs(nt);
        for (int i = 0; i < nt; i++) jobs[i] = {&v, bound[i], bound[i + 1], cmp};
        std::vector<std::thread> th;
        for (int i = 0; i < nt; i++) th.emplace_back(sort_job<T>, &jobs[i]);
        for (auto& t : th) t.join();
    }
    for (int w = 1; w < nt; w *= 2) {
        const int nm = (nt - w) / (2 * w);
        std::vector<MergeJob<T>> jobs(nm);
        std::vector<std::thread> th;
        for (int i = 0, k = 0; i + w < nt; i += 2 * w)
            jobs[k++] = {&v, bound[i], bound[i + w], bound[std::min(i + 2 * w, nt)], cmp};
        for (int k = 0; k < nm; k++) th.emplace_back(merge_job<T>, &jobs[k]);
        for (auto& t : th) t.join();
    }
}

static FILE* open_fastq(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) die("graph", "fopen %s (input must be uncompressed FASTQ)", path.c_str());
    return f;
}

static FILE* open_svb(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) die("graph", "fopen svb %s", path.c_str());
    if (std::fseek(f, 16, SEEK_SET) != 0) die("graph", "bad .svb %s", path.c_str());
    return f;
}

static void count_chain(int n, const uint64_t* st, const uint64_t* en, const uint64_t* sg,
                         uint64_t len, uint64_t mg, SigCntMap& lo, uint64_t& ncores) {
    for (int i = 0; i < n; i++) {
        if (st[i] < mg || en[i] + mg > len) continue;
        ncores++;
        bump(lo[sg[i]]);
    }
}

static void count_read(std::string& s, int level, int margin, SigCntMap& lo, uint64_t& ncores) {
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int nf, nr;
    const uint64_t *stf, *enf, *sgf, *strv, *enr, *sgr;
    extract2(s.data(), (int)len, level, &nf, &stf, &enf, &sgf, &nr, &strv, &enr, &sgr);
    const uint64_t mg = (uint64_t)margin;
    count_chain(nf, stf, enf, sgf, len, mg, lo, ncores);
    count_chain(nr, strv, enr, sgr, len, mg, lo, ncores);
}

void SigIndex::build(const std::vector<SigCount>& nodes) {
    v = nodes.data();
    off.assign((size_t(1) << BITS) + 1, 0);
    for (const SigCount& c : nodes) off[(size_t)(c.sig >> (64 - BITS)) + 1]++;
    for (size_t i = 1; i < off.size(); i++) off[i] += off[i - 1];
}

uint32_t SigIndex::find(uint64_t sig) const {
    const size_t b = (size_t)(sig >> (64 - BITS));
    for (uint32_t i = off[b], e = off[b + 1]; i < e; i++) if (v[i].sig == sig) return i;
    return UINT32_MAX;
}

size_t SigIndex::bytes() const { return off.capacity() * 4; }

size_t PostIdx::bytes() const { return w.capacity() * 8 + cnt.capacity() + worst.capacity() * 8; }

static inline uint64_t post_key(uint64_t w0) { return w0 & 0x3ffffffffull; }
constexpr uint64_t POST_READ_MAX = 1ull << 28;
constexpr uint64_t POST_ORD_MAX  = 1ull << 14;
constexpr uint64_t POST_LEN_MAX  = 1ull << 16;

static inline uint64_t post_w0(uint32_t samp, uint64_t read, uint32_t ord, uint32_t clen) {
    return read | (uint64_t)samp << 28 | (uint64_t)ord << 34 | (uint64_t)clen << 48;
}

static inline uint64_t post_w1(uint32_t st, uint32_t dl, uint32_t dr) {
    return (uint64_t)st | (uint64_t)dl << 32 | (uint64_t)dr << 48;
}

static void post_offer(PostIdx& px, uint32_t node, uint64_t w0, uint64_t w1) {
    const uint64_t key = post_key(w0);
    if (px.cnt[node] >= px.N && key >= px.worst[node]) return;
    std::lock_guard<std::mutex> lk(px.mtx[node & (NSHARD - 1)]);
    uint64_t* row = px.w.data() + (size_t)node * (size_t)px.N * 2;
    uint8_t& c = px.cnt[node];
    size_t j;
    if (c < px.N) {
        j = c;
    } else {
        if (key >= post_key(row[(size_t)(px.N - 1) * 2])) return;
        j = px.N - 1;
    }
    while (j > 0 && post_key(row[(j - 1) * 2]) > key) {
        row[j * 2] = row[(j - 1) * 2];
        row[j * 2 + 1] = row[(j - 1) * 2 + 1];
        j--;
    }
    row[j * 2] = w0;
    row[j * 2 + 1] = w1;
    if (c < px.N) c++;
    if (c == px.N) px.worst[node] = post_key(row[(size_t)(px.N - 1) * 2]);
}

static std::atomic<uint64_t> POST_SKIP{0};

static void post_write(const PostIdx& px, const std::string& path, uint64_t K) {
    std::vector<uint32_t> off(K + 1, 0);
    for (uint64_t i = 0; i < K; i++) off[i + 1] = off[i] + px.cnt[i];
    const uint64_t E = off[K];
    FILE* f = xopen((path + ".tmp").c_str(), "wb");
    wr(f, &px.pad, 4);
    wr(f, &px.N, 4);
    uint64_t Kw = K;
    wr(f, &Kw, 8);
    wr(f, off.data(), (K + 1) * 4);
    std::vector<uint64_t> buf;
    buf.reserve(1u << 20);
    for (uint64_t i = 0; i < K; i++) {
        const uint64_t* row = px.w.data() + i * (size_t)px.N * 2;
        for (uint32_t j = 0; j < px.cnt[i]; j++) {
            buf.push_back(row[j * 2]);
            buf.push_back(row[j * 2 + 1]);
        }
        if (buf.size() >= (1u << 20)) { wr(f, buf.data(), buf.size() * 8); buf.clear(); }
    }
    if (!buf.empty()) wr(f, buf.data(), buf.size() * 8);
    atomic_close(f, (path + ".tmp").c_str(), path.c_str());
    uint64_t full = 0;
    for (uint64_t i = 0; i < K; i++) full += (px.cnt[i] >= px.N);
    say("index", "records skipped (no room in slot): %llu", (unsigned long long)POST_SKIP.load());
    say("index", "%llu records, %.3f M nodes %.1f%% full (N=%u), %.2f GB -> %s",
        (unsigned long long)E, K / 1e6, K ? 100.0 * (double)full / (double)K : 0.0, px.N,
        (16.0 * (double)E + 4.0 * (double)(K + 1) + 16.0) / 1073741824.0, path.c_str());
}

static bool feed_take(Feed& fd, Batch& out) {
    std::unique_lock<std::mutex> lk(fd.qmtx);
    while (fd.q.empty() && fd.readers_left.load() != 0) fd.cv_empty.wait(lk);
    if (fd.q.empty()) return false;
    out = std::move(fd.q.front());
    fd.q.pop();
    fd.cv_full.notify_one();
    return true;
}

static void feed_push(Feed& fd, Batch& batch, size_t q_cap) {
    if (batch.seqs.empty()) return;
    std::unique_lock<std::mutex> lk(fd.qmtx);
    while (fd.q.size() >= q_cap) fd.cv_full.wait(lk);
    fd.q.push(std::move(batch));
    fd.cv_empty.notify_one();
}

struct ReaderJob {
    const std::vector<Shard>* shards;
    Feed* fd;
    size_t q_cap, batch_n;
};

static void reader_thread(ReaderJob* j) {
    constexpr unsigned BLOCK = 1u << 20;
    std::vector<char> buf(BLOCK);
    Batch batch;
    batch.seqs.reserve(j->batch_n);
    std::string carry;
    std::vector<uint64_t> mine(j->fd->reads_per.size(), 0);
    size_t si;
    std::string svseq;
    while ((si = j->fd->next.fetch_add(1)) < j->shards->size()) {
        const Shard& sh = (*j->shards)[si];
        const int col = sh.col;
        if (!batch.seqs.empty() && batch.col != col) feed_push(*j->fd, batch, j->q_cap);
        batch.col = col;
        if (batch.seqs.empty()) batch.first = mine[col];
        const bool svb = is_svb(sh.path.c_str());
        FILE* f = svb ? open_svb(sh.path) : open_fastq(sh.path);
        uint64_t bases = 0;
        if (svb) {
            while (svb_read_seq(f, svseq)) {
                if (batch.seqs.empty()) batch.first = mine[col];
                bases += svseq.size();
                mine[col]++;
                batch.seqs.emplace_back(std::move(svseq));
                if (batch.seqs.size() >= j->batch_n) feed_push(*j->fd, batch, j->q_cap);
                if (j->fd->max_bases && bases >= j->fd->max_bases) break;
            }
        } else {
            int line_no = 0;
            carry.clear();
            int got;
            bool done = false;
            while (!done && (got = (int)std::fread(buf.data(), 1, BLOCK, f)) > 0) {
                const char* base = buf.data();
                const char* end = base + got;
                const char* start = base;
                for (const char* p = base; p < end && !done; p++) {
                    if (*p != '\n') continue;
                    const size_t n = (size_t)(p - start);
                    if (line_no == 1) {
                        if (batch.seqs.empty()) batch.first = mine[col];
                        if (carry.empty()) batch.seqs.emplace_back(start, n);
                        else { carry.append(start, n); batch.seqs.push_back(std::move(carry)); carry.clear(); }
                        bases += batch.seqs.back().size();
                        mine[col]++;
                        if (batch.seqs.size() >= j->batch_n) feed_push(*j->fd, batch, j->q_cap);
                        done = j->fd->max_bases && bases >= j->fd->max_bases;
                    } else if (!carry.empty()) carry.clear();
                    line_no = (line_no + 1) & 3;
                    start = p + 1;
                }
                if (start < end && line_no == 1) carry.append(start, end - start);
            }
        }
        std::fclose(f);
        feed_push(*j->fd, batch, j->q_cap);
    }
    feed_push(*j->fd, batch, j->q_cap);
    {
        std::lock_guard<std::mutex> lk(j->fd->qmtx);
        for (size_t i = 0; i < mine.size(); i++) j->fd->reads_per[i] += mine[i];
    }
    if (j->fd->readers_left.fetch_sub(1) == 1) {
        std::lock_guard<std::mutex> lk(j->fd->qmtx);
        j->fd->cv_empty.notify_all();
    }
}

struct CountW {
    Feed* fd;
    int level, margin;
    ShTab* psn;
    std::mutex* nmtx;
    size_t flush;
    std::atomic<uint64_t>* emissions;
};

static void count_flush(CountW& w, SigCntMap& lo, std::vector<std::vector<std::pair<uint64_t, cov_t>>>& bn) {
    for (auto& v : bn) v.clear();
    for (const auto& kv : lo) bn[kv.first & (NSHARD - 1)].push_back(kv);
    for (int sh = 0; sh < NSHARD; sh++) {
        if (bn[sh].empty()) continue;
        std::lock_guard<std::mutex> lk(w.nmtx[sh]);
        for (const auto& kv : bn[sh]) addsat(w.psn[sh].m[kv.first], kv.second);
    }
    lo.clear();
}

static void count_worker(CountW* w) {
    uint64_t emis = 0;
    SigCntMap lo;
    lo.reserve(1u << 20);
    std::vector<std::vector<std::pair<uint64_t, cov_t>>> bn(NSHARD);
    Batch b;
    while (feed_take(*w->fd, b)) {
        for (std::string& s : b.seqs) count_read(s, w->level, w->margin, lo, emis);
        if (lo.size() >= w->flush) count_flush(*w, lo, bn);
    }
    count_flush(*w, lo, bn);
    *w->emissions += emis;
}

static void count_sample(const std::vector<Shard>& shards, int level, int margin, int threads, uint64_t max_bases,
                         ShTab* psn, std::mutex* nmtx,
                         std::atomic<uint64_t>& emissions, uint64_t& reads_here) {
    const int nw = std::max(1, threads - 1);
    Feed fd;
    fd.reads_per.assign(1, 0);
    fd.max_bases = max_bases;
    std::vector<Shard> mine = shards;
    for (Shard& s : mine) s.col = 0;

    ReaderJob rj{&mine, &fd, 256, 64};
    CountW w{&fd, level, margin, psn, nmtx, flush_env(), &emissions};
    std::vector<std::thread> th;
    th.emplace_back(reader_thread, &rj);
    for (int i = 0; i < nw; i++) th.emplace_back(count_worker, &w);
    for (auto& t : th) t.join();
    reads_here = fd.reads_per[0];
}

static cov_t* node_slot(NodeTab& t, uint64_t k, int S) {
    const auto r = t.m.try_emplace(k, (uint32_t)(t.a.size() / (size_t)S));
    if (r.second) t.a.resize(t.a.size() + (size_t)S, 0);
    return t.a.data() + (size_t)r.first->second * (size_t)S;
}

struct MergeW {
    ShTab* psn;
    NodeTab* gnode;
    int col, S, min_mult;
    std::atomic<uint64_t>* nseen, *npass, *n1, *n2;
    std::atomic<int>* next_sh;
};

static void merge_shard(MergeW* w, int sh) {
    uint64_t a = 0, b = 0, c1 = 0, c2 = 0;
    for (const auto& kv : w->psn[sh].m) {
        a++;
        c1 += (kv.second == 1);
        c2 += (kv.second == 2);
        if ((uint32_t)kv.second < (uint32_t)w->min_mult) continue;
        b++;
        node_slot(w->gnode[sh], kv.first, w->S)[w->col] = kv.second;
    }
    *w->nseen += a;
    *w->npass += b;
    *w->n1 += c1;
    *w->n2 += c2;
    SigCntMap().swap(w->psn[sh].m);
}

static void merge_worker(MergeW* w) {
    for (int sh = w->next_sh->fetch_add(1); sh < NSHARD; sh = w->next_sh->fetch_add(1))
        merge_shard(w, sh);
}

static void merge_sample(ShTab* psn, NodeTab* gnode, int col, int S, int min_mult, int threads,
                         uint64_t& seen, uint64_t& kept, uint64_t& one, uint64_t& two) {
    const int nt = std::max(1, threads);
    std::atomic<uint64_t> nseen{0}, npass{0}, n1{0}, n2{0};
    std::atomic<int> next_sh{0};
    MergeW w{psn, gnode, col, S, min_mult, &nseen, &npass, &n1, &n2, &next_sh};
    std::vector<std::thread> th;
    for (int i = 0; i < nt; i++) th.emplace_back(merge_worker, &w);
    for (auto& t : th) t.join();
    seen = nseen.load();
    kept = npass.load();
    one = n1.load();
    two = n2.load();
}

static void fill_chain(bool rev, int n, const uint64_t* st, const uint64_t* en, const uint64_t* sig,
                       uint64_t len, uint64_t mg, uint64_t rno, int col, const FillCtx& cx, LocalB& lo,
                       uint64_t& ncores) {
    static thread_local std::vector<int> ci;
    static thread_local std::vector<uint32_t> cnode;
    const bool wantpost = cx.px && cx.px->N != 0 && !rev;
    if (wantpost) { ci.clear(); cnode.clear(); }
    const int S = cx.S;
    const cov_t* cov = cx.cov->data();
    uint32_t prev = UINT32_MAX;
    bool have_prev = false;
    uint64_t prev_sig = 0;
    for (int i = 0; i < n; i++) {
        if (st[i] < mg || en[i] + mg > len) continue;
        ncores++;
        const uint32_t node = cx.idx->find(sig[i]);

        if (have_prev && prev_sig != sig[i] && prev != UINT32_MAX && node != UINT32_MAX &&
            prev != node && cov[(size_t)prev * S + col] != 0 && cov[(size_t)node * S + col] != 0) {
            lo.e.insert(((uint64_t)prev << 32) | (uint64_t)node);
        }
        if (wantpost) { ci.push_back(i); cnode.push_back(node); }
        prev = node;
        prev_sig = sig[i];
        have_prev = true;
    }
    if (!wantpost || ci.empty()) return;
    const int P = (int)cx.px->pad, m = (int)ci.size();
    for (int p = 0; p < m; p++) {
        const uint32_t node = cnode[p];
        if (node == UINT32_MAX) continue;
        const int la = p > P ? p - P : 0, rb = p + P < m ? p + P : m - 1;
        const uint64_t a = st[ci[p]], b = en[ci[p]];
        const uint64_t ws = st[ci[la]], we = en[ci[rb]];
        if (b <= a || b - a >= POST_LEN_MAX || a - ws >= POST_LEN_MAX ||
            we - b >= POST_LEN_MAX || rno >= POST_READ_MAX ||
            (uint64_t)p >= POST_ORD_MAX || a >= (1ull << 32)) {
            POST_SKIP.fetch_add(1, std::memory_order_relaxed);
            continue;
        }
        post_offer(*cx.px, node,
                   post_w0((uint32_t)col, rno, (uint32_t)p, (uint32_t)(b - a)),
                   post_w1((uint32_t)a, (uint32_t)(a - ws), (uint32_t)(we - b)));
    }
}

static void fill_read(std::string& s, int level, int margin, int col, uint64_t rno,
                      const FillCtx& cx, LocalB& lo, uint64_t& ncores) {
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int nf, nr;
    const uint64_t *stf, *enf, *sgf, *strv, *enr, *sgr;
    extract2(s.data(), (int)len, level, &nf, &stf, &enf, &sgf, &nr, &strv, &enr, &sgr);
    const uint64_t mg = (uint64_t)margin;
    fill_chain(false, nf, stf, enf, sgf, len, mg, rno, col, cx, lo, ncores);
    fill_chain(true, nr, strv, enr, sgr, len, mg, rno, col, cx, lo, ncores);
}

struct FillW {
    Feed* fd;
    int level, margin;
    const FillCtx* cx;
    EdgeKeySet* gedge;
    std::mutex* emtx;
    size_t flush;
    std::atomic<uint64_t>* emissions;
};

static void fill_flush(FillW& w, LocalB& lo, std::vector<std::vector<uint64_t>>& be) {
    for (auto& v : be) v.clear();
    for (uint64_t e : lo.e) be[(size_t)((e >> 32) ^ (e * 0x9e3779b97f4a7c15ULL)) & (NSHARD - 1)].push_back(e);
    for (int sh = 0; sh < NSHARD; sh++) {
        if (be[sh].empty()) continue;
        std::lock_guard<std::mutex> lk(w.emtx[sh]);
        for (uint64_t e : be[sh]) w.gedge[sh].insert(e);
    }
    lo.e.clear();
}

static void fill_worker(FillW* w) {
    uint64_t emis = 0;
    LocalB lo;
    lo.e.reserve(1u << 20);
    std::vector<std::vector<uint64_t>> be(NSHARD);
    Batch b;
    while (feed_take(*w->fd, b)) {
        for (size_t k = 0; k < b.seqs.size(); k++)
            fill_read(b.seqs[k], w->level, w->margin, b.col, b.first + k, *w->cx, lo, emis);
        if (lo.e.size() >= w->flush) fill_flush(*w, lo, be);
    }
    fill_flush(*w, lo, be);
    *w->emissions += emis;
}

static void fill_pass(const std::vector<Shard>& shards, int level, int margin, int threads, int S, uint64_t max_bases,
                      const SigIndex& idx, const std::vector<cov_t>& cov,
                      PostIdx& px, EdgeKeySet* gedge, std::mutex* emtx,
                      std::atomic<uint64_t>& emissions, std::vector<uint64_t>& reads_per) {
    const int nw = std::max(1, threads - 1);
    Feed fd;
    fd.reads_per.assign((size_t)S, 0);
    fd.max_bases = max_bases;
    FillCtx cx;
    cx.idx = &idx;
    cx.cov = &cov;
    cx.px = &px;
    cx.S = S;

    ReaderJob rj{&shards, &fd, 256, 64};
    FillW fw{&fd, level, margin, &cx, gedge, emtx, flush_env(), &emissions};
    std::vector<std::thread> th;
    th.emplace_back(reader_thread, &rj);
    for (int i = 0; i < nw; i++) th.emplace_back(fill_worker, &fw);
    for (auto& t : th) t.join();
    for (int i = 0; i < S; i++) reads_per[i] = fd.reads_per[i];
}

struct GatherJob {
    NodeTab* gnode;
    std::vector<std::vector<SigCount>>* part;
    std::vector<std::vector<cov_t>>* pcov;
    int sh, S;
};

static void gather_shard(GatherJob* j) {
    NodeTab& t = j->gnode[j->sh];
    std::vector<SigCount>& pv = (*j->part)[j->sh];
    std::vector<cov_t>& cv = (*j->pcov)[j->sh];
    pv.reserve(t.m.size());
    for (const auto& kv : t.m) {
        const cov_t* c = t.a.data() + (size_t)kv.second * (size_t)j->S;
        pv.push_back({kv.first, (uint32_t)(cv.size() / (size_t)j->S)});
        cv.insert(cv.end(), c, c + j->S);
    }
    SigNodeMap().swap(t.m);
    std::vector<cov_t>().swap(t.a);
}

static void gather_worker(std::atomic<int>* next_sh, std::vector<GatherJob>* jobs) {
    for (int sh = next_sh->fetch_add(1); sh < NSHARD; sh = next_sh->fetch_add(1))
        gather_shard(&(*jobs)[sh]);
}

static bool sig_cmp(const SigCount& a, const SigCount& b) { return a.sig < b.sig; }

static bool edge_cmp(const CEdge& a, const CEdge& b) {
    return a.src != b.src ? a.src < b.src : a.dst < b.dst;
}

static void run_pool_two_pass(const std::vector<Shard>& shards, int level, int margin, int threads, uint64_t max_bases,
                              int min_mult, int S, int index_n, int pad_cores, const std::string& dir,
                              std::vector<SigCount>& out_nodes,
                              std::vector<cov_t>& out_ncov, std::vector<CEdge>& out_edges,
                              std::vector<uint64_t>& reads_per) {
    NodeTab gnode[NSHARD];
    ShTab psn[NSHARD];
    std::mutex nmtx[NSHARD], emtx[NSHARD], smtx[NSHARD];
    std::atomic<uint64_t> emis_a{0};
    reads_per.assign((size_t)S, 0);
    uint64_t distinct = 0, kept = 0;

    for (int col = 0; col < S; col++) {
        std::vector<Shard> mine;
        for (const Shard& s : shards) if (s.col == col) mine.push_back(s);
        uint64_t rd = 0, seen = 0, pass = 0, one = 0, two = 0;
        count_sample(mine, level, margin, threads, max_bases, psn, nmtx, emis_a, rd);
        merge_sample(psn, gnode, col, S, min_mult, threads, seen, pass, one, two);
        reads_per[col] = rd;
        distinct += seen;
        kept += pass;
        say("graphA", "sample %d/%d: %llu reads, %.3f M distinct -> %.3f M >= %d (%.1f%%)",
            col + 1, S, (unsigned long long)rd, seen / 1e6, pass / 1e6, min_mult,
            seen ? 100.0 * pass / seen : 0.0);
        say("graphA", "  distinct seen once %.3f M (%.1f%%), twice %.3f M (%.1f%%)",
            one / 1e6, seen ? 100.0 * one / seen : 0.0,
            two / 1e6, seen ? 100.0 * two / seen : 0.0);
    }

    {
        std::vector<std::vector<SigCount>> part(NSHARD);
        std::vector<std::vector<cov_t>> pcov(NSHARD);
        std::atomic<int> next_sh{0};
        std::vector<GatherJob> jobs(NSHARD);
        for (int i = 0; i < NSHARD; i++)
            jobs[i] = GatherJob{gnode, &part, &pcov, i, S};
        const int nt = std::max(1, threads);
        std::vector<std::thread> th;
        for (int i = 0; i < nt; i++) th.emplace_back(gather_worker, &next_sh, &jobs);
        for (auto& t : th) t.join();
        size_t K = 0;
        for (int sh = 0; sh < NSHARD; sh++) K += part[sh].size();
        out_nodes.clear();
        out_nodes.reserve(K);
        std::vector<cov_t> tmpcov;
        tmpcov.reserve(K * (size_t)S);
        for (int sh = 0; sh < NSHARD; sh++) {
            const uint32_t base = (uint32_t)(tmpcov.size() / (size_t)S);
            for (SigCount c : part[sh]) { c.slot += base; out_nodes.push_back(c); }
            tmpcov.insert(tmpcov.end(), pcov[sh].begin(), pcov[sh].end());
            std::vector<SigCount>().swap(part[sh]);
            std::vector<cov_t>().swap(pcov[sh]);
        }
        psort(out_nodes, threads, sig_cmp);
        out_ncov.assign(K * (size_t)S, 0);
        for (size_t i = 0; i < K; i++) {
            const cov_t* c = tmpcov.data() + (size_t)out_nodes[i].slot * (size_t)S;
            std::copy(c, c + S, out_ncov.begin() + (std::ptrdiff_t)(i * (size_t)S));
            out_nodes[i].slot = (uint32_t)i;
        }
    }

    const size_t K = out_nodes.size();
    SigIndex idx;
    idx.build(out_nodes);
    say("graphA", "survivors=%.3f M, index %.2f GB -> pass B",
        K / 1e6, (idx.bytes() + K * 8) / 1073741824.0);

    EdgeKeySet gedge[NSHARD];
    std::atomic<uint64_t> emis_b{0};
    std::vector<uint64_t> rp_b((size_t)S, 0);

    PostIdx px;
    if (index_n > 0) {
        px.N = (uint32_t)index_n;
        px.pad = (uint32_t)pad_cores;
        px.w.assign(K * (size_t)index_n * 2, 0);
        px.cnt.assign(K, 0);
        px.worst.assign(K, 0);
        px.mtx = smtx;
        say("index", "N=%d pad=%d, buffer %.2f GB",
            index_n, pad_cores, px.bytes() / 1073741824.0);
    }
    fill_pass(shards, level, margin, threads, S, max_bases, idx, out_ncov, px,
              gedge, emtx, emis_b, rp_b);
    if (px.N != 0) {
        post_write(px, dir + "/index.bin", K);
        std::vector<uint64_t>().swap(px.w);
        std::vector<uint8_t>().swap(px.cnt);
        std::vector<uint64_t>().swap(px.worst);
    }

    {
        uint64_t ne = 0;
        for (int sh = 0; sh < NSHARD; sh++) ne += gedge[sh].size();
        out_edges.clear();
        out_edges.reserve(ne);
        for (int sh = 0; sh < NSHARD; sh++) {
            for (uint64_t e : gedge[sh]) out_edges.push_back({(uint32_t)(e >> 32), (uint32_t)e});
            EdgeKeySet().swap(gedge[sh]);
        }
    }
    say("graphA", "emissions=%.3f M", emis_a.load() / 1e6);
    say("graphB", "emissions=%.3f M", emis_b.load() / 1e6);
    say("graph1", "per-sample distinct sigs=%.3f M -> passed mm=%d: %.3f M (%.2f%%)",
        distinct / 1e6, min_mult, kept / 1e6, distinct ? 100.0 * kept / distinct : 0.0);
    say("graph1", "merged nodes=%.3f M | edges=%.3f M",
        K / 1e6, out_edges.size() / 1e6);
}

void run_graph_reads(const std::vector<Sample>& samples, const std::string& dir,
                     int level, int margin, int min_mult, int threads, uint64_t max_bases, int index_n, int pad_cores) {
    mkdir(dir.c_str(), 0755);
    const int S = (int)samples.size();
    std::vector<Shard> shards_all;
    for (int i = 0; i < S; i++)
        for (const auto& p : collect_shards(samples[i].path)) shards_all.push_back({p, i});
    if (max_bases && shards_all.size() != (size_t)S) die("graph", "--max-bases needs one file per sample");
    say("graph", "%d samples (%zu shards), threads=%d min_mult=%d max_bases=%llu",
        S, shards_all.size(), threads, min_mult, (unsigned long long)max_bases);

    std::vector<SigCount> nodes;
    std::vector<CEdge> ce;
    std::vector<cov_t> ncov;
    std::vector<uint64_t> nreads;
    run_pool_two_pass(shards_all, level, margin, threads, max_bases, min_mult, S, index_n, pad_cores, dir,
                      nodes, ncov, ce, nreads);
    const uint64_t K = nodes.size();
    say("graph", "  survivors=%lu", (unsigned long)K);

    psort(ce, threads, edge_cmp);

    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen((gp + ".tmp").c_str(), "wb");
    const uint32_t Sw = (uint32_t)S;
    wr(fg, &Sw, 4);
    uint64_t Kw = K;
    wr(fg, &Kw, 8);
    {
        std::vector<uint64_t> sig(K);
        for (uint64_t i = 0; i < K; i++) sig[i] = nodes[i].sig;
        wr(fg, sig.data(), K * 8);
    }
    std::vector<cov_t> plane(K);
    for (int j = 0; j < S; j++) {
        for (uint64_t i = 0; i < K; i++) plane[i] = ncov[(size_t)nodes[i].slot * (size_t)S + j];
        wr(fg, plane.data(), K * sizeof(cov_t));
    }

    std::vector<SigCount>().swap(nodes);
    std::vector<cov_t>().swap(ncov);
    const uint64_t E = ce.size();
    wr(fg, &E, 8);
    {
        std::vector<uint32_t> col(E);
        for (uint64_t i = 0; i < E; i++) col[i] = ce[i].src;
        wr(fg, col.data(), E * 4);
        for (uint64_t i = 0; i < E; i++) col[i] = ce[i].dst;
        wr(fg, col.data(), E * 4);
    }
    atomic_close(fg, (gp + ".tmp").c_str(), gp.c_str());
    for (int i = 0; i < S; i++)
        say("graph", "  %-8s %-6s %lu reads", samples[i].name.c_str(),
            samples[i].normal ? "normal" : "tumor", (unsigned long)nreads[i]);
    say("graph", "nodes=%lu edges=%lu S=%d -> %s/graph.bin",
        (unsigned long)K, (unsigned long)E, S, dir.c_str());
}
