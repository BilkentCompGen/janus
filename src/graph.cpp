#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "graph.h"
#include "samples.h"

#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include "unordered_dense.h"

#include <utility>
#include <vector>

FILE* xopen(const std::string& p, const char* m) { FILE* f = std::fopen(p.c_str(), m); if (!f) { std::fprintf(stderr, "[FATAL] fopen %s\n", p.c_str()); std::abort(); } return f; }
void wr(FILE* f, const void* p, size_t n) { if (n && std::fwrite(p, 1, n, f) != n) { std::fprintf(stderr, "[FATAL] fwrite\n"); std::abort(); } }
void atomic_close(FILE* f, const std::string& tmp, const std::string& dst) { std::fclose(f); if (std::rename(tmp.c_str(), dst.c_str())) { std::fprintf(stderr, "[FATAL] rename %s\n", dst.c_str()); std::abort(); } }

namespace {
static inline void bump(cov_t& c) { c += (c != COV_MAX); }

static inline void addsat(cov_t& d, cov_t v) { const uint32_t t = (uint32_t)d + v; d = t > COV_MAX ? COV_MAX : (cov_t)t; }

static bool is_gz(const std::string& p) { return p.size() >= 3 && p.compare(p.size() - 3, 3, ".gz") == 0; }

std::vector<std::string> collect_shards(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) { std::fprintf(stderr, "[FATAL] no input %s\n", path.c_str()); std::abort(); }
    if (is_gz(path)) { std::fprintf(stderr, "[FATAL] %s is compressed; snape reads uncompressed FASTQ only -- decompress it first\n", path.c_str()); std::abort(); }
    if (!S_ISDIR(st.st_mode)) return {path};
    std::vector<std::string> v;
    size_t skipped_gz = 0;
    DIR* d = opendir(path.c_str());
    if (!d) { std::fprintf(stderr, "[FATAL] opendir %s\n", path.c_str()); std::abort(); }
    for (struct dirent* e = readdir(d); e; e = readdir(d)) {
        const std::string n = e->d_name;
        if (n.find(".fastq") == std::string::npos && n.find(".fq") == std::string::npos && n.find(".svb") == std::string::npos) continue;
        if (is_gz(n)) { skipped_gz++; continue; }
        v.push_back(path + "/" + n);
    }
    closedir(d);
    if (skipped_gz) std::fprintf(stderr, "[graph] %zu compressed shard(s) in %s ignored (uncompressed FASTQ only)\n", skipped_gz, path.c_str());
    if (v.empty()) { std::fprintf(stderr, "[FATAL] no uncompressed .fastq/.fq/.svb shards in %s\n", path.c_str()); std::abort(); }
    std::sort(v.begin(), v.end());
    return v;
}

template <class T, class Cmp>
void psort(std::vector<T>& v, int threads, Cmp cmp) {
    const size_t n = v.size();
    int nt = 1; while (nt * 2 <= threads && (size_t)nt * 2 <= n / 65536) nt *= 2;
    if (nt < 2) { std::sort(v.begin(), v.end(), cmp); return; }
    std::vector<size_t> bound(nt + 1);
    for (int i = 0; i <= nt; i++) bound[i] = n * i / nt;
    { std::vector<std::thread> th;
      for (int i = 0; i < nt; i++) th.emplace_back([&, i]{ std::sort(v.begin() + bound[i], v.begin() + bound[i + 1], cmp); });
      for (auto& t : th) t.join(); }
    for (int w = 1; w < nt; w *= 2) {
        std::vector<std::thread> th;
        for (int i = 0; i + w < nt; i += 2 * w)
            th.emplace_back([&, i, w]{
                std::inplace_merge(v.begin() + bound[i], v.begin() + bound[i + w],
                                   v.begin() + bound[std::min(i + 2 * w, nt)], cmp);
            });
        for (auto& t : th) t.join();
    }
}

struct Shard { std::string path; int col; };
struct Batch { std::vector<std::string> seqs; int col = 0; uint64_t first = 0; };

static FILE* open_fastq(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[FATAL] fopen %s (input must be uncompressed FASTQ)\n", path.c_str()); std::abort(); }
    return f;
}

static bool is_svb(const std::string& p) { return p.size() >= 4 && p.compare(p.size() - 4, 4, ".svb") == 0; }

static const uint32_t* quad_table() {
    static uint32_t q[256];
    static const bool once = []{
        static const char L[4] = { 'A', 'C', 'G', 'T' };
        for (int b = 0; b < 256; b++) {
            const char c[4] = { L[(b >> 6) & 3], L[(b >> 4) & 3], L[(b >> 2) & 3], L[b & 3] };
            std::memcpy(&q[b], c, 4);
        }
        return true;
    }();
    (void)once;
    return q;
}

static FILE* open_svb(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[FATAL] fopen svb %s\n", path.c_str()); std::abort(); }
    if (std::fseek(f, 16, SEEK_SET) != 0) { std::fprintf(stderr, "[FATAL] bad .svb %s\n", path.c_str()); std::abort(); }
    return f;
}

static bool svb_next(FILE* f, std::string& seq) {
    uint32_t tag;
    if (std::fread(&tag, 4, 1, f) != 1) return false;
    const uint32_t L = tag & 0x7fffffffu;
    seq.resize(L);
    if (tag >> 31) {
        if (L && std::fread(&seq[0], 1, L, f) != L) { std::fprintf(stderr, "[FATAL] .svb raw truncated\n"); std::abort(); }
        return true;
    }
    const uint32_t nb = (L + 3) >> 2;
    static thread_local std::vector<uint8_t> pb;
    pb.resize(nb);
    if (nb && std::fread(pb.data(), 1, nb, f) != nb) { std::fprintf(stderr, "[FATAL] .svb packed truncated\n"); std::abort(); }

    static const uint32_t* QUAD = quad_table();
    static const char LUT[4] = { 'A', 'C', 'G', 'T' };
    const uint32_t full = L >> 2;
    char* out = &seq[0];
    for (uint32_t i = 0; i < full; i++) std::memcpy(out + 4 * i, &QUAD[pb[i]], 4);
    for (uint32_t i = full << 2; i < L; i++) seq[i] = LUT[(pb[i >> 2] >> (6 - 2 * (i & 3))) & 3];
    return true;
}

}

static bool rd_line(Reads& r, const char*& s, size_t& n) {
    r.carry.clear();
    for (;;) {
        if (r.pos >= r.end) {
            r.end = std::fread(r.buf.data(), 1, r.buf.size(), r.f);
            r.pos = 0;
            if (r.end == 0) {
                if (r.carry.empty()) return false;
                s = r.carry.data(); n = r.carry.size();
                while (n && (s[n - 1] == '\n' || s[n - 1] == '\r')) n--;
                return true;
            }
        }
        const char* base = r.buf.data() + r.pos;
        const char* nl = (const char*)std::memchr(base, '\n', r.end - r.pos);
        if (!nl) { r.carry.append(base, r.end - r.pos); r.pos = r.end; continue; }
        size_t len = (size_t)(nl - base);
        r.pos += len + 1;
        if (r.carry.empty()) { s = base; n = len; }
        else { r.carry.append(base, len); s = r.carry.data(); n = r.carry.size(); }
        while (n && (s[n - 1] == '\r')) n--;
        return true;
    }
}

Reads reads_open(const std::string& path, const std::string& sample) {
    Reads r;
    r.svb = is_svb(path);
    r.f = r.svb ? open_svb(path) : open_fastq(path);
    if (!r.svb) r.buf.resize(1u << 20);
    if (!sample.empty()) r.tag = sample + "|";
    return r;
}

bool reads_next(Reads& r, std::string& seq, std::string* name) {
    if (r.svb) {
        if (!svb_next(r.f, seq)) return false;
        if (name) { char b[32]; std::snprintf(b, sizeof b, "r%llu", (unsigned long long)r.idx); name->assign(r.tag); name->append(b); }
        r.idx++;
        return true;
    }
    const char* p; size_t n;
    if (!rd_line(r, p, n)) return false;
    if (name) {
        const char* sp = n > 1 ? (const char*)std::memchr(p + 1, ' ', n - 1) : nullptr;
        name->assign(r.tag);
        name->append(p + 1, sp ? (size_t)(sp - p - 1) : (n > 0 ? n - 1 : 0));
    }
    if (!rd_line(r, p, n)) return false;
    seq.assign(p, n);
    if (!rd_line(r, p, n)) return false;
    if (!rd_line(r, p, n)) return false;
    r.idx++;
    return true;
}

void reads_close(Reads& r) { if (r.f) std::fclose(r.f); r.f = nullptr; }

namespace {

using SigNodeMap = ankerl::unordered_dense::map<uint64_t, uint32_t>;
struct NodeTab { SigNodeMap m; std::vector<cov_t> a; };

static inline cov_t* node_slot(NodeTab& t, uint64_t k, int S) {
    const auto r = t.m.try_emplace(k, (uint32_t)(t.a.size() / (size_t)S));
    if (r.second) t.a.resize(t.a.size() + (size_t)S, 0);
    return t.a.data() + (size_t)r.first->second * (size_t)S;
}

using SigCntMap = ankerl::unordered_dense::map<uint64_t, cov_t>;
struct ShTab { SigCntMap m; };
using EdgeKeySet = ankerl::unordered_dense::set<uint64_t>;

constexpr int NSHARD = 64;



void count_read(std::string& s, int level, int margin, SigCntMap& lo, uint64_t& ncores) {
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int nf, nr; const uint64_t *stf, *enf, *sgf, *strv, *enr, *sgr;
    extract2(s.data(), (int)len, level, &nf, &stf, &enf, &sgf, &nr, &strv, &enr, &sgr);
    const uint64_t mg = (uint64_t)margin;
    auto chain = [&](int n, const uint64_t* st, const uint64_t* en, const uint64_t* sig) {
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            ncores++;
            bump(lo[sig[i]]);
        }
    };
    chain(nf, stf, enf, sgf);
    chain(nr, strv, enr, sgr);
}

struct SigIndex {
    static constexpr int BITS = 24;
    std::vector<uint32_t> off;
    const SigCount* v = nullptr;
    void build(const std::vector<SigCount>& nodes) {
        v = nodes.data();
        off.assign((size_t(1) << BITS) + 1, 0);
        for (const SigCount& c : nodes) off[(size_t)(c.sig >> (64 - BITS)) + 1]++;
        for (size_t i = 1; i < off.size(); i++) off[i] += off[i - 1];
    }
    uint32_t find(uint64_t sig) const {
        const size_t b = (size_t)(sig >> (64 - BITS));
        for (uint32_t i = off[b], e = off[b + 1]; i < e; i++) if (v[i].sig == sig) return i;
        return UINT32_MAX;
    }
    size_t bytes() const { return off.capacity() * 4; }
};

std::atomic<uint64_t> POST_SKIP{0};

struct PostIdx {
    uint32_t N = 0, pad = 0;
    std::vector<uint64_t> w;
    std::vector<uint8_t> cnt;
    std::vector<uint64_t> worst;
    std::mutex* mtx = nullptr;
    bool on() const { return N != 0; }
    size_t bytes() const { return w.capacity() * 8 + cnt.capacity() + worst.capacity() * 8; }
};

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

void post_offer(PostIdx& px, uint32_t node, uint64_t w0, uint64_t w1) {
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
        row[j * 2] = row[(j - 1) * 2]; row[j * 2 + 1] = row[(j - 1) * 2 + 1]; j--;
    }
    row[j * 2] = w0; row[j * 2 + 1] = w1;
    if (c < px.N) c++;
    if (c == px.N) px.worst[node] = post_key(row[(size_t)(px.N - 1) * 2]);
}

void post_write(const PostIdx& px, const std::string& path, uint64_t K) {
    std::vector<uint32_t> off(K + 1, 0);
    for (uint64_t i = 0; i < K; i++) off[i + 1] = off[i] + px.cnt[i];
    const uint64_t E = off[K];
    FILE* f = xopen(path + ".tmp", "wb");
    wr(f, &px.pad, 4); wr(f, &px.N, 4);
    uint64_t Kw = K; wr(f, &Kw, 8);
    wr(f, off.data(), (K + 1) * 4);
    std::vector<uint64_t> buf; buf.reserve(1u << 20);
    for (uint64_t i = 0; i < K; i++) {
        const uint64_t* row = px.w.data() + i * (size_t)px.N * 2;
        for (uint32_t j = 0; j < px.cnt[i]; j++) { buf.push_back(row[j * 2]); buf.push_back(row[j * 2 + 1]); }
        if (buf.size() >= (1u << 20)) { wr(f, buf.data(), buf.size() * 8); buf.clear(); }
    }
    if (!buf.empty()) wr(f, buf.data(), buf.size() * 8);
    atomic_close(f, path + ".tmp", path);
    uint64_t full = 0;
    for (uint64_t i = 0; i < K; i++) full += (px.cnt[i] >= px.N);
    std::fprintf(stderr, "[index] records skipped (no room in slot): %llu\n",
                 (unsigned long long)POST_SKIP.load());
    std::fprintf(stderr, "[index] %llu records, %.3f M nodes %.1f%% full (N=%u), %.2f GB -> %s\n",
                 (unsigned long long)E, K / 1e6, K ? 100.0 * (double)full / (double)K : 0.0, px.N,
                 (16.0 * (double)E + 4.0 * (double)(K + 1) + 16.0) / 1073741824.0, path.c_str());
}

struct FillCtx {
    const SigIndex* idx = nullptr;
    const std::vector<cov_t>* cov = nullptr;
    PostIdx* px = nullptr;
    int S = 0;
};

struct LocalB { EdgeKeySet e; };

void fill_read(std::string& s, int level, int margin, int col, uint64_t rno,
               const FillCtx& cx, LocalB& lo, uint64_t& ncores) {
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int nf, nr; const uint64_t *stf, *enf, *sgf, *strv, *enr, *sgr;
    extract2(s.data(), (int)len, level, &nf, &stf, &enf, &sgf, &nr, &strv, &enr, &sgr);
    const uint64_t mg = (uint64_t)margin;
    const int S = cx.S;
    const cov_t* cov = cx.cov->data();

    static thread_local std::vector<int> ci;
    static thread_local std::vector<uint32_t> cnode;
    auto chain = [&](bool rev, int n, const uint64_t* st, const uint64_t* en, const uint64_t* sig) {
        const bool wantpost = cx.px && cx.px->on() && !rev;
        if (wantpost) { ci.clear(); cnode.clear(); }
        uint32_t prev = UINT32_MAX; bool have_prev = false; uint64_t prev_sig = 0;
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            ncores++;
            const uint32_t node = cx.idx->find(sig[i]);

            if (have_prev && prev_sig != sig[i] && prev != UINT32_MAX && node != UINT32_MAX &&
                prev != node && cov[(size_t)prev * S + col] != 0 && cov[(size_t)node * S + col] != 0) {
                lo.e.insert(((uint64_t)prev << 32) | (uint64_t)node);
            }
            if (wantpost) { ci.push_back(i); cnode.push_back(node); }
            prev = node; prev_sig = sig[i]; have_prev = true;
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
    };
    chain(false, nf, stf, enf, sgf);
    chain(true, nr, strv, enr, sgr);
}

struct Feed {
    std::mutex qmtx;
    std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    std::atomic<int> readers_left{1};
    std::atomic<size_t> next{0};
    std::vector<uint64_t> reads_per;
    uint64_t max_bases = 0;
};

void read_shards(const std::vector<Shard>& shards, Feed& fd, size_t q_cap, size_t batch_n) {
    constexpr unsigned BLOCK = 1u << 20;
    std::vector<char> buf(BLOCK);
    Batch batch; batch.seqs.reserve(batch_n);
    std::string carry;
    std::vector<uint64_t> mine(fd.reads_per.size(), 0);
    auto push = [&]() {
        if (batch.seqs.empty()) return;
        std::unique_lock<std::mutex> lk(fd.qmtx);
        fd.cv_full.wait(lk, [&]{ return fd.q.size() < q_cap; });
        fd.q.push(std::move(batch)); fd.cv_empty.notify_one();
        batch.seqs.clear(); batch.seqs.reserve(batch_n);
    };
    size_t si; std::string svseq;
    while ((si = fd.next.fetch_add(1)) < shards.size()) {
        const int col = shards[si].col;
        if (!batch.seqs.empty() && batch.col != col) push();
        batch.col = col;
        if (batch.seqs.empty()) batch.first = mine[col];
        const bool svb = is_svb(shards[si].path);
        FILE* f = svb ? open_svb(shards[si].path) : open_fastq(shards[si].path);
        uint64_t bases = 0;
        if (svb) {
            while (svb_next(f, svseq)) {
                if (batch.seqs.empty()) batch.first = mine[col];
                bases += svseq.size();
                mine[col]++; batch.seqs.emplace_back(std::move(svseq));
                if (batch.seqs.size() >= batch_n) push();
                if (fd.max_bases && bases >= fd.max_bases) break;
            }
        } else {
            int line_no = 0; carry.clear(); int got; bool done = false;
            while (!done && (got = (int)std::fread(buf.data(), 1, BLOCK, f)) > 0) {
                const char* base = buf.data(); const char* end = base + got; const char* start = base;
                for (const char* p = base; p < end && !done; p++) if (*p == '\n') {
                    const size_t n = (size_t)(p - start);
                    if (line_no == 1) {
                        if (batch.seqs.empty()) batch.first = mine[col];
                        if (carry.empty()) batch.seqs.emplace_back(start, n);
                        else { carry.append(start, n); batch.seqs.push_back(std::move(carry)); carry.clear(); }
                        bases += batch.seqs.back().size();
                        mine[col]++;
                        if (batch.seqs.size() >= batch_n) push();
                        done = fd.max_bases && bases >= fd.max_bases;
                    } else if (!carry.empty()) carry.clear();
                    line_no = (line_no + 1) & 3; start = p + 1;
                }
                if (start < end && line_no == 1) carry.append(start, end - start);
            }
        }
        std::fclose(f);
        push();
    }
    push();
    { std::lock_guard<std::mutex> lk(fd.qmtx);
      for (size_t i = 0; i < mine.size(); i++) fd.reads_per[i] += mine[i]; }
    if (fd.readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(fd.qmtx); fd.cv_empty.notify_all(); }
}
bool feed_take(Feed& fd, Batch& out) {
    std::unique_lock<std::mutex> lk(fd.qmtx);
    fd.cv_empty.wait(lk, [&]{ return !fd.q.empty() || fd.readers_left.load() == 0; });
    if (fd.q.empty()) return false;
    out = std::move(fd.q.front()); fd.q.pop(); fd.cv_full.notify_one();
    return true;
}
void count_sample(const std::vector<Shard>& shards, int level, int margin, int threads, uint64_t max_bases,
                  ShTab* psn, std::mutex* nmtx,
                  std::atomic<uint64_t>& emissions, uint64_t& reads_here) {
    const int nw = std::max(1, threads - 1);
    Feed fd; fd.reads_per.assign(1, 0); fd.max_bases = max_bases;
    std::vector<Shard> mine = shards;
    for (Shard& s : mine) s.col = 0;

    auto worker = [&]() {
        uint64_t emis = 0;
        const char* fe = std::getenv("SNAPE_FLUSH");
        const size_t FLUSH = fe ? (size_t)std::atoll(fe) : (size_t)(1u << 20);
        SigCntMap lo; lo.reserve(1u << 20);
        std::vector<std::vector<std::pair<uint64_t, cov_t>>> bn(NSHARD);
        auto flush = [&]() {
            for (auto& v : bn) v.clear();
            for (const auto& kv : lo) bn[kv.first & (NSHARD - 1)].push_back(kv);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (bn[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(nmtx[sh]);
                for (const auto& kv : bn[sh]) addsat(psn[sh].m[kv.first], kv.second);
            }
            lo.clear();
        };
        Batch b;
        while (feed_take(fd, b)) {
            for (std::string& s : b.seqs) count_read(s, level, margin, lo, emis);
            if (lo.size() >= FLUSH) flush();
        }
        flush();
        emissions += emis;
    };

    { std::vector<std::thread> th;
      th.emplace_back([&]{ read_shards(mine, fd, 256, 64); });
      for (int i = 0; i < nw; i++) th.emplace_back(worker);
      for (auto& t : th) t.join(); }
    reads_here = fd.reads_per[0];
}

void merge_sample(ShTab* psn, NodeTab* gnode, int col, int S, int min_mult, int threads,
                  uint64_t& seen, uint64_t& kept, uint64_t& one, uint64_t& two) {
    const int nt = std::max(1, threads);
    std::atomic<uint64_t> nseen{0}, npass{0}, n1{0}, n2{0};
    std::atomic<int> next_sh{0};
    std::vector<std::thread> th;
    for (int i = 0; i < nt; i++) th.emplace_back([&]{
        for (int sh = next_sh.fetch_add(1); sh < NSHARD; sh = next_sh.fetch_add(1)) {
            uint64_t a = 0, b = 0, c1 = 0, c2 = 0;
            for (const auto& kv : psn[sh].m) {
                a++;
                c1 += (kv.second == 1); c2 += (kv.second == 2);
                if ((uint32_t)kv.second < (uint32_t)min_mult) continue;
                b++;
                node_slot(gnode[sh], kv.first, S)[col] = kv.second;
            }
            nseen += a; npass += b; n1 += c1; n2 += c2;
            SigCntMap().swap(psn[sh].m);
        }
    });
    for (auto& t : th) t.join();
    seen = nseen.load(); kept = npass.load(); one = n1.load(); two = n2.load();
}

void fill_pass(const std::vector<Shard>& shards, int level, int margin, int threads, int S, uint64_t max_bases,
               const SigIndex& idx, const std::vector<cov_t>& cov,
               PostIdx& px, EdgeKeySet* gedge, std::mutex* emtx,
               std::atomic<uint64_t>& emissions, std::vector<uint64_t>& reads_per) {
    const int nw = std::max(1, threads - 1);
    Feed fd; fd.reads_per.assign((size_t)S, 0); fd.max_bases = max_bases;
    FillCtx cx; cx.idx = &idx; cx.cov = &cov; cx.px = &px; cx.S = S;

    auto worker = [&]() {
        uint64_t emis = 0;
        const char* fe = std::getenv("SNAPE_FLUSH");
        const size_t FLUSH = fe ? (size_t)std::atoll(fe) : (size_t)(1u << 20);
        LocalB lo; lo.e.reserve(1u << 20);
        std::vector<std::vector<uint64_t>> be(NSHARD);
        auto flush = [&]() {
            for (auto& v : be) v.clear();
            for (uint64_t e : lo.e) be[(size_t)((e >> 32) ^ (e * 0x9e3779b97f4a7c15ULL)) & (NSHARD - 1)].push_back(e);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (be[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(emtx[sh]);
                for (uint64_t e : be[sh]) gedge[sh].insert(e);
            }
            lo.e.clear();
        };
        Batch b;
        while (feed_take(fd, b)) {
            for (size_t k = 0; k < b.seqs.size(); k++)
                fill_read(b.seqs[k], level, margin, b.col, b.first + k, cx, lo, emis);
            if (lo.e.size() >= FLUSH) flush();
        }
        flush();
        emissions += emis;
    };

    { std::vector<std::thread> th;
      th.emplace_back([&]{ read_shards(shards, fd, 256, 64); });
      for (int i = 0; i < nw; i++) th.emplace_back(worker);
      for (auto& t : th) t.join(); }
    for (int i = 0; i < S; i++) reads_per[i] = fd.reads_per[i];
}

void run_pool_two_pass(const std::vector<Shard>& shards, int level, int margin, int threads, uint64_t max_bases,
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
        distinct += seen; kept += pass;
        std::fprintf(stderr, "[graphA] sample %d/%d: %llu reads, %.3f M distinct -> %.3f M >= %d (%.1f%%)\n",
                     col + 1, S, (unsigned long long)rd, seen / 1e6, pass / 1e6, min_mult,
                     seen ? 100.0 * pass / seen : 0.0);
        std::fprintf(stderr, "[graphA]   distinct seen once %.3f M (%.1f%%), twice %.3f M (%.1f%%)\n",
                     one / 1e6, seen ? 100.0 * one / seen : 0.0,
                     two / 1e6, seen ? 100.0 * two / seen : 0.0);
    }

    {
        std::vector<std::vector<SigCount>> part(NSHARD);
        std::vector<std::vector<cov_t>> pcov(NSHARD);
        std::atomic<int> next_sh{0};
        std::vector<std::thread> th;
        for (int i = 0, nt = std::max(1, threads); i < nt; i++) th.emplace_back([&]{
            for (int sh = next_sh.fetch_add(1); sh < NSHARD; sh = next_sh.fetch_add(1)) {
                part[sh].reserve(gnode[sh].m.size());
                for (const auto& kv : gnode[sh].m) {
                    const cov_t* c = gnode[sh].a.data() + (size_t)kv.second * (size_t)S;
                    part[sh].push_back({kv.first, (uint32_t)(pcov[sh].size() / (size_t)S)});
                    pcov[sh].insert(pcov[sh].end(), c, c + S);
                }
                SigNodeMap().swap(gnode[sh].m);
                std::vector<cov_t>().swap(gnode[sh].a);
            }
        });
        for (auto& t : th) t.join();
        size_t K = 0;
        for (int sh = 0; sh < NSHARD; sh++) K += part[sh].size();
        out_nodes.clear(); out_nodes.reserve(K);
        std::vector<cov_t> tmpcov; tmpcov.reserve(K * (size_t)S);
        for (int sh = 0; sh < NSHARD; sh++) {
            const uint32_t base = (uint32_t)(tmpcov.size() / (size_t)S);
            for (SigCount c : part[sh]) { c.slot += base; out_nodes.push_back(c); }
            tmpcov.insert(tmpcov.end(), pcov[sh].begin(), pcov[sh].end());
            std::vector<SigCount>().swap(part[sh]);
            std::vector<cov_t>().swap(pcov[sh]);
        }
        psort(out_nodes, threads, [](const SigCount& a, const SigCount& b){ return a.sig < b.sig; });
        out_ncov.assign(K * (size_t)S, 0);
        for (size_t i = 0; i < K; i++) {
            const cov_t* c = tmpcov.data() + (size_t)out_nodes[i].slot * (size_t)S;
            std::copy(c, c + S, out_ncov.begin() + (std::ptrdiff_t)(i * (size_t)S));
            out_nodes[i].slot = (uint32_t)i;
        }
    }

    const size_t K = out_nodes.size();
    SigIndex idx; idx.build(out_nodes);
    std::fprintf(stderr, "[graphA] survivors=%.3f M, index %.2f GB -> pass B\n",
                 K / 1e6, (idx.bytes() + K * 8) / 1073741824.0);

    EdgeKeySet gedge[NSHARD];
    std::atomic<uint64_t> emis_b{0};
    std::vector<uint64_t> rp_b((size_t)S, 0);

    PostIdx px;
    if (index_n > 0) {
        px.N = (uint32_t)index_n; px.pad = (uint32_t)pad_cores;
        px.w.assign(K * (size_t)index_n * 2, 0);
        px.cnt.assign(K, 0);
        px.worst.assign(K, 0);
        px.mtx = smtx;
        std::fprintf(stderr, "[index] N=%d pad=%d, buffer %.2f GB\n",
                     index_n, pad_cores, px.bytes() / 1073741824.0);
    }
    fill_pass(shards, level, margin, threads, S, max_bases, idx, out_ncov, px,
              gedge, emtx, emis_b, rp_b);
    if (px.on()) {
        post_write(px, dir + "/index.bin", K);
        std::vector<uint64_t>().swap(px.w);
        std::vector<uint8_t>().swap(px.cnt);
        std::vector<uint64_t>().swap(px.worst);
    }

    {
        uint64_t ne = 0;
        for (int sh = 0; sh < NSHARD; sh++) ne += gedge[sh].size();
        out_edges.clear(); out_edges.reserve(ne);
        for (int sh = 0; sh < NSHARD; sh++) {
            for (uint64_t e : gedge[sh]) out_edges.push_back({(uint32_t)(e >> 32), (uint32_t)e});
            EdgeKeySet().swap(gedge[sh]);
        }
    }
    std::fprintf(stderr, "[graphA] emissions=%.3f M\n", emis_a.load() / 1e6);
    std::fprintf(stderr, "[graphB] emissions=%.3f M\n", emis_b.load() / 1e6);
    std::fprintf(stderr, "[graph1] per-sample distinct sigs=%.3f M -> passed mm=%d: %.3f M (%.2f%%)\n",
                 distinct / 1e6, min_mult, kept / 1e6, distinct ? 100.0 * kept / distinct : 0.0);
    std::fprintf(stderr, "[graph1] merged nodes=%.3f M | edges=%.3f M\n",
                 K / 1e6, out_edges.size() / 1e6);
}

}

void run_graph_reads(const std::vector<Sample>& samples, const std::string& dir,
                     int level, int margin, int min_mult, int threads, uint64_t max_bases, int index_n, int pad_cores) {
    mkdir(dir.c_str(), 0755);
    const int S = (int)samples.size();
    std::vector<Shard> shards_all;
    for (int i = 0; i < S; i++)
        for (const auto& p : collect_shards(samples[i].path)) shards_all.push_back({p, i});
    if (max_bases && shards_all.size() != (size_t)S) { std::fprintf(stderr, "[FATAL] --max-bases needs one file per sample\n"); std::abort(); }
    std::fprintf(stderr, "[graph] %d samples (%zu shards), threads=%d min_mult=%d max_bases=%llu\n",
                 S, shards_all.size(), threads, min_mult, (unsigned long long)max_bases);

    std::vector<SigCount> nodes; std::vector<CEdge> ce;
    std::vector<cov_t> ncov; std::vector<uint64_t> nreads;
    run_pool_two_pass(shards_all, level, margin, threads, max_bases, min_mult, S, index_n, pad_cores, dir,
                      nodes, ncov, ce, nreads);
    const uint64_t K = nodes.size();
    std::fprintf(stderr, "[graph]   survivors=%lu\n", (unsigned long)K);

    psort(ce, threads, [](const CEdge& x, const CEdge& y){ return x.src != y.src ? x.src < y.src : x.dst < y.dst; });

    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen(gp + ".tmp", "wb");
    const uint32_t Sw = (uint32_t)S;
    wr(fg, &Sw, 4);
    uint64_t Kw = K; wr(fg, &Kw, 8);
    { std::vector<uint64_t> sig(K);
      for (uint64_t i = 0; i < K; i++) sig[i] = nodes[i].sig;
      wr(fg, sig.data(), K * 8); }
    std::vector<cov_t> plane(K);
    for (int j = 0; j < S; j++) {
        for (uint64_t i = 0; i < K; i++) plane[i] = ncov[(size_t)nodes[i].slot * (size_t)S + j];
        wr(fg, plane.data(), K * sizeof(cov_t));
    }

    std::vector<SigCount>().swap(nodes); std::vector<cov_t>().swap(ncov);
    const uint64_t E = ce.size(); wr(fg, &E, 8);
    { std::vector<uint32_t> col(E);
      for (uint64_t i = 0; i < E; i++) col[i] = ce[i].src;
      wr(fg, col.data(), E * 4);
      for (uint64_t i = 0; i < E; i++) col[i] = ce[i].dst;
      wr(fg, col.data(), E * 4); }
    atomic_close(fg, gp + ".tmp", gp);
    for (int i = 0; i < S; i++)
        std::fprintf(stderr, "[graph]   %-8s %-6s %lu reads\n", samples[i].name.c_str(),
                     samples[i].normal ? "normal" : "tumor", (unsigned long)nreads[i]);
    std::fprintf(stderr, "[graph] nodes=%lu edges=%lu S=%d -> %s/graph.bin\n",
                 (unsigned long)K, (unsigned long)E, S, dir.c_str());
}
