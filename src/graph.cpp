#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "graph.h"
#include "samples.h"

#include <sys/stat.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <sched.h>
#include <dirent.h>
#include <algorithm>
#include <atomic>
#include <chrono>
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
void rd(FILE* f, void* p, size_t n) { if (n && std::fread(p, 1, n, f) != n) { std::fprintf(stderr, "[FATAL] truncated checkpoint\n"); std::abort(); } }
void atomic_close(FILE* f, const std::string& tmp, const std::string& dst) { std::fclose(f); if (std::rename(tmp.c_str(), dst.c_str())) { std::fprintf(stderr, "[FATAL] rename %s\n", dst.c_str()); std::abort(); } }

static inline uint64_t thr_ns() { struct timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec; }
static double self_cpu_s()  { struct rusage r; getrusage(RUSAGE_SELF, &r);     return r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1e6 + r.ru_stime.tv_sec + r.ru_stime.tv_usec / 1e6; }
static double peak_gb()     { struct rusage r; getrusage(RUSAGE_SELF, &r);     return r.ru_maxrss / 1048576.0; }
static double rss_gb() {
    FILE* f = std::fopen("/proc/self/statm", "r"); if (!f) return 0.0;
    long total = 0, res = 0;
    if (std::fscanf(f, "%ld %ld", &total, &res) != 2) res = 0;
    std::fclose(f);
    return (double)res * 4096.0 / 1073741824.0;
}

namespace {
using SEdge = std::pair<uint64_t, uint64_t>;

struct Tic { std::chrono::steady_clock::time_point t0, prev; const char* tag;
    explicit Tic(const char* tg) : t0(std::chrono::steady_clock::now()), prev(t0), tag(tg) {}
    void mark(const char* lbl) {
        const auto now = std::chrono::steady_clock::now();
        const double dprev = std::chrono::duration<double>(now - prev).count();
        const double dtot = std::chrono::duration<double>(now - t0).count();
        std::fprintf(stderr, "[%s/%s] +%.2fs  (total %.2fs)\n", tag, lbl, dprev, dtot);
        prev = now;
    }
};

struct EdgeHash { uint64_t operator()(const SEdge& e) const { return e.first * 0x9e3779b97f4a7c15ULL + e.second; } };
using EdgeSet = ankerl::unordered_dense::set<SEdge, EdgeHash>;

struct Prof {
    uint64_t ins_ns = 0;
    uint64_t seq_ns = 0;
    uint64_t seq_smp = 0;
    uint64_t seq_n = 0;
    uint64_t seq_bp = 0;
    uint64_t seq_rev = 0;
    uint64_t edge_n = 0;
    uint64_t wait_ns = 0;
    uint64_t read_ns = 0;
    uint64_t local_hw = 0;
};
Prof PROF;
std::mutex PROF_MTX;

uint64_t CLK_NS = 0;

constexpr uint64_t SEQ_SAMPLE = 1024;

void prof_merge(const Prof& p) {
    std::lock_guard<std::mutex> lk(PROF_MTX);
    PROF.ins_ns += p.ins_ns; PROF.seq_ns += p.seq_ns; PROF.seq_smp += p.seq_smp;
    PROF.seq_n  += p.seq_n;  PROF.seq_bp += p.seq_bp; PROF.seq_rev += p.seq_rev;
    PROF.edge_n += p.edge_n; PROF.wait_ns += p.wait_ns; PROF.read_ns += p.read_ns;
    PROF.local_hw += p.local_hw;
}

uint64_t clock_cost_ns() {
    uint64_t best = ~0ull;
    for (int i = 0; i < 1000; i++) {
        const uint64_t a = thr_ns(), b = thr_ns();
        best = std::min(best, b - a);
    }
    return best;
}

static inline void bump(cov_t& c) { c += (c != COV_MAX); }

struct SeqRef { uint64_t off = 0; uint32_t len = 0; };

static inline uint32_t b2(char c) {
    const uint32_t v = ((uint32_t)(uint8_t)c >> 1) & 3;
    return v ^ (v >> 1);
}

constexpr uint32_t ACGT_WORD = 0x54474341u;

static inline uint32_t acgt_bad(uint8_t c, uint32_t v) {
    return ((ACGT_WORD >> (v * 8)) & 0xff) ^ (uint32_t)(c & 0xdfu);
}

static SeqRef seq_put(std::vector<uint8_t>& a, const uint8_t* v, uint32_t n) {
    SeqRef r; r.off = a.size(); r.len = n;
    a.resize(r.off + (n + 3) / 4);
    uint8_t* o = a.data() + r.off;
    uint32_t i = 0;
    for (; i + 4 <= n; i += 4)
        o[i >> 2] = (uint8_t)((v[i] << 6) | (v[i + 1] << 4) | (v[i + 2] << 2) | v[i + 3]);
    if (i < n) {
        uint32_t last = 0;
        for (uint32_t j = i; j < n; j++) last |= (uint32_t)v[j] << (6 - 2 * (j & 3));
        o[i >> 2] = (uint8_t)last;
    }
    return r;
}

static SeqRef seq_put_raw(std::vector<uint8_t>& a, const uint8_t* v, uint32_t n) {
    SeqRef r; r.off = a.size(); r.len = n | SEQ_RAW;
    for (uint32_t i = 0; i < n; i++) a.push_back((uint8_t)"ACGT"[v[i]]);
    return r;
}

static SeqRef seq_copy(std::vector<uint8_t>& dst, const uint8_t* src, SeqRef r) {
    SeqRef o; o.off = dst.size(); o.len = r.len;
    const uint64_t nb = seq_bytes(r.len);
    dst.insert(dst.end(), src + r.off, src + r.off + nb);
    return o;
}
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
struct Batch { std::vector<std::string> seqs; int col = 0; };

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
    char hdr[32];
    if (std::fread(hdr, 1, 32, f) != 32 || std::memcmp(hdr, "SVB1", 4) != 0) {
        std::fprintf(stderr, "[FATAL] bad .svb header %s\n", path.c_str()); std::abort(); }
    uint32_t ver; std::memcpy(&ver, hdr + 4, 4);
    if (ver != 1) { std::fprintf(stderr, "[FATAL] .svb version %u != 1 (%s)\n", ver, path.c_str()); std::abort(); }
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


void report_bytes_a(int col, const ShTab* psn, const NodeTab* gnode) {
    size_t pm = 0, gm = 0, ga = 0;
    for (int sh = 0; sh < NSHARD; sh++) {
        pm += psn[sh].m.values().capacity() * sizeof(std::pair<uint64_t, cov_t>) + psn[sh].m.bucket_count() * 8;
        gm += gnode[sh].m.values().capacity() * sizeof(std::pair<uint64_t, uint32_t>) + gnode[sh].m.bucket_count() * 8;
        ga += gnode[sh].a.capacity() * sizeof(cov_t);
    }
    const double G = 1073741824.0;
    std::fprintf(stderr, "[bytes] A sample %d, before merge: sample-map %.2f || graph-map %.2f | "
                 "graph-cov %.2f || sum %.2f GB, RSS %.2f GB\n",
                 col, pm / G, gm / G, ga / G, (pm + gm + ga) / G, rss_gb());
}

void report_bytes_b(const EdgeKeySet* ge, const std::vector<uint8_t>* sq, size_t idx_bytes,
                    size_t cov_bytes, size_t sref_bytes) {
    size_t e = 0, s = 0;
    for (int sh = 0; sh < NSHARD; sh++) {
        e += ge[sh].values().capacity() * 8 + ge[sh].bucket_count() * 8;
        s += sq[sh].capacity();
    }
    const double G = 1073741824.0;
    std::fprintf(stderr, "[bytes] B: survivor-index %.2f | graph-cov %.2f | seqref %.2f | "
                 "edges %.2f | seq %.2f || sum %.2f GB, RSS %.2f GB\n",
                 idx_bytes / G, cov_bytes / G, sref_bytes / G, e / G, s / G,
                 (idx_bytes + cov_bytes + sref_bytes + e + s) / G, rss_gb());
}

void rss_poll(std::atomic<bool>* stop, std::string path, double* peak, double* peak_t) {
    FILE* f = xopen(path + ".tmp", "w");
    std::fprintf(f, "t_s\trss_gb\n");
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    while (!stop->load()) {
        const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const double r = rss_gb();
        std::fprintf(f, "%.2f\t%.3f\n", t, r);
        if (r > *peak) { *peak = r; *peak_t = t; }
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    atomic_close(f, path + ".tmp", path);
}

void count_read(std::string& s, int level, int margin, SigCntMap& lo,
                uint64_t& cpu_ext, uint64_t& ncores, uint64_t& ov_s, uint64_t& ov_n, Prof& pf) {
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int nf, nr; const uint64_t *stf, *enf, *sgf, *strv, *enr, *sgr;
    const uint64_t te = thr_ns();
    extract2(s.data(), (int)len, level, &nf, &stf, &enf, &sgf, &nr, &strv, &enr, &sgr);
    cpu_ext += thr_ns() - te;
    const uint64_t mg = (uint64_t)margin;
    auto chain = [&](int n, const uint64_t* st, const uint64_t* en, const uint64_t* sig) {
        for (int i = 0; i + 1 < n; i++) if (en[i] > st[i + 1]) { ov_s += en[i] - st[i + 1]; ov_n++; }
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            ncores++;
            bump(lo[sig[i]]);
        }
    };
    const uint64_t ti = thr_ns();
    chain(nf, stf, enf, sgf);
    chain(nr, strv, enr, sgr);
    pf.ins_ns += thr_ns() - ti;
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

struct FillCtx {
    const SigIndex* idx = nullptr;
    const std::vector<cov_t>* cov = nullptr;
    std::vector<std::atomic<uint8_t>>* have = nullptr;
    std::vector<SeqRef>* sref = nullptr;
    int S = 0;
};

struct LocalB {
    EdgeKeySet e;
    std::vector<uint8_t> a;
    std::vector<std::pair<uint32_t, SeqRef>> put;
};

void fill_read(std::string& s, int level, int margin, int col, const FillCtx& cx, LocalB& lo,
               uint64_t& cpu_ext, uint64_t& ncores, Prof& pf) {
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int nf, nr; const uint64_t *stf, *enf, *sgf, *strv, *enr, *sgr;
    const uint64_t te = thr_ns();
    extract2(s.data(), (int)len, level, &nf, &stf, &enf, &sgf, &nr, &strv, &enr, &sgr);
    cpu_ext += thr_ns() - te;
    const uint64_t mg = (uint64_t)margin;
    const int S = cx.S;
    const cov_t* cov = cx.cov->data();

    static thread_local std::vector<uint8_t> cd;
    auto chain = [&](bool rev, int n, const uint64_t* st, const uint64_t* en, const uint64_t* sig) {
        uint32_t prev = UINT32_MAX; bool have_prev = false; uint64_t prev_sig = 0;
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            ncores++;
            const uint32_t node = cx.idx->find(sig[i]);

            if (node != UINT32_MAX && en[i] > st[i] && en[i] <= len &&
                (*cx.have)[node].load(std::memory_order_relaxed) == 0 &&
                (*cx.have)[node].exchange(1, std::memory_order_relaxed) == 0) {
                const bool smp = (pf.seq_n % SEQ_SAMPLE) == 0;
                const uint64_t ts = smp ? thr_ns() : 0;
                const uint32_t w = (uint32_t)(en[i] - st[i]);
                const char* src = s.data() + (rev ? len - en[i] : st[i]);
                const uint32_t flip = rev ? 3u : 0u;
                cd.resize(w);
                uint32_t m = 0, bad = 0, pv = 4;
                for (uint32_t x = 0; x < w; x++) {
                    const uint8_t c = (uint8_t)src[rev ? w - 1 - x : x];
                    const uint32_t v0 = b2((char)c);
                    bad |= acgt_bad(c, v0);
                    const uint32_t v = v0 ^ flip;
                    cd[m] = (uint8_t)v;
                    m += (v != pv);
                    pv = v;
                }
                lo.put.emplace_back(node, bad ? seq_put_raw(lo.a, cd.data(), m)
                                              : seq_put(lo.a, cd.data(), m));
                if (smp) { pf.seq_ns += thr_ns() - ts; pf.seq_smp++; }
                pf.seq_n++; pf.seq_bp += m; pf.seq_rev += rev;
            }

            if (have_prev && prev_sig != sig[i] && prev != UINT32_MAX && node != UINT32_MAX &&
                prev != node && cov[(size_t)prev * S + col] != 0 && cov[(size_t)node * S + col] != 0) {
                lo.e.insert(((uint64_t)prev << 32) | (uint64_t)node);
                pf.edge_n++;
            }
            prev = node; prev_sig = sig[i]; have_prev = true;
        }
    };
    const uint64_t ti = thr_ns();
    chain(false, nf, stf, enf, sgf);
    chain(true, nr, strv, enr, sgr);
    pf.ins_ns += thr_ns() - ti;
}

struct Feed {
    std::mutex qmtx;
    std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    std::atomic<int> readers_left{1};
    std::atomic<size_t> next{0};
    std::vector<uint64_t> reads_per;
};

void read_shards(const std::vector<Shard>& shards, Feed& fd, size_t q_cap, size_t batch_n) {
    const uint64_t r0 = thr_ns();
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
        const bool svb = is_svb(shards[si].path);
        FILE* f = svb ? open_svb(shards[si].path) : open_fastq(shards[si].path);
        if (svb) {
            while (svb_next(f, svseq)) {
                mine[col]++; batch.seqs.emplace_back(std::move(svseq));
                if (batch.seqs.size() >= batch_n) push();
            }
        } else {
            int line_no = 0; carry.clear(); int got;
            while ((got = (int)std::fread(buf.data(), 1, BLOCK, f)) > 0) {
                const char* base = buf.data(); const char* end = base + got; const char* start = base;
                for (const char* p = base; p < end; p++) if (*p == '\n') {
                    const size_t n = (size_t)(p - start);
                    if (line_no == 1) {
                        if (carry.empty()) batch.seqs.emplace_back(start, n);
                        else { carry.append(start, n); batch.seqs.push_back(std::move(carry)); carry.clear(); }
                        mine[col]++;
                        if (batch.seqs.size() >= batch_n) push();
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
    Prof rp; rp.read_ns = thr_ns() - r0; prof_merge(rp);
}

bool feed_take(Feed& fd, Batch& out, Prof& pf) {
    std::unique_lock<std::mutex> lk(fd.qmtx);
    const std::chrono::steady_clock::time_point w0 = std::chrono::steady_clock::now();
    fd.cv_empty.wait(lk, [&]{ return !fd.q.empty() || fd.readers_left.load() == 0; });
    pf.wait_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(
                      std::chrono::steady_clock::now() - w0).count();
    if (fd.q.empty()) return false;
    out = std::move(fd.q.front()); fd.q.pop(); fd.cv_full.notify_one();
    return true;
}

void count_sample(const std::vector<Shard>& shards, int level, int margin, int threads,
                  ShTab* psn, std::mutex* nmtx,
                  std::atomic<uint64_t>& emissions, std::atomic<uint64_t>& extract_ns,
                  std::atomic<uint64_t>& ov_sum, std::atomic<uint64_t>& ov_cnt,
                  std::atomic<uint64_t>& flush_ns, uint64_t& reads_here) {
    const int nw = std::max(1, threads - 1);
    Feed fd; fd.reads_per.assign(1, 0);
    std::vector<Shard> mine = shards;
    for (Shard& s : mine) s.col = 0;

    auto worker = [&]() {
        uint64_t ext = 0, emis = 0, ovs = 0, ovn = 0;
        const char* fe = std::getenv("SNAPE_FLUSH");
        const size_t FLUSH = fe ? (size_t)std::atoll(fe) : (size_t)(1u << 20);
        Prof pf;
        SigCntMap lo; lo.reserve(1u << 20);
        std::vector<std::vector<std::pair<uint64_t, cov_t>>> bn(NSHARD);
        auto flush = [&]() {
            const uint64_t tf = thr_ns();
            for (auto& v : bn) v.clear();
            for (const auto& kv : lo) bn[kv.first & (NSHARD - 1)].push_back(kv);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (bn[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(nmtx[sh]);
                for (const auto& kv : bn[sh]) addsat(psn[sh].m[kv.first], kv.second);
            }
            lo.clear();
            flush_ns += thr_ns() - tf;
        };
        Batch b;
        while (feed_take(fd, b, pf)) {
            for (std::string& s : b.seqs) count_read(s, level, margin, lo, ext, emis, ovs, ovn, pf);
            pf.local_hw = std::max(pf.local_hw,
                lo.values().capacity() * sizeof(std::pair<uint64_t, cov_t>) + lo.bucket_count() * 8);
            if (lo.size() >= FLUSH) flush();
        }
        flush();
        emissions += emis; extract_ns += ext; ov_sum += ovs; ov_cnt += ovn;
        prof_merge(pf);
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

void fill_pass(const std::vector<Shard>& shards, int level, int margin, int threads, int S,
               const SigIndex& idx, const std::vector<cov_t>& cov,
               std::vector<std::atomic<uint8_t>>& have, std::vector<SeqRef>& sref,
               EdgeKeySet* gedge, std::vector<uint8_t>* gseq,
               std::mutex* emtx, std::mutex* smtx,
               std::atomic<uint64_t>& emissions, std::atomic<uint64_t>& extract_ns,
               std::atomic<uint64_t>& flush_ns, std::vector<uint64_t>& reads_per) {
    const int nw = std::max(1, threads - 1);
    Feed fd; fd.reads_per.assign((size_t)S, 0);
    FillCtx cx; cx.idx = &idx; cx.cov = &cov; cx.have = &have; cx.sref = &sref; cx.S = S;

    auto worker = [&]() {
        uint64_t ext = 0, emis = 0;
        const char* fe = std::getenv("SNAPE_FLUSH");
        const size_t FLUSH = fe ? (size_t)std::atoll(fe) : (size_t)(1u << 20);
        Prof pf;
        LocalB lo; lo.e.reserve(1u << 20);
        std::vector<std::vector<uint64_t>> be(NSHARD);
        std::vector<std::vector<std::pair<uint32_t, SeqRef>>> bs(NSHARD);
        auto flush = [&]() {
            const uint64_t tf = thr_ns();
            for (auto& v : be) v.clear();
            for (uint64_t e : lo.e) be[(size_t)((e >> 32) ^ (e * 0x9e3779b97f4a7c15ULL)) & (NSHARD - 1)].push_back(e);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (be[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(emtx[sh]);
                for (uint64_t e : be[sh]) gedge[sh].insert(e);
            }
            lo.e.clear();
            for (auto& v : bs) v.clear();
            for (const auto& p : lo.put) bs[idx.v[p.first].sig & (NSHARD - 1)].push_back(p);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (bs[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(smtx[sh]);
                for (const auto& p : bs[sh]) sref[p.first] = seq_copy(gseq[sh], lo.a.data(), p.second);
            }
            lo.put.clear(); lo.a.clear();
            flush_ns += thr_ns() - tf;
        };
        Batch b;
        while (feed_take(fd, b, pf)) {
            for (std::string& s : b.seqs) fill_read(s, level, margin, b.col, cx, lo, ext, emis, pf);
            const size_t lb = lo.e.values().capacity() * 8 + lo.e.bucket_count() * 8 + lo.a.capacity();
            pf.local_hw = std::max(pf.local_hw, lb);
            if (lo.e.size() >= FLUSH || lo.a.size() >= (size_t)(64u << 20)) flush();
        }
        flush();
        emissions += emis; extract_ns += ext;
        prof_merge(pf);
    };

    { std::vector<std::thread> th;
      th.emplace_back([&]{ read_shards(shards, fd, 256, 64); });
      for (int i = 0; i < nw; i++) th.emplace_back(worker);
      for (auto& t : th) t.join(); }
    for (int i = 0; i < S; i++) reads_per[i] = fd.reads_per[i];
}

void run_pool_two_pass(const std::vector<Shard>& shards, int level, int margin, int threads,
                       int min_mult, int S, std::vector<SigCount>& out_nodes,
                       std::vector<cov_t>& out_ncov, std::vector<CEdge>& out_edges,
                       std::vector<SeqRef>& out_sref, std::vector<uint8_t>& out_seq,
                       uint64_t& ov_s, uint64_t& ov_c, std::vector<uint64_t>& reads_per) {
    NodeTab gnode[NSHARD];
    ShTab psn[NSHARD];
    std::mutex nmtx[NSHARD], emtx[NSHARD], smtx[NSHARD];
    std::atomic<uint64_t> emis_a{0}, ext_a{0}, ov_sum{0}, ov_cnt{0}, flush_a{0};
    reads_per.assign((size_t)S, 0);
    uint64_t distinct = 0, kept = 0;

    for (int col = 0; col < S; col++) {
        std::vector<Shard> mine;
        for (const Shard& s : shards) if (s.col == col) mine.push_back(s);
        uint64_t rd = 0, seen = 0, pass = 0, one = 0, two = 0;
        count_sample(mine, level, margin, threads, psn, nmtx, emis_a, ext_a, ov_sum, ov_cnt, flush_a, rd);
        report_bytes_a(col, psn, gnode);
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

    std::vector<std::atomic<uint8_t>> have(K);
    for (size_t i = 0; i < K; i++) have[i].store(0, std::memory_order_relaxed);
    out_sref.assign(K, SeqRef{});
    EdgeKeySet gedge[NSHARD];
    std::vector<uint8_t> gseq[NSHARD];
    std::atomic<uint64_t> emis_b{0}, ext_b{0}, flush_b{0};
    std::vector<uint64_t> rp_b((size_t)S, 0);

    fill_pass(shards, level, margin, threads, S, idx, out_ncov, have, out_sref,
              gedge, gseq, emtx, smtx, emis_b, ext_b, flush_b, rp_b);
    report_bytes_b(gedge, gseq, idx.bytes() + K * 8, out_ncov.capacity(), out_sref.capacity() * sizeof(SeqRef));
    std::vector<std::atomic<uint8_t>>().swap(have);

    {
        uint64_t ne = 0;
        for (int sh = 0; sh < NSHARD; sh++) ne += gedge[sh].size();
        out_edges.clear(); out_edges.reserve(ne);
        for (int sh = 0; sh < NSHARD; sh++) {
            for (uint64_t e : gedge[sh]) out_edges.push_back({(uint32_t)(e >> 32), (uint32_t)e});
            EdgeKeySet().swap(gedge[sh]);
        }
    }
    {
        uint64_t base[NSHARD], tot = 0;
        for (int sh = 0; sh < NSHARD; sh++) { base[sh] = tot; tot += gseq[sh].size(); }
        out_seq.clear(); out_seq.reserve(tot);
        for (int sh = 0; sh < NSHARD; sh++) {
            out_seq.insert(out_seq.end(), gseq[sh].begin(), gseq[sh].end());
            std::vector<uint8_t>().swap(gseq[sh]);
        }
        for (size_t i = 0; i < K; i++)
            if (out_sref[i].len != 0) out_sref[i].off += base[out_nodes[i].sig & (NSHARD - 1)];
    }

    ov_s = ov_sum.load(); ov_c = ov_cnt.load();
    uint64_t stored = 0;
    for (size_t i = 0; i < K; i++) stored += (out_sref[i].len != 0);
    std::fprintf(stderr, "[graphA] emissions=%.3f M, extract %.1f CPU-s, flush %.1f CPU-s\n",
                 emis_a.load() / 1e6, ext_a.load() / 1e9, flush_a.load() / 1e9);
    std::fprintf(stderr, "[graphB] emissions=%.3f M, extract %.1f CPU-s, flush %.1f CPU-s\n",
                 emis_b.load() / 1e6, ext_b.load() / 1e9, flush_b.load() / 1e9);
    std::fprintf(stderr, "[graph1] per-sample distinct sigs=%.3f M -> passed mm=%d: %.3f M (%.2f%%)\n",
                 distinct / 1e6, min_mult, kept / 1e6, distinct ? 100.0 * kept / distinct : 0.0);
    std::fprintf(stderr, "[graph1] merged nodes=%.3f M | edges=%.3f M | seq stored for %.3f M nodes (%.1f%%)\n",
                 K / 1e6, out_edges.size() / 1e6, stored / 1e6, K ? 100.0 * stored / K : 0.0);
    std::fprintf(stderr, "[graph1] core sequence=%.3f GiB, ardisik cekirdek ortusme %.1f bp\n",
                 out_seq.size() / 1073741824.0, ov_c ? (double)ov_s / (double)ov_c : 0.0);
}

}

void run_graph_reads(const std::vector<Sample>& samples, const std::string& dir,
                     int level, int margin, int min_mult, int threads) {
    mkdir(dir.c_str(), 0755);
    const int S = (int)samples.size();
    std::vector<Shard> shards_all;
    for (int i = 0; i < S; i++)
        for (const auto& p : collect_shards(samples[i].path)) shards_all.push_back({p, i});
    std::fprintf(stderr, "[graph] %d samples (%zu shards), threads=%d min_mult=%d\n",
                 S, shards_all.size(), threads, min_mult);
    Tic tic("graph");
    CLK_NS = clock_cost_ns();
    std::atomic<bool> rss_stop{false};
    double rss_peak = 0, rss_peak_t = 0;
    std::thread rss_th(rss_poll, &rss_stop, dir + "/rss_trace.tsv", &rss_peak, &rss_peak_t);

    struct Snap { double wall, self, rss, peak; };
    const auto torig = std::chrono::steady_clock::now();
    auto snap = [&]() -> Snap { return { std::chrono::duration<double>(std::chrono::steady_clock::now() - torig).count(), self_cpu_s(), rss_gb(), peak_gb() }; };
    const Snap P0 = snap();

    std::vector<SigCount> nodes; std::vector<CEdge> ce;
    std::vector<cov_t> ncov; std::vector<uint64_t> nreads;
    std::vector<SeqRef> sref; std::vector<uint8_t> sbuf; uint64_t ov_s = 0, ov_c = 0;
    run_pool_two_pass(shards_all, level, margin, threads, min_mult, S, nodes, ncov, ce, sref, sbuf,
                    ov_s, ov_c, nreads);
    const uint64_t K = nodes.size();
    const Snap P1 = snap();
    tic.mark("pool_both");
    std::fprintf(stderr, "[graph]   survivors=%lu\n", (unsigned long)K);

    const Snap P2 = snap();

    psort(ce, threads, [](const CEdge& x, const CEdge& y){ return x.src != y.src ? x.src < y.src : x.dst < y.dst; });
    const Snap P4 = snap();
    tic.mark("sort_edges");

    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen(gp + ".tmp", "wb");
    const uint32_t ver = 5, Sw = (uint32_t)S, resv = 0;
    wr(fg, "SNPG", 4); wr(fg, &ver, 4); wr(fg, &Sw, 4); wr(fg, &resv, 4);
    wr(fg, &ov_s, 8); wr(fg, &ov_c, 8);
    uint64_t Kw = K; wr(fg, &Kw, 8);
    { std::vector<uint64_t> sig(K);
      for (uint64_t i = 0; i < K; i++) sig[i] = nodes[i].sig;
      wr(fg, sig.data(), K * 8); }
    std::vector<cov_t> plane(K);
    for (int j = 0; j < S; j++) {
        for (uint64_t i = 0; i < K; i++) plane[i] = ncov[(size_t)nodes[i].slot * (size_t)S + j];
        wr(fg, plane.data(), K * sizeof(cov_t));
    }

    { std::vector<uint32_t> sl(K);
      uint64_t tot = 0;
      for (uint64_t i = 0; i < K; i++) { sl[i] = sref[nodes[i].slot].len; tot += seq_bytes(sl[i]); }
      wr(fg, &tot, 8);
      wr(fg, sl.data(), K * 4);
      std::vector<uint8_t> buf; buf.reserve(1u << 22);
      for (uint64_t i = 0; i < K; i++) {
          const SeqRef& r = sref[nodes[i].slot];
          const uint64_t nb = seq_bytes(r.len);
          buf.insert(buf.end(), sbuf.data() + r.off, sbuf.data() + r.off + nb);
          if (buf.size() >= (1u << 22)) { wr(fg, buf.data(), buf.size()); buf.clear(); }
      }
      if (!buf.empty()) wr(fg, buf.data(), buf.size()); }
    std::vector<SigCount>().swap(nodes); std::vector<cov_t>().swap(ncov);
    std::vector<SeqRef>().swap(sref); std::vector<uint8_t>().swap(sbuf);
    const uint64_t E = ce.size(); wr(fg, &E, 8);
    { std::vector<uint32_t> col(E);
      for (uint64_t i = 0; i < E; i++) col[i] = ce[i].src;
      wr(fg, col.data(), E * 4);
      for (uint64_t i = 0; i < E; i++) col[i] = ce[i].dst;
      wr(fg, col.data(), E * 4); }
    atomic_close(fg, gp + ".tmp", gp);
    const Snap P5 = snap();
    tic.mark("write");
    for (int i = 0; i < S; i++)
        std::fprintf(stderr, "[graph]   %-8s %-6s %lu reads\n", samples[i].name.c_str(),
                     samples[i].normal ? "normal" : "tumor", (unsigned long)nreads[i]);
    std::fprintf(stderr, "[graph] nodes=%lu edges=%lu S=%d -> %s/graph.bin\n",
                 (unsigned long)K, (unsigned long)E, S, dir.c_str());

    rss_stop.store(true); rss_th.join();
    cpu_set_t cs; CPU_ZERO(&cs);
    const double ncpu = (sched_getaffinity(0, sizeof(cs), &cs) == 0 && CPU_COUNT(&cs) > 0) ? (double)CPU_COUNT(&cs) : (double)sysconf(_SC_NPROCESSORS_ONLN);
    auto row = [&](const char* name, const Snap& lo, const Snap& hi) {
        const double wall = hi.wall - lo.wall, cpu = hi.self - lo.self;
        std::fprintf(stderr, "[eff] %-16s wall=%8.1fs  CPU=%9.1f  eff=%5.1f%%  RSS=%6.2fGB  peak=%6.2fGB\n",
                     name, wall, cpu, 100.0 * cpu / (ncpu * wall + 1e-9), hi.rss, hi.peak);
    };
    std::fprintf(stderr, "\n===== graph build CPU-seconds & execution efficiency (ncpu=%.0f) =====\n", ncpu);
    std::fprintf(stderr, "[eff] %-16s %9s  %9s  %6s  %9s %9s\n", "phase", "wall_s", "CPU_s", "eff", "RSS_end", "peak");
    row("count(pool_both)", P0, P1);
    row("sig arrays",       P1, P2);
    row("sort_edges",       P2, P4);
    row("write",            P4, P5);
    row("TOTAL",            P0, P5);
    std::fprintf(stderr, "=====================================================================\n");

    const double seq_net = (double)PROF.seq_ns - (double)(PROF.seq_smp * CLK_NS);
    const double seq_est = PROF.seq_smp ? seq_net * (double)PROF.seq_n / (double)PROF.seq_smp : 0.0;
    std::fprintf(stderr, "[prof] per-core bookkeeping  %9.1f CPU-s\n", PROF.ins_ns / 1e9);
    std::fprintf(stderr, "[prof]   of which packing    %9.1f CPU-s (estimated from %llu of %llu blocks, "
                 "%.2f Gbp packed, %.1f%% of them reverse-strand)\n",
                 seq_est / 1e9, (unsigned long long)PROF.seq_smp, (unsigned long long)PROF.seq_n,
                 PROF.seq_bp / 1e9, PROF.seq_n ? 100.0 * PROF.seq_rev / PROF.seq_n : 0.0);
    std::fprintf(stderr, "[prof] worker-local tables    %9.2f GB (sum of per-worker high water) | "
                 "clock pair %llu ns, subtracted\n",
                 PROF.local_hw / 1073741824.0, (unsigned long long)CLK_NS);
    std::fprintf(stderr, "[prof] edge insertions       %9.3f G\n", PROF.edge_n / 1e9);
    std::fprintf(stderr, "[prof] reader thread CPU     %9.1f s | workers blocked on empty queue %9.1f s (wall, summed)\n",
                 PROF.read_ns / 1e9, PROF.wait_ns / 1e9);
    std::fprintf(stderr, "[prof] RSS trace peak %.2f GB at t=%.1fs -> %s/rss_trace.tsv\n",
                 rss_peak, rss_peak_t, dir.c_str());
}
