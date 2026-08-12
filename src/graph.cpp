#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1   // expose sched_getaffinity under strict -std=c++17
#endif
#include "graph.h"
// extract is declared in snape.h (provided by src/extract.cpp over biolcp).

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
#include "unordered_dense.h"   // ankerl::unordered_dense: open-addressing flat hash map/set keyed by u64
                               // (replaces node-based std::unordered_map/set — kills per-entry malloc + pointer-chase)
#include <utility>
#include <vector>

// ---- fail-fast binary checkpoint IO (declared in snape.h) ----
FILE* xopen(const std::string& p, const char* m) { FILE* f = std::fopen(p.c_str(), m); if (!f) { std::fprintf(stderr, "[FATAL] fopen %s\n", p.c_str()); std::abort(); } return f; }
void wr(FILE* f, const void* p, size_t n) { if (n && std::fwrite(p, 1, n, f) != n) { std::fprintf(stderr, "[FATAL] fwrite\n"); std::abort(); } }
void rd(FILE* f, void* p, size_t n) { if (n && std::fread(p, 1, n, f) != n) { std::fprintf(stderr, "[FATAL] truncated checkpoint\n"); std::abort(); } }
void atomic_close(FILE* f, const std::string& tmp, const std::string& dst) { std::fclose(f); if (std::rename(tmp.c_str(), dst.c_str())) { std::fprintf(stderr, "[FATAL] rename %s\n", dst.c_str()); std::abort(); } }   // never leave a half-written checkpoint resume would trust

// ---- CPU-seconds accounting (for the per-component cost/efficiency breakdown) ----
// thr_ns(): on-CPU time of the CALLING thread only (advances only while scheduled, so
// blocked/idle time is excluded -> true CPU-seconds, contention-immune).
// self_cpu_s(): kernel-accounted user+sys CPU for all in-process threads. With the input read
// directly there are no child processes left to account for.
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
using SEdge = std::pair<uint64_t, uint64_t>;        // an adjacency, keyed by its two signatures

// ---- NO REFERENCE ----
// The graph is built from reads alone. There is no reference index, so no strand vote and no repeat
// merge: every read enters BOTH ways (its bytes as sequenced and their reverse complement) and every
// kept core enters as extract() produced it. A locus therefore has a mirror twin in the other strand
// family; the twins are two nodes with different sigs and no core in common, and pairing them is the
// scorer's job, not the graph's. Repeats are ordinary nodes, told apart by COVERAGE.

// Node identity = core_sig (MurmurHash64A over the core's homopolymer-compressed bytes), computed in
// extract(). That sig is the only thing stored per node; no sequence is kept anywhere.

// ---- stage timing (wall-clock; one Tic per stage entry) ----
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

// Per-edge (prev_sig, cur_sig) -> (normal_count, tumor_count). Hash mixes the sig pair; not tagged
// is_avalanching so unordered_dense re-mixes it (prev*golden+cur alone is not fully avalanched).
struct EdgeHash { uint64_t operator()(const SEdge& e) const { return e.first * 0x9e3779b97f4a7c15ULL + e.second; } };
using EdgeMap = ankerl::unordered_dense::map<SEdge, std::pair<uint32_t, uint32_t>, EdgeHash>;

static bool is_gz(const std::string& p) { return p.size() >= 3 && p.compare(p.size() - 3, 3, ".gz") == 0; }

// A sample input is a directory of FASTQ/.svb shards (read in parallel) or a single file.
// Compressed input is rejected, not decompressed: `.fastq.gz` still contains ".fastq", so without
// this check it would be fopen'd as plain text and parsed into silent garbage.
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

// Parallel sort: split, sort the pieces on `threads` threads, then merge pairwise in place. Only
// <thread> and std::sort/std::inplace_merge, so it behaves the same on any compiler and any
// architecture -- and the result is the same order std::sort would have produced.
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
struct Batch { std::vector<std::string> seqs; int col; };

// Plain, uncompressed FASTQ. Decompression is not part of this program: measured on this box,
// gunzipping the 80 GB whole-genome .gz costs ~10 min while an assembler reads the resulting
// 188 GB plain FASTQ, parses it AND builds its graph in ~5 min. Paying that per run -- twice,
// once per pass -- bought nothing, so the input contract is uncompressed FASTQ, full stop.
static FILE* open_fastq(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[FATAL] fopen %s (input must be uncompressed FASTQ)\n", path.c_str()); std::abort(); }
    return f;
}

// ---- .svb compact reads (2-bit packed; snape uses only the ACGT sequence, no qual/name) ------
// A sample is a dir of shard.<k>.svb; each file = 32B header then sequential records:
//   len_tag u32 (low 31 = length L, bit31 = RAW); payload = RAW ? L ASCII bytes verbatim :
//   ceil(L/4) packed bytes, 2 bits/base (A=00 C=01 G=10 T=11), base i in byte i/4, high 2 bits
//   first. fopen'd (a real seekable file, NOT a pipe), so records fread exact-size -- no 1 MB
//   block/carry juggling. Sequential -> pass 1 and pass 2 re-read in identical (shard,rid) order.
static bool is_svb(const std::string& p) { return p.size() >= 4 && p.compare(p.size() - 4, 4, ".svb") == 0; }

// The four ASCII bases each of the 256 packed bytes stands for, as one word.
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
    return f;                                            // positioned at the first record
}

// Decode ONE record into `seq` (ACGT, or verbatim for a RAW record). false at clean EOF.
static bool svb_next(FILE* f, std::string& seq) {
    uint32_t tag;
    if (std::fread(&tag, 4, 1, f) != 1) return false;    // clean EOF
    const uint32_t L = tag & 0x7fffffffu;
    seq.resize(L);
    if (tag >> 31) {                                     // RAW escape (N / non-ACGT read)
        if (L && std::fread(&seq[0], 1, L, f) != L) { std::fprintf(stderr, "[FATAL] .svb raw truncated\n"); std::abort(); }
        for (char& c : seq) if (c >= 'a' && c <= 'z') c -= 32;   // the only .svb bytes not already ACGT
        return true;
    }
    const uint32_t nb = (L + 3) >> 2;
    static thread_local std::vector<uint8_t> pb;
    pb.resize(nb);
    if (nb && std::fread(pb.data(), 1, nb, f) != nb) { std::fprintf(stderr, "[FATAL] .svb packed truncated\n"); std::abort(); }
    // One store per PACKED BYTE: QUAD[b] is the four bases byte b encodes, written as a 4-byte word,
    // instead of a shift and a mask per base over 391 Gbp of reads.
    static const uint32_t* QUAD = quad_table();
    static const char LUT[4] = { 'A', 'C', 'G', 'T' };
    const uint32_t full = L >> 2;
    char* out = &seq[0];
    for (uint32_t i = 0; i < full; i++) std::memcpy(out + 4 * i, &QUAD[pb[i]], 4);
    for (uint32_t i = full << 2; i < L; i++) seq[i] = LUT[(pb[i >> 2] >> (6 - 2 * (i & 3))) & 3];
    return true;
}

}  // namespace   (the reader below is snape.h's, and the rest of this file is not)

// ---- the shared reads input (snape.h) ---------------------------------------------------------
// The same two readers above, behind one interface, so bubbles and junctions can take a .svb
// directly instead of being handed a FASTQ stream by a shell.
Reads reads_open(const std::string& path) {
    Reads r;
    r.svb = is_svb(path);
    r.f = r.svb ? open_svb(path) : open_fastq(path);
    return r;
}

bool reads_next(Reads& r, std::string& seq, std::string* name) {
    if (r.svb) {
        if (!svb_next(r.f, seq)) return false;
        if (name) { char b[32]; std::snprintf(b, sizeof b, "r%llu", (unsigned long long)r.idx); *name = b; }
        r.idx++;
        return true;
    }
    static thread_local char* buf = nullptr;                 // one line buffer per calling thread
    static thread_local size_t cap = 0;
    ssize_t n = getline(&buf, &cap, r.f);
    if (n <= 0) return false;
    while (n && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) n--;
    if (name) {
        const char* sp = (const char*)std::memchr(buf + 1, ' ', n > 0 ? (size_t)n - 1 : 0);
        name->assign(buf + 1, sp ? (size_t)(sp - buf - 1) : (size_t)(n > 0 ? n - 1 : 0));
    }
    n = getline(&buf, &cap, r.f);
    if (n <= 0) return false;
    while (n && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) n--;
    seq.assign(buf, (size_t)n);
    for (char& c : seq) if (c >= 'a' && c <= 'z') c -= 32;
    if (getline(&buf, &cap, r.f) <= 0) return false;         // '+'
    if (getline(&buf, &cap, r.f) <= 0) return false;         // quality, discarded
    r.idx++;
    return true;
}

void reads_close(Reads& r) { if (r.f) std::fclose(r.f); r.f = nullptr; }

namespace {

// ---- the build: ONE pass over the reads, exact maps -----------------------------------------
// One pass and not two, because an exact map knows a signature's MULTIPLICITY and its NEIGHBOURS at
// the same time. It costs holding adjacencies for error cores as well, since survivorship is not
// decided until the end, but the pruning rule stays exact: an edge is kept iff both endpoints
// survive. With homopolymer normalisation the map fits -- distinct is 78.7 M per colour and the two
// colours share one record.
using SigNodeMap = ankerl::unordered_dense::map<uint64_t, std::pair<uint32_t, uint32_t>>;

// One read -> node counts AND raw consecutive adjacencies, both strands, in a single walk.
// Both strands come out of ONE extract2() call, so the read is homopolymer-compressed once instead
// of twice and its reverse complement is never materialised. Cores within `margin` of either end of
// their strand are dropped, exactly as the two-call version dropped them.
void process_read_single(std::string& s, int level, int margin, int col,
                         SigNodeMap& nmap, EdgeMap& emap, uint64_t& cpu_ext, uint64_t& ncores) {
    // Reads arrive upper-cased: the FASTQ reader does it as it parses, the .svb decoder emits ACGT
    // and upper-cases the RAW records it does not decode. Doing it again here was a second pass over
    // every base of both samples -- 391 Gbp on hcc1954 -- to find nothing.
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int nf, nr; const uint64_t *stf, *enf, *sgf, *strv, *enr, *sgr;
    const uint64_t te = thr_ns();
    extract2(s.data(), (int)len, level, &nf, &stf, &enf, &sgf, &nr, &strv, &enr, &sgr);
    cpu_ext += thr_ns() - te;
    const uint64_t mg = (uint64_t)margin;
    auto chain = [&](int n, const uint64_t* st, const uint64_t* en, const uint64_t* sig) {
        uint64_t prev = 0; bool have = false;
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            ncores++;
            auto& c = nmap[sig[i]]; (col == 0 ? c.first : c.second)++;
            // Every consecutive pair, survivor or not -- the mm gate is applied after the pass.
            if (have && prev != sig[i]) { auto& x = emap[SEdge(prev, sig[i])]; (col == 0 ? x.first : x.second)++; }
            prev = sig[i]; have = true;
        }
    };
    chain(nf, stf, enf, sgf);
    chain(nr, strv, enr, sgr);
}

// Survivor lookup. out_nodes is sorted by signature, so the top SIGBITS bits of a signature index
// straight into the run of the array that could hold it: one offset table of 2^SIGBITS + 1 uint32s
// (67 MB) replaces a 70 M-entry hash set (is this core alive?) and a 70 M-entry hash map (which node
// id is it?), 2.4 GB between them, and answers both in one probe of about four entries.
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
    uint32_t find(uint64_t sig) const {                       // node id, or UINT32_MAX if dead
        const size_t b = (size_t)(sig >> (64 - BITS));
        for (uint32_t i = off[b], e = off[b + 1]; i < e; i++) if (v[i].sig == sig) return i;
        return UINT32_MAX;
    }
};


// The single pass itself: readers stream both samples, workers accumulate exact local node/edge maps
// and periodically merge them into 64 sig-sharded globals.
void run_pool_single(const std::vector<Shard>& shards, int level, int margin, int threads, int min_mult,
                     std::vector<SigCount>& out_nodes, std::vector<CEdge>& out_edges,
                     uint64_t& reads_n, uint64_t& reads_t) {
    constexpr int NSHARD = 64;
    SigNodeMap gnode[NSHARD]; EdgeMap gedge[NSHARD];
    std::mutex nmtx[NSHARD], emtx[NSHARD];
    const int nshards = (int)shards.size();
    const int nrd = nshards;
    const int nw = std::max(1, threads - nrd);
    std::atomic<size_t> next{0}, done{0};
    std::atomic<int> readers_left{nrd};
    std::atomic<uint64_t> rn{0}, rt{0}, emissions{0}, extract_ns{0};
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    constexpr size_t Q_CAP = 256, BATCH = 64;
    constexpr unsigned BLOCK = 1u << 20;

    auto reader = [&]() {
        std::vector<char> buf(BLOCK);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string carry; int col_local = 0; uint64_t shard_reads = 0;
        auto flush = [&]() {
            if (batch.empty()) return;
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(Batch{std::move(batch), col_local}); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        size_t si; std::string svseq;
        while ((si = next.fetch_add(1)) < shards.size()) {
            col_local = shards[si].col; shard_reads = 0;
            const bool svb = is_svb(shards[si].path);
            FILE* f = svb ? open_svb(shards[si].path) : open_fastq(shards[si].path);
            if (svb) {
                while (svb_next(f, svseq)) { shard_reads++; batch.emplace_back(std::move(svseq)); if (batch.size() >= BATCH) flush(); }
                std::fclose(f);
            } else {
                int line_no = 0; carry.clear(); int got;
                while ((got = (int)std::fread(buf.data(), 1, BLOCK, f)) > 0) {
                    const char* base = buf.data(); const char* end = base + got; const char* start = base;
                    for (const char* p = base; p < end; p++) if (*p == '\n') {
                        const size_t n = (size_t)(p - start);
                        if (line_no == 1) {
                            if (carry.empty()) batch.emplace_back(start, n);
                            else { carry.append(start, n); batch.push_back(std::move(carry)); carry.clear(); }
                            shard_reads++;
                            if (batch.size() >= BATCH) flush();
                        } else if (!carry.empty()) carry.clear();
                        line_no = (line_no + 1) & 3; start = p + 1;
                    }
                    if (start < end && line_no == 1) carry.append(start, end - start);
                }
                std::fclose(f);
            }
            flush();
            (col_local == 0 ? rn : rt) += shard_reads;
            std::fprintf(stderr, "[graph1] shard %zu/%d done (col=%d, %lu reads)\n",
                         (size_t)done.fetch_add(1) + 1, nshards, col_local, (unsigned long)shard_reads);
        }
        flush();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
    };

    auto worker = [&]() {
        uint64_t ext = 0, emis = 0;
        // 512 k entries, swept on the 1.5 M-read pair (wall, peak RSS, same graph.bin):
        // 4 M 95.9 s / 18.9 GB, 2 M 98.3 s / 15.9 GB, 512 k 94.1 s / 10.6 GB, 128 k 100.7 s / 9.9 GB
        // (the 2 M and 128 k runs shared the box with another job, so read their walls as +-4%).
        // A worker's own tables and the per-shard staging vectors are both bounded by this, so it is
        // the largest single term in peak memory; below ~512 k the extra merges start costing time.
        const char* fe = std::getenv("SNAPE_FLUSH");   // test hook: force mid-run merges at smoke scale
        const size_t FLUSH = fe ? (size_t)std::atoll(fe) : (size_t)(1u << 19);
        SigNodeMap ln; EdgeMap le; ln.reserve(1u << 20); le.reserve(1u << 20);
        // The per-shard staging vectors are reused: at this threshold they are the second largest
        // allocation a worker makes, and rebuilding them per flush was pure allocator traffic.
        std::vector<std::vector<std::pair<uint64_t, std::pair<uint32_t, uint32_t>>>> bn(NSHARD);
        std::vector<std::vector<std::pair<SEdge, std::pair<uint32_t, uint32_t>>>> be(NSHARD);
        auto flush = [&]() {
            for (auto& v : bn) v.clear();
            for (const auto& kv : ln) bn[kv.first & (NSHARD - 1)].push_back(kv);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (bn[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(nmtx[sh]);
                for (const auto& kv : bn[sh]) { auto& c = gnode[sh][kv.first]; c.first += kv.second.first; c.second += kv.second.second; }
            }
            ln.clear();
            for (auto& v : be) v.clear();
            for (const auto& kv : le) be[EdgeHash{}(kv.first) & (NSHARD - 1)].push_back(kv);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (be[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(emtx[sh]);
                for (const auto& kv : be[sh]) { auto& c = gedge[sh][kv.first]; c.first += kv.second.first; c.second += kv.second.second; }
            }
            le.clear();
        };
        while (true) {
            Batch b;
            {
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || readers_left.load() == 0; });
                if (q.empty()) break;
                b = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            const int c = b.col & 1;
            for (std::string& s : b.seqs) process_read_single(s, level, margin, c, ln, le, ext, emis);
            if (ln.size() >= FLUSH || le.size() >= FLUSH) flush();
        }
        flush();
        emissions += emis; extract_ns += ext;
    };

    { std::vector<std::thread> th;
      for (int i = 0; i < nrd; i++) th.emplace_back(reader);
      for (int i = 0; i < nw; i++) th.emplace_back(worker);
      for (auto& t : th) t.join(); }
    reads_n = rn; reads_t = rt;

    // ---- the mm gate, applied once, on exact counts ----
    // Shard by shard in parallel: the 64 maps are independent, so the gate, the id assignment and
    // the edge resolve are all per-shard work that was being done on one core.
    uint64_t distinct = 0, kept = 0;
    for (int sh = 0; sh < NSHARD; sh++) distinct += gnode[sh].size();
    {
        std::vector<std::vector<SigCount>> part(NSHARD);
        std::atomic<int> next_sh{0};
        std::vector<std::thread> th;
        for (int i = 0; i < threads; i++) th.emplace_back([&]{
            for (int sh = next_sh.fetch_add(1); sh < NSHARD; sh = next_sh.fetch_add(1)) {
                part[sh].reserve(gnode[sh].size() / 2);
                for (const auto& kv : gnode[sh]) {
                    const uint32_t n = kv.second.first, t = kv.second.second;
                    if (n < (uint32_t)min_mult && t < (uint32_t)min_mult) continue;   // survives in NEITHER colour
                    part[sh].push_back({kv.first, n, t});
                }
                SigNodeMap().swap(gnode[sh]);                                         // free as we go
            }
        });
        for (auto& t : th) t.join();
        for (int sh = 0; sh < NSHARD; sh++) kept += part[sh].size();
        out_nodes.clear(); out_nodes.reserve(kept);
        for (int sh = 0; sh < NSHARD; sh++) {
            out_nodes.insert(out_nodes.end(), part[sh].begin(), part[sh].end());
            std::vector<SigCount>().swap(part[sh]);       // released as it is consumed, so the two
        }                                                 // copies never both exist in full
    }
    psort(out_nodes, threads, [](const SigCount& a, const SigCount& b){ return a.sig < b.sig; });

    SigIndex idx; idx.build(out_nodes);
    uint64_t eraw = 0;
    for (int sh = 0; sh < NSHARD; sh++) eraw += gedge[sh].size();
    {
        std::vector<std::vector<CEdge>> part(NSHARD);
        std::atomic<int> next_sh{0};
        std::vector<std::thread> th;
        for (int i = 0; i < threads; i++) th.emplace_back([&]{
            for (int sh = next_sh.fetch_add(1); sh < NSHARD; sh = next_sh.fetch_add(1)) {
                part[sh].reserve(gedge[sh].size() / 3);
                for (const auto& kv : gedge[sh]) {
                    // Pass 2 breaks the chain at a dead core, so an edge exists iff BOTH ends
                    // survived; resolving to node ids here is the same lookup, so it costs nothing.
                    const uint32_t a = idx.find(kv.first.first), b = idx.find(kv.first.second);
                    if (a == UINT32_MAX || b == UINT32_MAX || a == b) continue;
                    part[sh].push_back({a, b, kv.second.first, kv.second.second});
                }
                EdgeMap().swap(gedge[sh]);
            }
        });
        for (auto& t : th) t.join();
        size_t ne = 0;
        for (int sh = 0; sh < NSHARD; sh++) ne += part[sh].size();
        out_edges.clear(); out_edges.reserve(ne);
        for (int sh = 0; sh < NSHARD; sh++) {
            out_edges.insert(out_edges.end(), part[sh].begin(), part[sh].end());
            std::vector<CEdge>().swap(part[sh]);
        }
    }
    std::fprintf(stderr, "[graph1] core emissions=%.3f M, of which %.1f CPU-s in extract\n",
                 emissions.load() / 1e6, extract_ns.load() / 1e9);
    std::fprintf(stderr, "[graph1] distinct sigs=%.3f M -> kept=%.3f M (%.2f%%; discarded as read error %.3f M)\n",
                 distinct / 1e6, kept / 1e6, distinct ? 100.0 * kept / distinct : 0.0, (distinct - kept) / 1e6);
    std::fprintf(stderr, "[graph1] raw adjacencies=%.3f M -> both-ends-survive=%.3f M (%.2f%%)\n",
                 eraw / 1e6, out_edges.size() / 1e6, eraw ? 100.0 * out_edges.size() / eraw : 0.0);
}



}  // namespace

// ---- the stage: reads -> counted cores and adjacencies -> self-contained graph.bin ----
// graph.bin schema: u64 K, u64 sig[K], u32 ncov[K], u32 tcov[K], u64 E, CEdge[E].
//   sig[K]  = MurmurHash64A over the core's homopolymer-compressed bytes = the node's identity;
//             colour is derived from ncov/tcov on demand and never stored.
//   CEdge   = {u32 src, u32 dst, u32 ncov, u32 tcov}: directed adjacency with per-colour traversal
//             counts (how many normal / tumor reads cross it). No core sequence is stored anymore.
void run_graph_reads(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir,
                     int level, int margin, int min_mult, int threads) {
    mkdir(dir.c_str(), 0755);
    const std::vector<std::string> ns = collect_shards(normal_fq), ts = collect_shards(tumor_fq);
    std::vector<Shard> shards_n, shards_t, shards_all;
    for (const auto& p : ns) shards_n.push_back({p, 0});
    for (const auto& p : ts) shards_t.push_back({p, 1});
    shards_all = shards_n; shards_all.insert(shards_all.end(), shards_t.begin(), shards_t.end());
    std::fprintf(stderr, "[graph] normal=%zu shards, tumor=%zu shards, threads=%d min_mult=%d (merged: bloom-count + fill + emit)\n",
                 ns.size(), ts.size(), threads, min_mult);
    Tic tic("graph");

    // ---- per-phase CPU-second + wall accounting. wall via steady_clock; CPU via getrusage
    // (SELF = all in-process threads incl. readers/workers/main). Snapshots taken at each phase
    // boundary; deltas attribute cost to that phase. ----
    struct Snap { double wall, self, rss, peak; };
    const auto torig = std::chrono::steady_clock::now();
    auto snap = [&]() -> Snap { return { std::chrono::duration<double>(std::chrono::steady_clock::now() - torig).count(), self_cpu_s(), rss_gb(), peak_gb() }; };
    const Snap P0 = snap();

    // ---- the build: ONE pass over the reads, exact maps, nodes and edges together ----
    std::vector<SigCount> nodes; std::vector<CEdge> ce; uint64_t nr_n = 0, nr_t = 0;
    run_pool_single(shards_all, level, margin, threads, min_mult, nodes, ce, nr_n, nr_t);
    const uint64_t K = nodes.size();
    const Snap P1 = snap();
    tic.mark("pool_both");
    std::fprintf(stderr, "[graph]   survivors=%lu\n", (unsigned long)K);

    // ---- node sigs (node id = sig-sorted index; the sig IS the stored node identity now, colour is
    // derived from ncov/tcov on demand and never materialised). The edges came back already
    // resolved to node ids -- the survivor index did it with the same probe that tested liveness. ----
    std::vector<uint64_t> sig(K); std::vector<uint32_t> ncov(K), tcov(K);
    for (uint64_t i = 0; i < K; i++) { const SigCount& c = nodes[i]; sig[i] = c.sig; ncov[i] = c.n; tcov[i] = c.t; }
    std::vector<SigCount>().swap(nodes);
    const Snap P2 = snap();

    psort(ce, threads, [](const CEdge& x, const CEdge& y){ return x.src != y.src ? x.src < y.src : x.dst < y.dst; });
    const Snap P4 = snap();
    tic.mark("sort_edges");

    // ---- write self-contained graph.bin: u64 K, u64 sig[K], u32 ncov[K], u32 tcov[K], u64 E, CEdge[E] ----
    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen(gp + ".tmp", "wb");
    uint64_t Kw = K; wr(fg, &Kw, 8);
    wr(fg, sig.data(), K * 8); wr(fg, ncov.data(), K * 4); wr(fg, tcov.data(), K * 4);
    uint64_t E = ce.size(); wr(fg, &E, 8); wr(fg, ce.data(), E * sizeof(CEdge));
    atomic_close(fg, gp + ".tmp", gp);
    const Snap P5 = snap();
    tic.mark("write");
    std::fprintf(stderr, "[graph] normal_reads=%lu tumor_reads=%lu nodes=%lu edges=%lu -> %s/graph.bin\n",
                 (unsigned long)nr_n, (unsigned long)nr_t, (unsigned long)K, (unsigned long)E, dir.c_str());

    // ---- per-component CPU-second + execution-efficiency table -------------------------------
    // CPU = Δself (all in-process threads). wall = phase wall.
    // efficiency = CPU / (ncpu * wall) = fraction of the box's cores the phase actually kept busy
    // (1.0 = every core saturated for the whole phase; <1.0 = idle cores / serial sections).
    cpu_set_t cs; CPU_ZERO(&cs);   // efficiency denominator = cores ACTUALLY allotted (taskset/affinity), not the box's total
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
}
