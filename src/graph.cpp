#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1   // expose sched_getaffinity under strict -std=c++17
#endif
#include "graph.h"
#include "samples.h"
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

// Per-edge key (prev_sig, cur_sig). Hash mixes the sig pair; not tagged
// is_avalanching so unordered_dense re-mixes it (prev*golden+cur alone is not fully avalanched).
struct EdgeHash { uint64_t operator()(const SEdge& e) const { return e.first * 0x9e3779b97f4a7c15ULL + e.second; } };
using EdgeSet = ankerl::unordered_dense::set<SEdge, EdgeHash>;

// Edges carry NO per-sample counts. Nothing downstream ever read them: `ecov` appears in no
// expression in bubbles.cpp, junctions.cpp, call.cpp or pick.cpp, and both graph.bin readers stop
// after `edst`. At S=21 that was 21 bytes of arena on every edge -- about 29% of the edge table and
// 17% of peak memory -- kept alive for a consumer that does not exist. An edge is a fact about the
// graph's shape, and its shape is all bubbles and junctions ask for.

// Counts saturate: no stage reads their magnitude, only zero / non-zero and the min_mult floor.
static inline void bump(cov_t& c) { c += (c != COV_MAX); }

// A surviving core's RAW bytes, 2-bit packed (A=00 C=01 G=10 T=11). bubbles needs a core's bytes to
// measure an arm in bases and to derive its minimizers, and the only other place those bytes exist
// is the reads -- which is why bubbles used to re-read every FASTQ: measured, 2.4 TB of scanning to
// recover 49 MB of sequence. Carrying them in graph.bin costs ~3.7 GiB and turns that scan into a
// 6 GiB read. The bytes and not the minimizers, so k and w stay parameters of bubbles.
//
// A core holding anything outside ACGT is stored verbatim instead, flagged in the high bit of `len`:
// packing an N as an A would change the homopolymer compression and every minimizer downstream.
struct SeqRef { uint64_t off = 0; uint32_t len = 0; };     // len == 0: nothing stored yet

static inline int b2(char c) { return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1; }

static SeqRef seq_put(std::vector<uint8_t>& a, const char* s, uint32_t n) {
    SeqRef r; r.off = a.size(); r.len = n;
    for (uint32_t i = 0; i < n; i++) if (b2(s[i]) < 0) { r.len |= SEQ_RAW; break; }
    if (r.len & SEQ_RAW) { a.insert(a.end(), s, s + n); return r; }
    a.resize(r.off + (n + 3) / 4, 0);
    for (uint32_t i = 0; i < n; i++) a[r.off + (i >> 2)] |= (uint8_t)(b2(s[i]) << (6 - 2 * (i & 3)));
    return r;
}

// Move an already-packed core between arenas: worker -> sample -> graph, each time only for cores
// that have earned the space by clearing the gate.
static SeqRef seq_copy(std::vector<uint8_t>& dst, const uint8_t* src, SeqRef r) {
    SeqRef o; o.off = dst.size(); o.len = r.len;
    const uint64_t nb = seq_bytes(r.len);
    dst.insert(dst.end(), src + r.off, src + r.off + nb);
    return o;
}
static inline void addsat(cov_t& d, cov_t v) { const uint32_t t = (uint32_t)d + v; d = t > COV_MAX ? COV_MAX : (cov_t)t; }

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
Reads reads_open(const std::string& path, const std::string& sample) {
    Reads r;
    r.svb = is_svb(path);
    r.f = r.svb ? open_svb(path) : open_fastq(path);
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
    static thread_local char* buf = nullptr;                 // one line buffer per calling thread
    static thread_local size_t cap = 0;
    ssize_t n = getline(&buf, &cap, r.f);
    if (n <= 0) return false;
    while (n && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) n--;
    if (name) {
        const char* sp = (const char*)std::memchr(buf + 1, ' ', n > 0 ? (size_t)n - 1 : 0);
        name->assign(r.tag);
        name->append(buf + 1, sp ? (size_t)(sp - buf - 1) : (size_t)(n > 0 ? n - 1 : 0));
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
using SigNodeMap = ankerl::unordered_dense::map<uint64_t, uint32_t>;
struct NodeTab { SigNodeMap m; std::vector<cov_t> a; std::vector<SeqRef> r; std::vector<uint8_t> s; };

static inline cov_t* node_slot(NodeTab& t, uint64_t k, int S) {
    const auto r = t.m.try_emplace(k, (uint32_t)(t.a.size() / (size_t)S));
    if (r.second) { t.a.resize(t.a.size() + (size_t)S, 0); t.r.emplace_back(); }
    return t.a.data() + (size_t)r.first->second * (size_t)S;
}

// One sample's own tables: a count per signature, its bytes, and the bare existence of an adjacency.
// No colour, because inside one sample there is only one colour.
struct Core { cov_t n = 0; SeqRef s; };
using SigCntMap = ankerl::unordered_dense::map<uint64_t, Core>;
struct ShTab { SigCntMap m; std::vector<uint8_t> a; };
struct Local { SigCntMap n; EdgeSet e; std::vector<uint8_t> a; };

// One read -> its cores' counts AND its raw consecutive adjacencies, both strands, in a single walk.
// Both strands come out of ONE extract2() call, so the read is homopolymer-compressed once instead
// of twice and its reverse complement is never materialised. Cores within `margin` of either end of
// their strand are dropped, exactly as the two-call version dropped them.
void process_read_single(std::string& s, int level, int margin, Local& lo,
                         uint64_t& cpu_ext, uint64_t& ncores, uint64_t& ov_s, uint64_t& ov_n) {
    // Reads arrive upper-cased: the FASTQ reader does it as it parses, the .svb decoder emits ACGT
    // and upper-cases the RAW records it does not decode.
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int nf, nr; const uint64_t *stf, *enf, *sgf, *strv, *enr, *sgr;
    const uint64_t te = thr_ns();
    extract2(s.data(), (int)len, level, &nf, &stf, &enf, &sgf, &nr, &strv, &enr, &sgr);
    cpu_ext += thr_ns() - te;
    const uint64_t mg = (uint64_t)margin;
    // extract2 gives the reverse strand's cores in the REVERSE COMPLEMENT's own frame, and never
    // materialises it -- that is the point of doing both strands from one compression. Storing a
    // core's bytes needs the string those offsets index, so the reverse complement is built once per
    // read, and only when a reverse-strand core actually has to be stored. It is a linear pass at
    // ~0.3 ns/base against extract2's ~15, so it does not show up next to the parse.
    static thread_local std::string cb;                 // the core's own bytes, built per stored core
    auto chain = [&](bool rev, int n, const uint64_t* st, const uint64_t* en, const uint64_t* sig) {
        // Consecutive cores of one read overlap in that read's own bytes, so summing core lengths
        // over an arm overcounts by (cores-1) x overlap. Measured here, where the read is in hand,
        // and carried in the header: bubbles used to re-read every FASTQ to get this one number.
        for (int i = 0; i + 1 < n; i++) if (en[i] > st[i + 1]) { ov_s += en[i] - st[i + 1]; ov_n++; }
        uint64_t prev = 0; bool have = false;
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            ncores++;
            const auto r = lo.n.try_emplace(sig[i], Core{});
            // The LAST core of a strand has no trustworthy `en`: extract's emit() reads one byte
            // past the core, and for that one core it lands on whatever the buffer held before, so
            // the span can even come out inverted. The count still goes in -- that stale `en` is
            // what the margin test has always used -- but its BYTES are skipped, and at 18x another
            // read carries the same core in its interior and supplies them.
            if (r.second && en[i] > st[i] && en[i] <= len) {
                const uint32_t w = (uint32_t)(en[i] - st[i]);
                cb.assign(s, (size_t)st[i], w);
                // The reverse strand's cores are spans of the FORWARD read whose sequence is that
                // span reverse-complemented -- that is the string whose compression was hashed.
                if (rev) {
                    for (uint32_t a = 0, b = w - 1; a < b; a++, b--) std::swap(cb[a], cb[b]);
                    for (char& c : cb) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : c;
                }
                r.first->second.s = seq_put(lo.a, cb.data(), w);
            }
            bump(r.first->second.n);
            // Every consecutive pair, survivor or not -- the mm gate fires when the sample ends.
            if (have && prev != sig[i]) lo.e.insert(SEdge(prev, sig[i]));
            prev = sig[i]; have = true;
        }
    };
    chain(false, nf, stf, enf, sgf);
    chain(true, nr, strv, enr, sgr);
}

// Survivor lookup. out_nodes is sorted by signature, so the top BITS bits of a signature index
// straight into the run of the array that could hold it: one offset table of 2^BITS + 1 uint32s
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


constexpr int NSHARD = 64;

// One sample's reads -> that sample's own 64-sharded exact tables.
//
// ONE reader thread and threads-1 workers. Measured on this box, on a real 114 GiB read set: a
// single reader does 1622 Mbp/s -- block read, line split and upper-case included -- while one core
// of extract2() does 65 Mbp/s. One reader feeds about 25 workers, so the old `threads/2` reader cap
// left half the machine copying bytes: at --threads 32 it ran 16 readers and 16 workers and the
// process sat at 1650% CPU on a 32-core box.
void count_sample(const std::vector<Shard>& shards, int level, int margin, int threads, int min_mult,
                  ShTab* psn, EdgeSet* pse, std::mutex* nmtx, std::mutex* emtx,
                  std::atomic<uint64_t>& emissions, std::atomic<uint64_t>& extract_ns,
                  std::atomic<uint64_t>& ov_sum, std::atomic<uint64_t>& ov_cnt, uint64_t& reads_here) {
    const int nw = std::max(1, threads - 1);
    std::atomic<size_t> next{0};
    std::atomic<uint64_t> rcnt{0};
    std::atomic<int> readers_left{1};
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<std::vector<std::string>> q;
    constexpr size_t Q_CAP = 256, BATCH = 64;
    constexpr unsigned BLOCK = 1u << 20;

    auto reader = [&]() {
        std::vector<char> buf(BLOCK);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string carry;
        auto push = [&]() {
            if (batch.empty()) return;
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(std::move(batch)); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        size_t si; std::string svseq;
        while ((si = next.fetch_add(1)) < shards.size()) {
            uint64_t shard_reads = 0;
            const bool svb = is_svb(shards[si].path);
            FILE* f = svb ? open_svb(shards[si].path) : open_fastq(shards[si].path);
            if (svb) {
                while (svb_next(f, svseq)) { shard_reads++; batch.emplace_back(std::move(svseq)); if (batch.size() >= BATCH) push(); }
            } else {
                int line_no = 0; carry.clear(); int got;
                while ((got = (int)std::fread(buf.data(), 1, BLOCK, f)) > 0) {
                    const char* base = buf.data(); const char* end = base + got; const char* start = base;
                    for (const char* p = base; p < end; p++) if (*p == '\n') {
                        const size_t n = (size_t)(p - start);
                        if (line_no == 1) {
                            if (carry.empty()) batch.emplace_back(start, n);
                            else { carry.append(start, n); batch.push_back(std::move(carry)); carry.clear(); }
                            // Upper-case HERE, where a FASTQ sequence line is formed: a soft-masked
                            // FASTQ would otherwise reach extract() as lower-case bytes, hpc() splits
                            // at every case boundary and the cores come out with other signatures.
                            for (char& c : batch.back()) if (c >= 'a' && c <= 'z') c -= 32;
                            shard_reads++;
                            if (batch.size() >= BATCH) push();
                        } else if (!carry.empty()) carry.clear();
                        line_no = (line_no + 1) & 3; start = p + 1;
                    }
                    if (start < end && line_no == 1) carry.append(start, end - start);
                }
            }
            std::fclose(f);
            push();
            rcnt += shard_reads;
        }
        push();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
    };

    auto worker = [&]() {
        uint64_t ext = 0, emis = 0;
        // How many entries a worker accumulates before merging into this sample's shard tables. It
        // bounds a worker's own two tables AND its per-shard staging vectors, and it is a speed dial
        // too: a bigger table dedupes more repeats before the locked, cache-missing merge sees them.
        const char* fe = std::getenv("SNAPE_FLUSH");
        const size_t FLUSH = fe ? (size_t)std::atoll(fe) : (size_t)(1u << 20);
        uint64_t ovs = 0, ovn = 0;
        Local lo; lo.n.reserve(1u << 20); lo.e.reserve(1u << 20);
        std::vector<std::vector<std::pair<uint64_t, Core>>> bn(NSHARD);
        std::vector<std::vector<SEdge>> be(NSHARD);
        auto flush = [&]() {
            for (auto& v : bn) v.clear();
            for (const auto& kv : lo.n) bn[kv.first & (NSHARD - 1)].push_back(kv);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (bn[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(nmtx[sh]);
                for (const auto& kv : bn[sh]) {
                    Core& e = psn[sh].m[kv.first];
                    addsat(e.n, kv.second.n);
                    // The bytes are promoted the moment the count clears the gate, and never before:
                    // half of a sample's distinct signatures are read errors and cost nothing here.
                    if (e.s.len == 0 && (uint32_t)e.n >= (uint32_t)min_mult)
                        e.s = seq_copy(psn[sh].a, lo.a.data(), kv.second.s);
                }
            }
            lo.n.clear(); lo.a.clear();
            for (auto& v : be) v.clear();
            for (const auto& e : lo.e) be[EdgeHash{}(e) & (NSHARD - 1)].push_back(e);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (be[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(emtx[sh]);
                for (const auto& e : be[sh]) pse[sh].insert(e);
            }
            lo.e.clear();
        };
        while (true) {
            std::vector<std::string> b;
            {
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || readers_left.load() == 0; });
                if (q.empty()) break;
                b = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            for (std::string& s : b) process_read_single(s, level, margin, lo, ext, emis, ovs, ovn);
            if (lo.n.size() >= FLUSH || lo.e.size() >= FLUSH) flush();
        }
        flush();
        emissions += emis; extract_ns += ext; ov_sum += ovs; ov_cnt += ovn;
    };

    { std::vector<std::thread> th;
      th.emplace_back(reader);
      for (int i = 0; i < nw; i++) th.emplace_back(worker);
      for (auto& t : th) t.join(); }
    reads_here = rcnt.load();
}

static inline bool survives(const ShTab* psn, uint64_t sig, int min_mult) {
    const SigCntMap& m = psn[sig & (NSHARD - 1)].m;
    const auto it = m.find(sig);
    return it != m.end() && (uint32_t)it->second.n >= (uint32_t)min_mult;
}

// The mm gate, applied to ONE sample on that sample's own exact counts, and the merge into the
// S-plane graph. Two phases, because phase 2 reads every shard of `psn` while phase 1 only writes
// its own: freeing inside phase 1 would race with phase 2's lookups.
void merge_sample(ShTab* psn, EdgeSet* pse, NodeTab* gnode, EdgeSet* gedge,
                  int col, int S, int min_mult, int threads, uint64_t& seen, uint64_t& kept) {
    const int nt = std::max(1, threads);
    std::atomic<uint64_t> nseen{0}, npass{0};
    {
        std::atomic<int> next_sh{0};
        std::vector<std::thread> th;
        for (int i = 0; i < nt; i++) th.emplace_back([&]{
            for (int sh = next_sh.fetch_add(1); sh < NSHARD; sh = next_sh.fetch_add(1)) {
                uint64_t a = 0, b = 0;
                for (const auto& kv : psn[sh].m) {
                    a++;
                    if ((uint32_t)kv.second.n < (uint32_t)min_mult) continue;  // read error in THIS sample
                    b++;
                    node_slot(gnode[sh], kv.first, S)[col] = kv.second.n;
                    // The first sample to promote a core donates its bytes; the rest are the same
                    // 217 bases of the same locus, so the graph keeps one copy.
                    const uint32_t slot = gnode[sh].m.find(kv.first)->second;
                    if (gnode[sh].r[slot].len == 0)
                        gnode[sh].r[slot] = seq_copy(gnode[sh].s, psn[sh].a.data(), kv.second.s);
                }
                nseen += a; npass += b;
            }
        });
        for (auto& t : th) t.join();
    }
    {   // An edge is kept iff BOTH ends cleared the gate in this sample: the pair was seen in this
        // sample's reads, so this sample's counts are the ones entitled to decide it. The per-sample
        // and the global edge tables shard on the same hash, so a thread only ever writes its own.
        std::atomic<int> next_sh{0};
        std::vector<std::thread> th;
        for (int i = 0; i < nt; i++) th.emplace_back([&]{
            for (int sh = next_sh.fetch_add(1); sh < NSHARD; sh = next_sh.fetch_add(1))
                for (const auto& e : pse[sh])
                    if (survives(psn, e.first, min_mult) && survives(psn, e.second, min_mult))
                        gedge[sh].insert(e);
        });
        for (auto& t : th) t.join();
    }
    for (int sh = 0; sh < NSHARD; sh++) {
        SigCntMap().swap(psn[sh].m); std::vector<uint8_t>().swap(psn[sh].a); EdgeSet().swap(pse[sh]);
    }
    seen = nseen.load(); kept = npass.load();
}

// The build: ONE SAMPLE AT A TIME. Each sample's exact tables are gated at min_mult and merged into
// the S-plane graph, then thrown away, so the merged map only ever holds signatures some sample saw
// at least min_mult times.
//
// Gating after the merge -- which is what this used to do -- meant holding every sample's read
// errors simultaneously, and a read error is seen once, so they are most of what a sample's table
// contains. Measured at S=21 over 1.2 Tbp: RSS grew 3.72 GiB/min, linearly, never flattening, and
// the kernel killed the process at 141 GiB before it had finished reading.
//
// What this changes in the RESULT: `cov[j] == 0` now means "sample j did not reach min_mult", not
// "sample j never saw it". The somatic test reads plane 0, so somatic now means the normal did not
// reach the threshold -- the same rule the tumour planes were always held to. At 18x a true germline
// core has ~1e-6 chance of landing under three reads, so one or two normal reads at a somatic locus
// is a read error, and vetoing the call on it was accepting one sample's noise as another's
// evidence. The germline question is asked again, properly, on the normal's own alignment in `call`.
void run_pool_single(const std::vector<Shard>& shards, int level, int margin, int threads, int min_mult,
                     int S, std::vector<SigCount>& out_nodes, std::vector<cov_t>& out_ncov,
                     std::vector<CEdge>& out_edges, std::vector<SeqRef>& out_sref,
                     std::vector<uint8_t>& out_seq, uint64_t& ov_s, uint64_t& ov_c,
                     std::vector<uint64_t>& reads_per) {
    NodeTab gnode[NSHARD]; EdgeSet gedge[NSHARD];
    ShTab psn[NSHARD]; EdgeSet pse[NSHARD];
    std::mutex nmtx[NSHARD], emtx[NSHARD];
    std::atomic<uint64_t> emissions{0}, extract_ns{0}, ov_sum{0}, ov_cnt{0};
    reads_per.assign((size_t)S, 0);
    uint64_t distinct = 0, kept = 0;

    for (int col = 0; col < S; col++) {
        std::vector<Shard> mine;
        for (const Shard& s : shards) if (s.col == col) mine.push_back(s);
        uint64_t rd = 0, seen = 0, pass = 0;
        count_sample(mine, level, margin, threads, min_mult, psn, pse, nmtx, emtx, emissions,
                     extract_ns, ov_sum, ov_cnt, rd);
        merge_sample(psn, pse, gnode, gedge, col, S, min_mult, threads, seen, pass);
        reads_per[col] = rd;
        distinct += seen; kept += pass;
        std::fprintf(stderr, "[graph1] sample %d/%d: %llu reads, %.3f M distinct -> %.3f M >= %d (%.1f%%)\n",
                     col + 1, S, (unsigned long long)rd, seen / 1e6, pass / 1e6, min_mult,
                     seen ? 100.0 * pass / seen : 0.0);
    }

    // Flatten. Every entry in gnode already cleared the gate in some sample, so there is no second
    // gate to apply here -- only the id assignment.
    {
        std::vector<std::vector<SigCount>> part(NSHARD);
        std::vector<std::vector<cov_t>> pcov(NSHARD);
        std::vector<std::vector<SeqRef>> pref(NSHARD);
        std::atomic<int> next_sh{0};
        std::vector<std::thread> th;
        for (int i = 0, nt = std::max(1, threads); i < nt; i++) th.emplace_back([&]{
            for (int sh = next_sh.fetch_add(1); sh < NSHARD; sh = next_sh.fetch_add(1)) {
                part[sh].reserve(gnode[sh].m.size());
                for (const auto& kv : gnode[sh].m) {
                    const cov_t* c = gnode[sh].a.data() + (size_t)kv.second * (size_t)S;
                    part[sh].push_back({kv.first, (uint32_t)(pcov[sh].size() / (size_t)S)});
                    pcov[sh].insert(pcov[sh].end(), c, c + S);
                    pref[sh].push_back(gnode[sh].r[kv.second]);
                }
                SigNodeMap().swap(gnode[sh].m);
                std::vector<cov_t>().swap(gnode[sh].a);
                std::vector<SeqRef>().swap(gnode[sh].r);
            }
        });
        for (auto& t : th) t.join();
        size_t K = 0;
        for (int sh = 0; sh < NSHARD; sh++) K += part[sh].size();
        out_nodes.clear(); out_nodes.reserve(K);
        out_ncov.clear(); out_ncov.reserve(K * (size_t)S);
        out_sref.clear(); out_sref.reserve(K);
        out_seq.clear();
        for (int sh = 0; sh < NSHARD; sh++) {
            const uint32_t base = (uint32_t)(out_ncov.size() / (size_t)S);
            const uint64_t sbase = out_seq.size();
            for (SigCount c : part[sh]) { c.slot += base; out_nodes.push_back(c); }
            for (SeqRef r : pref[sh]) { r.off += sbase; out_sref.push_back(r); }
            out_ncov.insert(out_ncov.end(), pcov[sh].begin(), pcov[sh].end());
            out_seq.insert(out_seq.end(), gnode[sh].s.begin(), gnode[sh].s.end());
            std::vector<SigCount>().swap(part[sh]);
            std::vector<cov_t>().swap(pcov[sh]);
            std::vector<SeqRef>().swap(pref[sh]);
            std::vector<uint8_t>().swap(gnode[sh].s);
        }
    }
    psort(out_nodes, threads, [](const SigCount& a, const SigCount& b){ return a.sig < b.sig; });
    SigIndex idx; idx.build(out_nodes);

    uint64_t eraw = 0;
    for (int sh = 0; sh < NSHARD; sh++) eraw += gedge[sh].size();
    {
        std::vector<std::vector<CEdge>> part(NSHARD);
        std::atomic<int> next_sh{0};
        std::vector<std::thread> th;
        for (int i = 0, nt = std::max(1, threads); i < nt; i++) th.emplace_back([&]{
            for (int sh = next_sh.fetch_add(1); sh < NSHARD; sh = next_sh.fetch_add(1)) {
                part[sh].reserve(gedge[sh].size());
                for (const auto& e : gedge[sh]) {
                    const uint32_t a = idx.find(e.first), b = idx.find(e.second);
                    if (a == UINT32_MAX || b == UINT32_MAX || a == b) continue;
                    part[sh].push_back({a, b});
                }
                EdgeSet().swap(gedge[sh]);
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
    std::fprintf(stderr, "[graph1] per-sample distinct sigs=%.3f M -> passed mm=%d: %.3f M (%.2f%%)\n",
                 distinct / 1e6, min_mult, kept / 1e6, distinct ? 100.0 * kept / distinct : 0.0);
    ov_s = ov_sum.load(); ov_c = ov_cnt.load();
    std::fprintf(stderr, "[graph1] merged nodes=%.3f M | raw adjacencies=%.3f M -> kept=%.3f M (%.2f%%)\n",
                 out_nodes.size() / 1e6, eraw / 1e6, out_edges.size() / 1e6,
                 eraw ? 100.0 * out_edges.size() / eraw : 0.0);
    std::fprintf(stderr, "[graph1] core sequence=%.3f GiB, ardisik cekirdek ortusme %.1f bp\n",
                 out_seq.size() / 1073741824.0, ov_c ? (double)ov_s / (double)ov_c : 0.0);
}



}  // namespace

// ---- the stage: reads -> counted cores and adjacencies -> self-contained graph.bin ----
// graph.bin v2 schema (little-endian throughout):
//   char[4] "SNPG" | u32 version=4 | u32 S | u32 reserved | u64 ov_sum | u64 ov_cnt
//   u64 K | u64 sig[K] | cov_t cov[S][K]
//   u64 sbytes | u32 slen[K] | u8 seq[sbytes]
//   u64 E | u32 esrc[E] | u32 edst[E]
//
//   slen[i]  the core's RAW length, high bit set when its bytes are stored verbatim rather than
//            2-bit packed (a core holding an N). seq[] is those bytes, concatenated in node order.
//   ov_*     the summed overlap of consecutive cores within a read, and how many pairs it came
//            from. bubbles needs the mean to compare arms of different core counts, and used to
//            re-read every FASTQ to measure it.
//
// v1 had NO magic and NO version -- it began with K -- so a v1 file and a v2 file were
// indistinguishable except by size arithmetic, and `run` resumes on nothing but the presence of the
// file. The header exists so that a stale graph in a work dir is an abort and not silent garbage.
//
//   sig[K]     MurmurHash64A over the core's homopolymer-compressed bytes = the node's identity.
//   cov[s][i]  reads of sample s carrying node i. Plane 0 is ALWAYS the normal, planes 1..S-1 the
//              tumour lines in sheet order. COLOUR IS STILL DERIVED, NEVER STORED: cov[0][i]==0 is
//              somatic, and which lines carry it is which planes are non-zero.
//   Nodes are SoA and so are edges, because v1's AoS edge record forced junctions.cpp to slurp the
//   whole edge section as a flat u32 array with a hard-coded stride.
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

    // ---- per-phase CPU-second + wall accounting. wall via steady_clock; CPU via getrusage
    // (SELF = all in-process threads incl. readers/workers/main). Snapshots taken at each phase
    // boundary; deltas attribute cost to that phase. ----
    struct Snap { double wall, self, rss, peak; };
    const auto torig = std::chrono::steady_clock::now();
    auto snap = [&]() -> Snap { return { std::chrono::duration<double>(std::chrono::steady_clock::now() - torig).count(), self_cpu_s(), rss_gb(), peak_gb() }; };
    const Snap P0 = snap();

    // ---- the build: ONE pass over the reads, exact maps, nodes and edges together ----
    std::vector<SigCount> nodes; std::vector<CEdge> ce;
    std::vector<cov_t> ncov; std::vector<uint64_t> nreads;
    std::vector<SeqRef> sref; std::vector<uint8_t> sbuf; uint64_t ov_s = 0, ov_c = 0;
    run_pool_single(shards_all, level, margin, threads, min_mult, S, nodes, ncov, ce, sref, sbuf,
                    ov_s, ov_c, nreads);
    const uint64_t K = nodes.size();
    const Snap P1 = snap();
    tic.mark("pool_both");
    std::fprintf(stderr, "[graph]   survivors=%lu\n", (unsigned long)K);

    // ---- node sigs (node id = sig-sorted index; the sig IS the stored node identity now, colour is
    // derived from ncov/tcov on demand and never materialised). The edges came back already
    // resolved to node ids -- the survivor index did it with the same probe that tested liveness. ----
    const Snap P2 = snap();          // nodes come back sig-sorted: SigIndex needs that inside the pool

    psort(ce, threads, [](const CEdge& x, const CEdge& y){ return x.src != y.src ? x.src < y.src : x.dst < y.dst; });
    const Snap P4 = snap();
    tic.mark("sort_edges");

    // ---- write self-contained graph.bin v2 ----
    // The count planes are GATHERED one at a time through the sorted nodes' slots, so only one plane
    // (K * sizeof(cov_t)) is ever materialised on top of the arena -- at 21 samples the alternative,
    // building the whole SoA block first, would hold a second 3.5 GB copy.
    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen(gp + ".tmp", "wb");
    const uint32_t ver = 4, Sw = (uint32_t)S, resv = 0;
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
    // Core sequence, in node order, WHILE `nodes` is still alive -- it is released just below, and
    // gathering through its slots afterwards is what a first version of this did and segfaulted on.
    // Written through a buffer: one fwrite per core would be K of them, and K is 68 million.
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
