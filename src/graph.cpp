#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1   // expose sched_getaffinity under strict -std=c++17
#endif
#include "graph.h"
#include "btllib/bloom_filter.hpp"   // vendored: plain (atomic, lock-free) BloomFilter
// extract is declared in svp.h (provided by src/extract.cpp over biolcp).
#include "posmap.h"       // SVP_POSAWARE: reference position map (sig -> loci), consumer
#include "placement.h"    // SVP_POSAWARE: lookup-based read placement (no DP alignment)

#include <sys/stat.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <sched.h>
#include <dirent.h>
#include <algorithm>
#include <array>
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

// ---- fail-fast binary checkpoint IO (declared in svp.h; reused by compact.cpp) ----
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

namespace {
using SEdge = std::pair<uint64_t, uint64_t>;        // signature edge (edges.bin record)
// Per-emission sigs are collected as bare uint64_t: each pool pass is a single colour
// (normal OR tumor), so a per-record colour byte is pure waste. The run length of equal
// sigs after sort IS the n-count (normal pass) / t-count (tumor pass); the other stays 0.

// ASCII complement lookup. Default 'N' (non-ACGT collapses to 'N'); A/T and C/G swapped.
// Used by the remap+reverse RC pass and by canonical-window hashing.
static uint8_t COMP_LUT[256];
static void init_comp_lut() {
    for (int i = 0; i < 256; i++) COMP_LUT[i] = 'N';
    COMP_LUT[(uint8_t)'A'] = 'T'; COMP_LUT[(uint8_t)'T'] = 'A';
    COMP_LUT[(uint8_t)'C'] = 'G'; COMP_LUT[(uint8_t)'G'] = 'C';
}

// ---- NO REFERENCE ----
// The graph is built from reads alone. There is no reference index, so no strand vote and no repeat
// merge: every read enters BOTH ways (its bytes as sequenced and their reverse complement) and every
// kept core enters as extract() produced it. A locus therefore has a mirror twin in the other strand
// family; the twins are two nodes with different sigs and no core in common, and pairing them is the
// scorer's job, not the graph's.
// Repeats are handled by COVERAGE, not by a reference list (opt-in SVP_SPLIT_REPEATS, below).

// Node identity = core_sig (MurmurHash64A over the core's bytes), computed in extract(). The sig itself
// is what gets stored per node -- the core bytes are no longer kept.

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

// One parse of one strand: extract and drop the cores within `margin` of an end. The extract buffer
// stays valid in (n,st,en,sig) until the NEXT extract call, so a parse must be consumed before the
// other strand is parsed.
struct Parse { int n; const uint64_t *st, *en, *sig; };
Parse parse_strand(const char* data, uint64_t len, int level, int margin,
                   std::vector<int>& keep, uint64_t& cpu_ext) {
    Parse p;
    const uint64_t te = thr_ns();
    extract(data, (int)len, 0, level, &p.n, &p.st, &p.en, &p.sig);
    cpu_ext += thr_ns() - te;
    keep.clear();
    const uint64_t mg = (uint64_t)margin;
    for (int i = 0; i < p.n; i++) if (p.st[i] >= mg && p.en[i] + mg <= len) keep.push_back(i);
    return p;
}

// The kept cores of a parse, in read order, as extract() produced them.
void spans_of(const Parse& p, const std::vector<int>& keep, std::vector<Span>& out) {
    out.clear();
    out.reserve(keep.size());
    for (int i : keep) out.push_back({(uint32_t)p.st[i], (uint32_t)p.en[i], p.sig[i]});
}

void rc_into(const std::string& s, std::string& rc_buf) {
    const uint64_t len = s.size();
    rc_buf.resize(len);
    for (uint64_t j = 0; j < len; j++) rc_buf[j] = (char)COMP_LUT[(uint8_t)s[len - 1 - j]];
}

// One read -> interior signatures (NO edges; edges are emitted in pass 2 over the surviving-sig subset
// to avoid carrying 40+ GB of raw edges through the pool stage).
//
// The read goes in BOTH ways. Nothing here can tell which orientation a read was sequenced in without a
// reference, and the parse is not RC-symmetric (the two strands cut at different boundaries), so the two
// families are simply both kept.
//
// cap_st/cap_w (per strand): the pass-2 boundary cache. Caller guarantees len<=65535.
void process_read(std::string& s, int level, int margin,
                  std::vector<uint64_t>& lsigs, std::string& rc_buf, uint64_t& cpu_ext,
                  std::vector<uint16_t>* cfs, std::vector<uint16_t>* cfw,
                  std::vector<uint16_t>* crs, std::vector<uint16_t>* crw) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    std::vector<int> keep; std::vector<Span> eff;

    auto emit = [&](const Parse& p, std::vector<uint16_t>* cs, std::vector<uint16_t>* cw) {
        spans_of(p, keep, eff);
        for (const Span& e : eff) {
            if (cs) { cs->push_back((uint16_t)e.st); cw->push_back((uint16_t)(e.en - e.st)); }
            lsigs.push_back(e.sig);
        }
    };

    Parse p = parse_strand(s.data(), len, level, margin, keep, cpu_ext);
    emit(p, cfs, cfw);
    rc_into(s, rc_buf);
    p = parse_strand(rc_buf.data(), len, level, margin, keep, cpu_ext);
    emit(p, crs, crw);
}

// Per-edge (prev_sig, cur_sig) -> (normal_count, tumor_count). Hash mixes the sig pair; not tagged
// is_avalanching so unordered_dense re-mixes it (prev*golden+cur alone is not fully avalanched).
struct EdgeHash { uint64_t operator()(const SEdge& e) const { return e.first * 0x9e3779b97f4a7c15ULL + e.second; } };
using EdgeMap = ankerl::unordered_dense::map<SEdge, std::pair<uint32_t, uint32_t>, EdgeHash>;

// Pass 2: walk each read's surviving cores and COUNT edge traversals per colour. No byte-fill (payload
// dropped) -- only the sig identifies a node now. An adjacency (prev,cur) is counted once per
// occurrence, into the caller's per-colour slot (col 0 = normal -> .first, col 1 = tumor -> .second),
// for every consecutive pair where BOTH cores survived the mm gate (are in idx).
//
// Identical parse to pass 1 -- both strands, every kept core -- so the counted adjacencies match the
// node set.
// Repeat-run capture (SVP_SPLIT_REPEATS). A read physically threads  A -> X.. -> B  through a collapsed
// repeat run X. That pairing of the in-flank A to the out-flank B exists ONLY here, in the read: once the
// graph is built the repeat node has N in-flanks and N out-flanks and which goes with which is gone. So the
// thread is recorded now, and the X-cut (tools/trace.cpp apply_xcut) splits the run per (in,out) pair later.
void process_read_edges(std::string& s, int level, int margin, int col,
                        const ankerl::unordered_dense::map<uint64_t, uint32_t>& idx,
                        EdgeMap& emap, std::string& rc_buf, uint64_t& cpu_ext,
                        const uint8_t* is_repeat, std::vector<std::array<uint32_t, 3>>* trip) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    std::vector<int> keep; std::vector<Span> eff;

    auto chain = [&](const Parse& p) {
        spans_of(p, keep, eff);
        uint64_t prev = 0; bool have = false;
        // repeat-run state: the node id of the last non-repeat core seen (the in-flank), and the run itself
        uint32_t in_flank = UINT32_MAX, run_mid = UINT32_MAX; bool in_run = false;
        for (const Span& e : eff) {
            auto it = idx.find(e.sig);
            if (it == idx.end()) {                                         // core died at the mm gate -> chain breaks
                if (in_run && trip) { trip->push_back({run_mid, in_flank, UINT32_MAX}); in_run = false; }
                have = false; in_flank = UINT32_MAX; continue;
            }
            const uint32_t id = it->second;
            if (have && prev != e.sig) { auto& c = emap[SEdge(prev, e.sig)]; (col == 0 ? c.first : c.second)++; }
            prev = e.sig; have = true;
            if (is_repeat) {
                if (is_repeat[id]) {
                    if (!in_run) { in_run = true; run_mid = id; }          // run opens; anchor it at its first core
                } else {
                    if (in_run) { if (trip) trip->push_back({run_mid, in_flank, id}); in_run = false; }  // A -> X.. -> B
                    in_flank = id;
                }
            }
        }
        if (in_run && trip) trip->push_back({run_mid, in_flank, UINT32_MAX});   // read ended inside the run
    };

    Parse p = parse_strand(s.data(), len, level, margin, keep, cpu_ext);
    chain(p);
    rc_into(s, rc_buf);
    p = parse_strand(rc_buf.data(), len, level, margin, keep, cpu_ext);
    chain(p);
}

static bool is_gz(const std::string& p) { return p.size() >= 3 && p.compare(p.size() - 3, 3, ".gz") == 0; }

// A sample input is a directory of FASTQ/.svb shards (read in parallel) or a single file.
// Compressed input is rejected, not decompressed: `.fastq.gz` still contains ".fastq", so without
// this check it would be fopen'd as plain text and parsed into silent garbage.
std::vector<std::string> collect_shards(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) { std::fprintf(stderr, "[FATAL] no input %s\n", path.c_str()); std::abort(); }
    if (is_gz(path)) { std::fprintf(stderr, "[FATAL] %s is compressed; svp reads uncompressed FASTQ only -- decompress it first\n", path.c_str()); std::abort(); }
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

struct Shard { std::string path; int col; };
struct Batch { std::vector<std::string> seqs; int col; uint32_t shard = 0; uint32_t rid0 = 0; };   // rid0 = global read index (within shard) of seqs[0]

// Plain, uncompressed FASTQ. Decompression is not part of this program: measured on this box,
// gunzipping the 80 GB whole-genome .gz costs ~10 min while an assembler reads the resulting
// 188 GB plain FASTQ, parses it AND builds its graph in ~5 min. Paying that per run -- twice,
// once per pass -- bought nothing, so the input contract is uncompressed FASTQ, full stop.
static FILE* open_fastq(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[FATAL] fopen %s (input must be uncompressed FASTQ)\n", path.c_str()); std::abort(); }
    return f;
}

// ---- .svb compact reads (2-bit packed; svp uses only the ACGT sequence, no qual/name) ------
// A sample is a dir of shard.<k>.svb; each file = 32B header then sequential records:
//   len_tag u32 (low 31 = length L, bit31 = RAW); payload = RAW ? L ASCII bytes verbatim :
//   ceil(L/4) packed bytes, 2 bits/base (A=00 C=01 G=10 T=11), base i in byte i/4, high 2 bits
//   first. fopen'd (a real seekable file, NOT a pipe), so records fread exact-size -- no 1 MB
//   block/carry juggling. Sequential -> pass 1 and pass 2 re-read in identical (shard,rid) order.
static bool is_svb(const std::string& p) { return p.size() >= 4 && p.compare(p.size() - 4, 4, ".svb") == 0; }

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
        return true;
    }
    const uint32_t nb = (L + 3) >> 2;
    static thread_local std::vector<uint8_t> pb;
    pb.resize(nb);
    if (nb && std::fread(pb.data(), 1, nb, f) != nb) { std::fprintf(stderr, "[FATAL] .svb packed truncated\n"); std::abort(); }
    static const char LUT[4] = { 'A', 'C', 'G', 'T' };
    for (uint32_t i = 0; i < L; i++) seq[i] = LUT[(pb[i >> 2] >> (6 - 2 * (i & 3))) & 3];
    return true;
}

// ---- cross-pass core-boundary cache --------------------------------------------------------
// Pass 1 records every read's filtered LCP core boundaries (uint16 st + uint16 width, fwd and
// rc streams) so pass 2 replays them instead of re-running extract (~4.5k core-s, 18.6%
// of the build). Key = (shard, rid=read index within its shard): deterministic across passes
// (same file -> same byte order -> same parse order), so the out-of-order worker queue is
// irrelevant. Per-worker append-only arenas (lock-free); a (shard,rid)->location index is
// scattered after pass 1. Reads > 65535 bp can't be uint16-encoded -> flagged reext (pass 2
// falls back to extract for them). ~22.8 GB held pass1->pass2 at ~5.69 G emissions.
struct CoreCache {
    bool enabled = false;
    int nw = 0, nshard = 0;
    std::vector<std::vector<uint16_t>> fs, fw, rs, rw;       // [wid] fwd-start, fwd-width, rc-start, rc-width
    struct Loc { uint32_t foff = 0, roff = 0; uint16_t fn = 0, rn = 0, wid = 0; uint8_t reext = 1; };
    struct Rec { uint32_t shard, rid, foff, roff; uint16_t fn, rn; uint8_t reext; };
    std::vector<std::vector<Rec>> recs;                     // [wid] per-read records (transient, pass 1)
    std::vector<std::vector<Loc>> loc;                      // [shard][rid] -> arena location
    std::vector<uint64_t> reads_per_shard;                  // set by pass-1 readers
    void init(int nworkers, int nshards) {
        enabled = true; nw = nworkers; nshard = nshards;
        fs.resize(nw); fw.resize(nw); rs.resize(nw); rw.resize(nw); recs.resize(nw);
        loc.resize(nshard); reads_per_shard.assign(nshard, 0);
    }
    void finalize() {                                        // scatter per-worker recs -> loc[shard][rid]; free recs
        if (!enabled) return;
        for (int s = 0; s < nshard; s++) loc[s].assign(reads_per_shard[s], Loc{});
        for (int w = 0; w < nw; w++) for (const Rec& r : recs[w])
            loc[r.shard][r.rid] = Loc{r.foff, r.roff, r.fn, r.rn, (uint16_t)w, r.reext};
        std::vector<std::vector<Rec>>().swap(recs);
    }
    uint64_t bytes() const {                                 // arena footprint (the 22.8 GB)
        uint64_t b = 0;
        for (int w = 0; w < nw; w++) b += (fs[w].size() + fw[w].size() + rs[w].size() + rw[w].size()) * sizeof(uint16_t);
        for (int s = 0; s < nshard; s++) b += loc[s].size() * sizeof(Loc);
        return b;
    }
};

// ---- bloom cascade counting (replaces the exact per-shard hashmap) ----
// Per shard: (min_mult-1) plain bloom filters + an exact survivor map. Each occurrence of
// a sig advances it one filter deeper; once present in all (min_mult-1) filters, every
// further occurrence increments surv[sig]. So surv holds exactly the sigs seen >= min_mult
// times, with count = surv[sig] + (min_mult-1) (the first min_mult-1 hits are absorbed by
// the filters). Bloom collisions can promote a sig early (false survivor); they never drop
// one (FN=0). A shard is advanced only under its gmtx lock, so within a shard the cascade
// is serial and counts are exact -- the only deviation from the exact pipeline is bloom FP.
constexpr int BLOOM_HASHES = 4;
constexpr size_t BLOOM_BYTES_PER_FILTER = 4u << 20;   // 4 MB/shard/level (~19 bits/elem at ~110M distinct / 64 shards)

inline void bloom_hashes(uint64_t x, uint64_t* h) {
    for (int i = 0; i < BLOOM_HASHES; i++) {           // splitmix64 -> BLOOM_HASHES independent hashes from one sig
        x += 0x9e3779b97f4a7c15ULL;
        uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        h[i] = z ^ (z >> 31);
    }
}
// Plain (non-atomic) bloom. The cascade runs ONLY under the per-shard gmtx lock, so btllib's
// atomic fetch_or per bit was pure overhead (a lock-prefixed RMW + full barrier, ~15-30x a plain
// |=). Bit layout matches btllib exactly (bytes rounded up to 8; bit 1<<(h%8) at byte h/8) so the
// FP-bounded cascade stays bit-identical to the old path.
struct PlainBloom {
    std::vector<uint8_t> a;
    uint64_t bits = 0;
    void init(size_t bytes) { const size_t b = ((bytes + 7) / 8) * 8; a.assign(b, 0); bits = (uint64_t)b * 8; }
    bool contains_insert(const uint64_t* h) {            // true iff all bits already set; sets all bits
        uint8_t found = 1;
        for (int i = 0; i < BLOOM_HASHES; i++) {
            const uint64_t n = h[i] % bits; uint8_t& byte = a[n >> 3];
            found &= (uint8_t)((byte >> (n & 7)) & 1u); byte |= (uint8_t)(1u << (n & 7));
        }
        return (bool)found;
    }
};
struct ShardCascade {
    std::vector<PlainBloom> filt;
    ankerl::unordered_dense::map<uint64_t, uint32_t> surv;
    void init(int depth) { filt.resize(depth); for (auto& f : filt) f.init(BLOOM_BYTES_PER_FILTER); }
    void add(uint64_t s, const uint64_t* h, uint32_t cnt) {   // caller holds the shard lock
        auto it = surv.find(s);
        if (it != surv.end()) { it->second += cnt; return; }   // already promoted: skip the filter walk entirely
        for (uint32_t j = 0; j < cnt; j++) {
            bool promoted = true;
            for (auto& f : filt) if (!f.contains_insert(h)) { promoted = false; break; }
            if (promoted) { surv[s] += (cnt - j); break; }     // promoted now: absorb remaining occurrences, stop walking
        }
    }
};

// Both samples in ONE pool: shards_all (normal+tumor) -> 2 readers stream in parallel
// (overlapping the two input streams instead of reading the samples sequentially),
// feeding workers that route each colour-tagged batch into its own per-shard cascade
// (casc[0]=normal, casc[1]=tumor). Survivors are merged into one (n,t) record per sig.
// Affordable now only because the bloom cascades are RAM-bounded -- two exact hashmaps of
// all 103M distinct sigs would not coexist; the whole point of the bloom is to buy this.
void run_pool_both(const std::vector<Shard>& shards, int level, int margin, int threads,
                   int min_mult, std::vector<SigCount>& out_nodes, uint64_t& reads_n, uint64_t& reads_t,
                   CoreCache& cache, bool use_cc) {
    constexpr int NSHARD = 64;
    ShardCascade casc[2][NSHARD];                          // [0]=normal, [1]=tumor; one mutex each per colour
    for (int c = 0; c < 2; c++) for (int sh = 0; sh < NSHARD; sh++) casc[c][sh].init(min_mult - 1);
    std::mutex gmtx[2][NSHARD];
    const int nshards = (int)shards.size();
    const int nrd = nshards;                                       // one reader per file (cheap: fread + newline scan)
    const int nw  = std::max(1, threads - nrd);
    if (use_cc) cache.init(nw, nshards);          // pass-2 core-boundary cache (per-worker arenas filled below)
    std::atomic<size_t> next{0}, done{0};
    std::atomic<uint64_t> greads[2]{ {0}, {0} };
    std::atomic<int> readers_left{nrd};
    std::atomic<uint64_t> tot_rd_ns{0}, tot_w_busy_ns{0}, tot_w_idle_ns{0}, tot_q_full_ns{0};
    std::atomic<uint64_t> tot_count_ns{0}, tot_lock_ns{0};   // drain+flush (counting) vs gmtx acquire (contention)
    // CPU-seconds (CLOCK_THREAD_CPUTIME_ID, on-CPU only): reader parse, worker LCP (extract+sig),
    // the extract subset of LCP, and counting (drain+flush). Sum of worker LCP+count ~= worker on-CPU.
    std::atomic<uint64_t> cpu_rd_ns{0}, cpu_lcp_ns{0}, cpu_ext_ns{0}, cpu_cnt_ns{0}, cpu_w_ns{0};
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    constexpr size_t Q_CAP = 256, BATCH = 64;
    constexpr unsigned BLOCK = 1u << 20;
    auto reader = [&]() {
        const auto t_start = std::chrono::steady_clock::now();
        const uint64_t cpu0 = thr_ns();
        uint64_t local_q_full_ns = 0;
        std::vector<char> buf(BLOCK);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string carry; int col_local = 0; uint32_t cur_shard = 0; uint64_t shard_reads = 0;
        auto flush = [&]() {
            if (batch.empty()) return;
            const uint32_t r0 = (uint32_t)(shard_reads - batch.size());   // rid of batch[0] within this shard
            const auto t_w0 = std::chrono::steady_clock::now();
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            local_q_full_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_w0).count();
            q.push(Batch{std::move(batch), col_local, cur_shard, r0}); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        size_t si; std::string svseq;
        while ((si = next.fetch_add(1)) < shards.size()) {
            col_local = shards[si].col; cur_shard = (uint32_t)si; shard_reads = 0;
            const bool svb = is_svb(shards[si].path);
            FILE* f = svb ? open_svb(shards[si].path) : open_fastq(shards[si].path);
            if (svb) {
                while (svb_next(f, svseq)) {
                    shard_reads++; batch.emplace_back(std::move(svseq));
                    if (batch.size() >= BATCH) flush();
                }
                std::fclose(f);
            } else {
            int line_no = 0; carry.clear();
            int got;
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
            flush();                                 // never let a batch span shards (keeps rid/col single-shard)
            if (use_cc) cache.reads_per_shard[si] = shard_reads;
            std::fprintf(stderr, "[graph] shard %zu/%d done (col=%d)\n", (size_t)done.fetch_add(1) + 1, nshards, col_local);
        }
        flush();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
        tot_rd_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_start).count();
        tot_q_full_ns += local_q_full_ns;
        cpu_rd_ns += thr_ns() - cpu0;
    };
    auto worker = [&](int wid) {
        const uint64_t wcpu0 = thr_ns();
        uint64_t local_busy_ns = 0, local_idle_ns = 0, local_count_ns = 0, local_lock_ns = 0;
        uint64_t local_lcp_cpu = 0, local_ext_cpu = 0, local_cnt_cpu = 0;
        std::vector<CoreCache::Rec>* myrec = use_cc ? &cache.recs[wid] : nullptr;
        std::vector<uint16_t> *afs = nullptr, *afw = nullptr, *ars = nullptr, *arw = nullptr;
        if (use_cc) { afs = &cache.fs[wid]; afw = &cache.fw[wid]; ars = &cache.rs[wid]; arw = &cache.rw[wid]; }
        std::vector<uint64_t> lsig[2]; lsig[0].reserve(1u << 20); lsig[1].reserve(1u << 20);
        ankerl::unordered_dense::map<uint64_t, uint32_t> lmap[2]; lmap[0].reserve(1u << 20); lmap[1].reserve(1u << 20);
        std::string rc_buf; rc_buf.reserve(1u << 16);
        uint64_t lreads[2] = {0, 0};
        auto drain = [&](int c) { for (uint64_t g : lsig[c]) lmap[c][g]++; lsig[c].clear(); };
        auto flush = [&](int c) {
            drain(c);
            std::vector<std::vector<std::pair<uint64_t, uint32_t>>> by(NSHARD);
            for (const auto& kv : lmap[c]) by[kv.first & (NSHARD - 1)].push_back(kv);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (by[sh].empty()) continue;
                const auto tl0 = std::chrono::steady_clock::now();
                std::lock_guard<std::mutex> lk(gmtx[c][sh]);
                local_lock_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tl0).count();
                for (const auto& kv : by[sh]) { uint64_t h[BLOOM_HASHES]; bloom_hashes(kv.first, h); casc[c][sh].add(kv.first, h, kv.second); }
            }
            lmap[c].clear();
        };
        while (true) {
            Batch b;
            {
                const auto t_i0 = std::chrono::steady_clock::now();
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || readers_left.load() == 0; });
                local_idle_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_i0).count();
                if (q.empty()) break;
                b = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            const int c = b.col & 1;
            const auto t_b0 = std::chrono::steady_clock::now();
            const uint64_t lc0 = thr_ns();
            for (size_t j = 0; j < b.seqs.size(); j++) {
                std::string& s = b.seqs[j];
                const bool docap = use_cc && s.size() <= 65535;     // uint16 positions need len<=65535
                uint32_t foff = 0, roff = 0;
                if (docap) { foff = (uint32_t)afs->size(); roff = (uint32_t)ars->size(); }
                process_read(s, level, margin, lsig[c], rc_buf, local_ext_cpu,
                             docap ? afs : nullptr, docap ? afw : nullptr, docap ? ars : nullptr, docap ? arw : nullptr);
                lreads[c]++;
                if (myrec) {
                    CoreCache::Rec r; r.shard = b.shard; r.rid = b.rid0 + (uint32_t)j;
                    r.foff = foff; r.roff = roff;
                    r.fn = docap ? (uint16_t)(afs->size() - foff) : 0;
                    r.rn = docap ? (uint16_t)(ars->size() - roff) : 0;
                    r.reext = docap ? 0 : 1;
                    myrec->push_back(r);
                }
            }
            local_lcp_cpu += thr_ns() - lc0;
            local_busy_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_b0).count();
            const auto tc0 = std::chrono::steady_clock::now();
            const uint64_t cc0 = thr_ns();
            if (lsig[c].size() >= (1u << 20)) drain(c);
            if (lmap[c].size() >= (1u << 22)) flush(c);
            local_cnt_cpu += thr_ns() - cc0;
            local_count_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tc0).count();
        }
        const auto tcf = std::chrono::steady_clock::now();
        const uint64_t cf0 = thr_ns();
        flush(0); flush(1);
        local_cnt_cpu += thr_ns() - cf0;
        local_count_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tcf).count();
        greads[0] += lreads[0]; greads[1] += lreads[1];
        tot_w_busy_ns += local_busy_ns;
        tot_w_idle_ns += local_idle_ns;
        tot_count_ns += local_count_ns;
        tot_lock_ns += local_lock_ns;
        cpu_lcp_ns += local_lcp_cpu; cpu_ext_ns += local_ext_cpu; cpu_cnt_ns += local_cnt_cpu;
        cpu_w_ns += thr_ns() - wcpu0;
    };
    std::vector<std::thread> rds, wks;
    for (int i = 0; i < nrd; i++) rds.emplace_back(reader);
    for (int i = 0; i < nw;  i++) wks.emplace_back(worker, i);
    for (auto& t : rds) t.join();
    for (auto& t : wks) t.join();
    if (use_cc) {
        cache.finalize();
        std::fprintf(stderr, "[graph]   core-boundary cache: %.1f GB held for pass 2 (skips fill-pass extract)\n", cache.bytes() / 1e9);
    }
    std::fprintf(stderr, "[graph] pool(both): %d reader(s) + %d worker(s) on %d shard(s); reader-wall=%.1fs | per-worker-sum: LCP=%.1fs count(drain+flush)=%.1fs (gmtx lock-wait=%.1fs) idle=%.1fs\n",
                 nrd, nw, nshards, tot_rd_ns.load() / 1e9, tot_w_busy_ns.load() / 1e9, tot_count_ns.load() / 1e9, tot_lock_ns.load() / 1e9, tot_w_idle_ns.load() / 1e9);
    std::fprintf(stderr, "[graph] pool(both) CPU-sec (on-CPU, contention-immune): reader-parse=%.1f | worker-LCP=%.1f (of which extract=%.1f, sig/murmur=%.1f) | counting=%.1f | worker-total=%.1f\n",
                 cpu_rd_ns.load() / 1e9, cpu_lcp_ns.load() / 1e9, cpu_ext_ns.load() / 1e9, (cpu_lcp_ns.load() - cpu_ext_ns.load()) / 1e9, cpu_cnt_ns.load() / 1e9, cpu_w_ns.load() / 1e9);
    reads_n = greads[0].load(); reads_t = greads[1].load();
    // ---- merge the two cascades into one (n,t) record per surviving sig, per shard ----
    const uint32_t add = (uint32_t)(min_mult - 1);
    out_nodes.clear();
    const auto tm0 = std::chrono::steady_clock::now();
    for (int sh = 0; sh < NSHARD; sh++) {
        ankerl::unordered_dense::map<uint64_t, SigCount> m; m.reserve(casc[0][sh].surv.size() + casc[1][sh].surv.size());
        for (const auto& kv : casc[0][sh].surv) { SigCount& x = m[kv.first]; x.sig = kv.first; x.n = kv.second + add; }
        for (const auto& kv : casc[1][sh].surv) { SigCount& x = m[kv.first]; x.sig = kv.first; x.t = kv.second + add; }
        for (const auto& p : m) out_nodes.push_back(p.second);
    }
    std::sort(out_nodes.begin(), out_nodes.end(), [](const SigCount& a, const SigCount& b){ return a.sig < b.sig; });
    const double merge_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - tm0).count();
    std::fprintf(stderr, "[graph]   bloom cascade (both) -> %zu survive min_mult=%d (merge=%.1fs)\n", out_nodes.size(), min_mult, merge_s);
}

struct SCEdge { uint64_t a, b; uint32_t n, t; };   // sig-pair adjacency + per-colour counts (before node-id resolve)

// Pass-2 pool: re-LCP the same shards and COUNT per-colour edge traversals between surviving cores.
// No byte-fill, no payload -- a node is identified only by its sig now. Each worker accumulates a
// local (prev,cur)->(n,t) map, periodically merged into 64 sig-sharded global maps (one mutex each);
// after the join the shards are flattened into `out` (sig-pair edges + counts). Node-id resolution
// and the on-disk CEdge conversion happen in the caller.
void edge_count_pool(const std::vector<Shard>& shards, int level, int margin, int threads,
                     const ankerl::unordered_dense::map<uint64_t, uint32_t>& idx,
                     std::vector<SCEdge>& out,
                     const uint8_t* is_repeat = nullptr,
                     std::vector<std::array<uint32_t, 3>>* all_trip = nullptr, std::mutex* tmtx = nullptr) {
    constexpr int NSHARD = 64;
    EdgeMap gmap[NSHARD]; std::mutex gmtx[NSHARD];
    const int nshards = (int)shards.size();
    const int nrd = nshards;                                       // one reader per file
    const char* nwe = std::getenv("SVP_FILL_NW");
    const int nw  = nwe ? std::max(1, std::atoi(nwe)) : std::max(1, threads - nrd);
    std::atomic<size_t> next{0}, done{0};
    std::atomic<int> readers_left{nrd};
    std::atomic<uint64_t> cpu_rd_ns{0}, cpu_w_ns{0}, cpu_ext_ns{0};
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    constexpr size_t Q_CAP = 256, BATCH = 64;
    constexpr unsigned BLOCK = 1u << 20;
    auto reader = [&]() {
        const uint64_t cpu0 = thr_ns();
        std::vector<char> buf(BLOCK);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string carry; int col_local = 0; uint32_t cur_shard = 0; uint64_t shard_reads = 0;
        auto flush = [&]() {
            if (batch.empty()) return;
            const uint32_t r0 = (uint32_t)(shard_reads - batch.size());
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(Batch{std::move(batch), col_local, cur_shard, r0}); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        size_t si; std::string svseq;
        while ((si = next.fetch_add(1)) < shards.size()) {
            col_local = shards[si].col; cur_shard = (uint32_t)si; shard_reads = 0;
            const bool svb = is_svb(shards[si].path);
            FILE* f = svb ? open_svb(shards[si].path) : open_fastq(shards[si].path);
            if (svb) {
                while (svb_next(f, svseq)) {
                    shard_reads++; batch.emplace_back(std::move(svseq));
                    if (batch.size() >= BATCH) flush();
                }
                std::fclose(f);
            } else {
            int line_no = 0; carry.clear();
            int got;
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
            flush();                                 // shard boundary: keep batches single-shard
            std::fprintf(stderr, "[edge] shard %zu/%d done (col=%d)\n", (size_t)done.fetch_add(1) + 1, nshards, col_local);
        }
        flush();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
        cpu_rd_ns += thr_ns() - cpu0;
    };
    auto worker = [&]() {
        const uint64_t wcpu0 = thr_ns();
        uint64_t local_ext = 0;
        std::string rc_buf; rc_buf.reserve(1u << 16);
        EdgeMap local; local.reserve(1u << 20);
        constexpr size_t FLUSH = 1u << 22;                                     // merge local -> shards past ~4M distinct edges
        auto flush = [&]() {
            std::vector<std::vector<std::pair<SEdge, std::pair<uint32_t, uint32_t>>>> by(NSHARD);
            for (const auto& kv : local) by[EdgeHash{}(kv.first) & (NSHARD - 1)].push_back(kv);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (by[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(gmtx[sh]);
                for (const auto& kv : by[sh]) { auto& c = gmap[sh][kv.first]; c.first += kv.second.first; c.second += kv.second.second; }
            }
            local.clear();
        };
        std::vector<std::array<uint32_t, 3>> local_trip;
        while (true) {
            Batch b;
            {
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || readers_left.load() == 0; });
                if (q.empty()) break;
                b = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            const int c = b.col & 1;
            for (std::string& s : b.seqs) process_read_edges(s, level, margin, c, idx, local, rc_buf, local_ext, is_repeat, is_repeat ? &local_trip : nullptr);
            if (local.size() >= FLUSH) flush();
        }
        flush();
        if (is_repeat && all_trip && tmtx && !local_trip.empty()) {
            std::lock_guard<std::mutex> lk(*tmtx);
            all_trip->insert(all_trip->end(), local_trip.begin(), local_trip.end());
        }
        cpu_w_ns += thr_ns() - wcpu0; cpu_ext_ns += local_ext;
    };
    std::vector<std::thread> rds, wks;
    for (int i = 0; i < nrd; i++) rds.emplace_back(reader);
    for (int i = 0; i < nw;  i++) wks.emplace_back(worker);
    for (auto& t : rds) t.join();
    for (auto& t : wks) t.join();
    size_t tot = 0; for (int sh = 0; sh < NSHARD; sh++) tot += gmap[sh].size();
    out.clear(); out.reserve(tot);
    for (int sh = 0; sh < NSHARD; sh++) for (const auto& kv : gmap[sh]) out.push_back({kv.first.first, kv.first.second, kv.second.first, kv.second.second});
    std::fprintf(stderr, "[edge] pool CPU-sec (on-CPU): reader-parse=%.1f | worker-LCP extract=%.1f | worker-total=%.1f | distinct edges=%zu\n",
                 cpu_rd_ns.load() / 1e9, cpu_ext_ns.load() / 1e9, cpu_w_ns.load() / 1e9, tot);
}

// ============================================================================================
// SVP_POSAWARE: node identity is (sig, locus). Pass 1 (the bloom count above) is UNCHANGED and its
// survivor sigs are the survival gate. Placement runs ONCE here in pass 2 -- pass 1 never places, so
// there is no pass1/pass2 placement-consistency requirement; determinism holds within this pass alone.
// A repeat sig that pass 1 kept is split here into one node per reference copy the reads actually
// landed on; the read-order edge between two distant-locus nodes is the SV junction.
// ============================================================================================
struct SL { uint64_t sig, locus; };
inline bool operator==(const SL& a, const SL& b) { return a.sig == b.sig && a.locus == b.locus; }
struct SLHash { uint64_t operator()(const SL& k) const { return k.sig * 0x9e3779b97f4a7c15ULL + k.locus; } };
struct SLP { SL a, b; };
inline bool operator==(const SLP& x, const SLP& y) { return x.a == y.a && x.b == y.b; }
struct SLPHash { uint64_t operator()(const SLP& k) const {
    uint64_t h = k.a.sig * 0x9e3779b97f4a7c15ULL + k.a.locus;
    return (h ^ (k.b.sig * 0x100000001b3ULL)) + k.b.locus; } };
using NodeMap = ankerl::unordered_dense::map<SL, std::pair<uint32_t, uint32_t>, SLHash>;
using PairMap = ankerl::unordered_dense::map<SLP, std::pair<uint32_t, uint32_t>, SLPHash>;
struct PANode { uint64_t sig, locus; uint32_t n, t; };
struct PAEdge { uint64_t sa, la, sb, lb; uint32_t n, t; };

// One read -> place cores against the posmap, count per-(sig,locus) node occurrences and per-colour
// (sig,locus)->(sig,locus) edge traversals between surviving cores. Same parse as pass 1: both strands,
// raw cores.
void process_read_edges_pa(std::string& s, int level, int margin, int col,
                           const ankerl::unordered_dense::map<uint64_t, uint32_t>& surv,
                           const posmap::Map& pm, NodeMap& nmap, PairMap& emap,
                           std::string& rc_buf, uint64_t& cpu_ext) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    std::vector<int> keep; std::vector<Span> eff; std::vector<uint64_t> locus;
    auto go = [&](const Parse& p) {
        spans_of(p, keep, eff);
        place(eff, pm, locus);
        bool have = false; SL prev{0, 0};
        for (size_t i = 0; i < eff.size(); i++) {
            if (surv.find(eff[i].sig) == surv.end()) { have = false; continue; }   // survival gate (sig cleared mm in pass 1)
            const SL cur{eff[i].sig, locus[i]};
            { auto& c = nmap[cur]; (col == 0 ? c.first : c.second)++; }
            if (have && !(prev == cur)) { auto& c = emap[SLP{prev, cur}]; (col == 0 ? c.first : c.second)++; }
            prev = cur; have = true;
        }
    };
    Parse p = parse_strand(s.data(), len, level, margin, keep, cpu_ext);
    go(p);
    rc_into(s, rc_buf);
    p = parse_strand(rc_buf.data(), len, level, margin, keep, cpu_ext);
    go(p);
}

// Pass-2 pool for posaware: same reader/worker skeleton as edge_count_pool, but each worker accumulates
// (sig,locus) node counts AND (sig,locus)-pair edge counts, flushed into 64 sharded global maps.
void edge_count_pool_pa(const std::vector<Shard>& shards, int level, int margin, int threads,
                        const ankerl::unordered_dense::map<uint64_t, uint32_t>& surv, const posmap::Map& pm,
                        std::vector<PANode>& out_nodes, std::vector<PAEdge>& out_edges) {
    constexpr int NSHARD = 64;
    NodeMap gnode[NSHARD]; std::mutex gnmtx[NSHARD];
    PairMap gedge[NSHARD]; std::mutex gemtx[NSHARD];
    const int nshards = (int)shards.size();
    const int nrd = nshards;
    const char* nwe = std::getenv("SVP_FILL_NW");
    const int nw  = nwe ? std::max(1, std::atoi(nwe)) : std::max(1, threads - nrd);
    std::atomic<size_t> next{0}, done{0};
    std::atomic<int> readers_left{nrd};
    std::atomic<uint64_t> cpu_rd_ns{0}, cpu_w_ns{0}, cpu_ext_ns{0};
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    constexpr size_t Q_CAP = 256, BATCH = 64;
    constexpr unsigned BLOCK = 1u << 20;
    auto reader = [&]() {
        const uint64_t cpu0 = thr_ns();
        std::vector<char> buf(BLOCK);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string carry; int col_local = 0; uint32_t cur_shard = 0; uint64_t shard_reads = 0;
        auto flush = [&]() {
            if (batch.empty()) return;
            const uint32_t r0 = (uint32_t)(shard_reads - batch.size());
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(Batch{std::move(batch), col_local, cur_shard, r0}); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        size_t si; std::string svseq;
        while ((si = next.fetch_add(1)) < shards.size()) {
            col_local = shards[si].col; cur_shard = (uint32_t)si; shard_reads = 0;
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
            std::fprintf(stderr, "[edge-pa] shard %zu/%d done (col=%d)\n", (size_t)done.fetch_add(1) + 1, nshards, col_local);
        }
        flush();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
        cpu_rd_ns += thr_ns() - cpu0;
    };
    auto worker = [&]() {
        const uint64_t wcpu0 = thr_ns();
        uint64_t local_ext = 0;
        std::string rc_buf; rc_buf.reserve(1u << 16);
        NodeMap lnode; lnode.reserve(1u << 20);
        PairMap ledge; ledge.reserve(1u << 20);
        constexpr size_t FLUSH = 1u << 22;
        auto flush = [&]() {
            { std::vector<std::vector<std::pair<SL, std::pair<uint32_t, uint32_t>>>> by(NSHARD);
              for (const auto& kv : lnode) by[SLHash{}(kv.first) & (NSHARD - 1)].push_back(kv);
              for (int sh = 0; sh < NSHARD; sh++) { if (by[sh].empty()) continue; std::lock_guard<std::mutex> lk(gnmtx[sh]);
                  for (const auto& kv : by[sh]) { auto& c = gnode[sh][kv.first]; c.first += kv.second.first; c.second += kv.second.second; } }
              lnode.clear(); }
            { std::vector<std::vector<std::pair<SLP, std::pair<uint32_t, uint32_t>>>> by(NSHARD);
              for (const auto& kv : ledge) by[SLPHash{}(kv.first) & (NSHARD - 1)].push_back(kv);
              for (int sh = 0; sh < NSHARD; sh++) { if (by[sh].empty()) continue; std::lock_guard<std::mutex> lk(gemtx[sh]);
                  for (const auto& kv : by[sh]) { auto& c = gedge[sh][kv.first]; c.first += kv.second.first; c.second += kv.second.second; } }
              ledge.clear(); }
        };
        while (true) {
            Batch b;
            { std::unique_lock<std::mutex> lk(qmtx);
              cv_empty.wait(lk, [&]{ return !q.empty() || readers_left.load() == 0; });
              if (q.empty()) break;
              b = std::move(q.front()); q.pop(); cv_full.notify_one(); }
            const int c = b.col & 1;
            for (std::string& s : b.seqs) process_read_edges_pa(s, level, margin, c, surv, pm, lnode, ledge, rc_buf, local_ext);
            if (lnode.size() + ledge.size() >= FLUSH) flush();
        }
        flush();
        cpu_w_ns += thr_ns() - wcpu0; cpu_ext_ns += local_ext;
    };
    std::vector<std::thread> rds, wks;
    for (int i = 0; i < nrd; i++) rds.emplace_back(reader);
    for (int i = 0; i < nw;  i++) wks.emplace_back(worker);
    for (auto& t : rds) t.join();
    for (auto& t : wks) t.join();
    size_t tn = 0, te = 0;
    for (int sh = 0; sh < NSHARD; sh++) { tn += gnode[sh].size(); te += gedge[sh].size(); }
    out_nodes.clear(); out_nodes.reserve(tn);
    out_edges.clear(); out_edges.reserve(te);
    for (int sh = 0; sh < NSHARD; sh++) for (const auto& kv : gnode[sh]) out_nodes.push_back({kv.first.sig, kv.first.locus, kv.second.first, kv.second.second});
    for (int sh = 0; sh < NSHARD; sh++) for (const auto& kv : gedge[sh]) out_edges.push_back({kv.first.a.sig, kv.first.a.locus, kv.first.b.sig, kv.first.b.locus, kv.second.first, kv.second.second});
    std::fprintf(stderr, "[edge-pa] pool CPU-sec: reader=%.1f | worker=%.1f | distinct nodes=%zu edges=%zu\n",
                 cpu_rd_ns.load() / 1e9, cpu_w_ns.load() / 1e9, tn, te);
}

}  // namespace

// ---- merged stage: reads -> bloom-count -> edge-count -> self-contained graph.bin ----
// graph.bin v3 schema: u64 K, u64 sig[K], u32 ncov[K], u32 tcov[K], u64 E, CEdge[E].
//   sig[K]  = MurmurHash64A node identity (replaces the old raw-byte payload; colour is derived on
//             demand from ncov/tcov and no longer stored).
//   CEdge   = {u32 src, u32 dst, u32 ncov, u32 tcov}: directed adjacency with per-colour traversal
//             counts (how many normal / tumor reads cross it). No core sequence is stored anymore.
void run_graph_reads(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir,
                     int level, int margin, int min_mult, int threads) {
    init_comp_lut();
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
    struct Snap { double wall, self; };
    const auto torig = std::chrono::steady_clock::now();
    auto snap = [&]() -> Snap { return { std::chrono::duration<double>(std::chrono::steady_clock::now() - torig).count(), self_cpu_s() }; };
    const Snap P0 = snap();

    // ---- pass 1: bloom-cascade count BOTH samples in one parallel pool -> survivors with (n,t) ----
    // The cache records every read's core boundaries so pass 2 skips extract (SVP_NOCACHE disables).
    CoreCache cache; const bool use_cc = !std::getenv("SVP_NOCACHE");
    std::vector<SigCount> nodes; uint64_t nr_n = 0, nr_t = 0;
    run_pool_both(shards_all, level, margin, threads, min_mult, nodes, nr_n, nr_t, cache, use_cc);
    const uint64_t K = nodes.size();
    const Snap P1 = snap();
    tic.mark("pool_both");
    std::fprintf(stderr, "[graph]   survivors=%lu\n", (unsigned long)K);

    // ---- SVP_COUNT_ONLY: dump (sig,n,t) for every observed core and STOP (analysis mode; skips the
    // memory-heavy fill). counts.bin = u64 K + SigCount[K]. Lets us study what any mm-floor would drop. ----
    if (std::getenv("SVP_COUNT_ONLY")) {
        FILE* fc = xopen(dir + "/counts.bin.tmp", "wb");
        uint64_t Kw = K; wr(fc, &Kw, 8); wr(fc, nodes.data(), K * sizeof(SigCount));
        atomic_close(fc, dir + "/counts.bin.tmp", dir + "/counts.bin");
        std::fprintf(stderr, "[graph] COUNT_ONLY: wrote %lu (sig,n,t) -> %s/counts.bin (skipped fill)\n", (unsigned long)K, dir.c_str());
        return;
    }

    // ---- id map + node sigs (node id = sig-sorted index; the sig IS the stored node identity now,
    // colour is derived from ncov/tcov on demand and no longer materialised) ----
    std::vector<uint64_t> sig(K); std::vector<uint32_t> ncov(K), tcov(K);
    ankerl::unordered_dense::map<uint64_t, uint32_t> idx; idx.reserve(K);
    for (uint64_t i = 0; i < K; i++) {
        const SigCount& c = nodes[i];
        idx.emplace(c.sig, (uint32_t)i);
        sig[i] = c.sig; ncov[i] = c.n; tcov[i] = c.t;
    }
    std::vector<SigCount>().swap(nodes);
    const Snap P2 = snap();

    // ---- SVP_POSAWARE: (sig,locus) node model. Place cores against the posmap in pass 2, key nodes by
    // (sig,locus), write graph.bin (byte-sig identity PRESERVED -> eval_sig/refine stay valid) + a
    // locus.bin sidecar (u64 K, u64 locus[K]). Wholly separate from the sig-only path below; unset =>
    // that path runs byte-for-byte as before. `idx` (sig->id, above) is the survival gate. ----
    if (std::getenv("SVP_POSAWARE")) {
        const char* pmp = std::getenv("SVP_POSMAP");
        posmap::Map pm;
        if (!pm.open(pmp ? pmp : "outputs/svp/refpos_l5_all.bin")) {
            std::fprintf(stderr, "[FATAL] posaware: cannot open posmap (set SVP_POSMAP)\n"); std::abort();
        }
        std::vector<PANode> pn; std::vector<PAEdge> pe;
        edge_count_pool_pa(shards_all, level, margin, threads, idx, pm, pn, pe);
        pm.close();
        std::fprintf(stderr, "[posaware] observed %zu (sig,locus) nodes, %zu (sig,locus) edges\n", pn.size(), pe.size());

        const char* cme = std::getenv("SVP_POSAWARE_COPYMIN");
        const uint32_t copymin = cme ? (uint32_t)std::atoi(cme) : 1;   // per-copy read floor (default keep all observed)
        std::sort(pn.begin(), pn.end(), [](const PANode& a, const PANode& b){ return a.sig != b.sig ? a.sig < b.sig : a.locus < b.locus; });
        std::vector<uint64_t> psig, ploc; std::vector<uint32_t> pncov, ptcov;
        ankerl::unordered_dense::map<SL, uint32_t, SLHash> nid; nid.reserve(pn.size());
        for (const PANode& x : pn) {
            if (x.n + x.t < copymin) continue;
            nid.emplace(SL{x.sig, x.locus}, (uint32_t)psig.size());
            psig.push_back(x.sig); ploc.push_back(x.locus); pncov.push_back(x.n); ptcov.push_back(x.t);
        }
        std::vector<PANode>().swap(pn);
        const uint64_t Kp = psig.size();
        std::vector<CEdge> ce; ce.reserve(pe.size());
        for (const PAEdge& e : pe) {
            auto a = nid.find(SL{e.sa, e.la}), b = nid.find(SL{e.sb, e.lb});
            if (a != nid.end() && b != nid.end() && a->second != b->second) ce.push_back({a->second, b->second, e.n, e.t});
        }
        std::vector<PAEdge>().swap(pe);
        std::sort(ce.begin(), ce.end(), [](const CEdge& x, const CEdge& y){ return x.src != y.src ? x.src < y.src : x.dst < y.dst; });

        const std::string gp = dir + "/graph.bin";
        FILE* fg = xopen(gp + ".tmp", "wb");
        uint64_t Kw = Kp; wr(fg, &Kw, 8);
        wr(fg, psig.data(), Kp * 8); wr(fg, pncov.data(), Kp * 4); wr(fg, ptcov.data(), Kp * 4);
        uint64_t E = ce.size(); wr(fg, &E, 8); wr(fg, ce.data(), E * sizeof(CEdge));
        atomic_close(fg, gp + ".tmp", gp);
        const std::string lp = dir + "/locus.bin";
        FILE* flx = xopen(lp + ".tmp", "wb");
        wr(flx, &Kw, 8); wr(flx, ploc.data(), Kp * 8);
        atomic_close(flx, lp + ".tmp", lp);

        uint64_t nnov = 0, nunr = 0;
        for (uint64_t i = 0; i < Kp; i++) { if (ploc[i] == posmap::NOVEL) nnov++; else if (ploc[i] == posmap::UNRESOLVED) nunr++; }
        const double kd = (double)(Kp ? Kp : 1);
        std::fprintf(stderr, "[posaware] nodes=%lu  placed=%lu (%.1f%%)  novel=%lu (%.1f%%)  unresolved-repeat=%lu (%.1f%%)  copymin=%u  edges=%lu -> %s/graph.bin + locus.bin\n",
                     (unsigned long)Kp, (unsigned long)(Kp - nnov - nunr), 100.0 * (Kp - nnov - nunr) / kd,
                     (unsigned long)nnov, 100.0 * nnov / kd, (unsigned long)nunr, 100.0 * nunr / kd,
                     copymin, (unsigned long)E, dir.c_str());
        return;
    }

    // ---- repeat-resolution capture (opt-in SVP_SPLIT_REPEATS) ------------------------------------
    // A collapsed repeat run keeps its nodes, and the read's in->run->out thread is recorded alongside;
    // tools/trace.cpp apply_xcut then splits the run into one copy per (in,out) pair, which removes the
    // tangle WITHOUT removing the sequence. A repeat is flagged by COVERAGE -- a run collapsed from C
    // loci carries ~C x the single-copy depth -- so no reference list is involved.
    const bool split_rep = std::getenv("SVP_SPLIT_REPEATS");
    std::vector<uint8_t> is_repeat;
    std::vector<std::array<uint32_t, 3>> triples; std::mutex tmtx;
    if (split_rep) {
        std::vector<uint32_t> zn, zt;
        for (uint64_t i = 0; i < K; i++) { if (ncov[i]) zn.push_back(ncov[i]); if (tcov[i]) zt.push_back(tcov[i]); }
        auto med = [](std::vector<uint32_t>& v) -> uint32_t { if (v.empty()) return 0; std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end()); return v[v.size() / 2]; };
        const uint32_t dn = med(zn), dt = med(zt);
        const char* rme = std::getenv("SVP_REPEAT_MULT");
        const double rmult = rme ? std::atof(rme) : 1.5;        // flag = coverage >= rmult x median single-copy depth
        const uint32_t tn = (uint32_t)(rmult * dn), tt = (uint32_t)(rmult * dt);
        is_repeat.assign(K, 0); uint64_t nf = 0;
        for (uint64_t i = 0; i < K; i++) if ((dn && ncov[i] >= tn) || (dt && tcov[i] >= tt)) { is_repeat[i] = 1; nf++; }
        std::fprintf(stderr, "[split] candidate-repeat nodes: %lu/%lu (%.2f%%) | single-copy depth n=%u t=%u; flag >= %.2fx -> ncov>=%u or tcov>=%u (SVP_REPEAT_MULT)\n",
                     (unsigned long)nf, (unsigned long)K, 100.0 * nf / (double)K, dn, dt, rmult, tn, tt);
    }

    // ---- pass 2: count per-colour edge traversals between surviving cores (no payload fill) ----
    std::vector<SCEdge> sedges;
    edge_count_pool(shards_all, level, margin, threads, idx, sedges,
                    split_rep ? is_repeat.data() : nullptr, split_rep ? &triples : nullptr, split_rep ? &tmtx : nullptr);
    const Snap P3 = snap();
    tic.mark("edge_count");
    if (split_rep) {
        ankerl::unordered_dense::map<uint64_t, ankerl::unordered_dense::map<uint64_t, uint32_t>> agg;
        for (const auto& t : triples) agg[t[0]][((uint64_t)t[1] << 32) | t[2]]++;
        std::vector<RepTriple> recs;
        for (const auto& m : agg) for (const auto& p : m.second) recs.push_back({(uint32_t)m.first, (uint32_t)(p.first >> 32), (uint32_t)p.first, p.second});
        FILE* ft = xopen(dir + "/threads.bin.tmp", "wb");
        uint64_t nr = recs.size(); wr(ft, &nr, 8); wr(ft, recs.data(), nr * sizeof(RepTriple));
        atomic_close(ft, dir + "/threads.bin.tmp", dir + "/threads.bin");
        std::fprintf(stderr, "[split] %lu raw in->run->out threads -> %lu distinct triples over %lu runs -> threads.bin\n",
                     (unsigned long)triples.size(), (unsigned long)nr, (unsigned long)agg.size());
        std::vector<std::array<uint32_t, 3>>().swap(triples);
    }

    // ---- resolve sig-pair edges -> node-id CEdge with counts (sig<->id is a bijection => distinct) ----
    std::vector<CEdge> ce; ce.reserve(sedges.size());
    for (const SCEdge& e : sedges) {
        auto a = idx.find(e.a), b = idx.find(e.b);
        if (a != idx.end() && b != idx.end() && a->second != b->second) ce.push_back({a->second, b->second, e.n, e.t});
    }
    std::vector<SCEdge>().swap(sedges);
    std::sort(ce.begin(), ce.end(), [](const CEdge& x, const CEdge& y){ return x.src != y.src ? x.src < y.src : x.dst < y.dst; });
    const Snap P4 = snap();
    tic.mark("resolve_edges");

    // ---- write self-contained graph.bin v3: u64 K, u64 sig[K], u32 ncov[K], u32 tcov[K], u64 E, CEdge[E] ----
    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen(gp + ".tmp", "wb");
    uint64_t Kw = K; wr(fg, &Kw, 8);
    wr(fg, sig.data(), K * 8); wr(fg, ncov.data(), K * 4); wr(fg, tcov.data(), K * 4);
    uint64_t E = ce.size(); wr(fg, &E, 8); wr(fg, ce.data(), E * sizeof(CEdge));
    atomic_close(fg, gp + ".tmp", gp);
    const Snap P5 = snap();
    tic.mark("write");
    std::fprintf(stderr, "[graph] normal_reads=%lu tumor_reads=%lu nodes=%lu edges=%lu -> %s/graph.bin (v3: sig + edge counts, no payload)\n",
                 (unsigned long)nr_n, (unsigned long)nr_t, (unsigned long)K, (unsigned long)E, dir.c_str());

    // ---- per-component CPU-second + execution-efficiency table -------------------------------
    // CPU = Δself (all in-process threads). wall = phase wall.
    // efficiency = CPU / (ncpu * wall) = fraction of the box's cores the phase actually kept busy
    // (1.0 = every core saturated for the whole phase; <1.0 = idle cores / serial sections).
    cpu_set_t cs; CPU_ZERO(&cs);   // efficiency denominator = cores ACTUALLY allotted (taskset/affinity), not the box's total
    const double ncpu = (sched_getaffinity(0, sizeof(cs), &cs) == 0 && CPU_COUNT(&cs) > 0) ? (double)CPU_COUNT(&cs) : (double)sysconf(_SC_NPROCESSORS_ONLN);
    auto row = [&](const char* name, const Snap& lo, const Snap& hi) {
        const double wall = hi.wall - lo.wall, cpu = hi.self - lo.self;
        std::fprintf(stderr, "[eff] %-16s wall=%8.1fs  CPU=%9.1f  eff=%5.1f%%\n",
                     name, wall, cpu, 100.0 * cpu / (ncpu * wall + 1e-9));
    };
    std::fprintf(stderr, "\n===== graph build CPU-seconds & execution efficiency (ncpu=%.0f) =====\n", ncpu);
    std::fprintf(stderr, "[eff] %-16s %9s  %9s  %6s\n", "phase", "wall_s", "CPU_s", "eff");
    row("count(pool_both)", P0, P1);
    row("sig+idx",          P1, P2);
    row("edge_count",       P2, P3);
    row("resolve_edges",    P3, P4);
    row("write",            P4, P5);
    row("TOTAL",            P0, P5);
    std::fprintf(stderr, "=====================================================================\n");
}
