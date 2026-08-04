#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1   // expose sched_getaffinity under strict -std=c++17
#endif
#include "graph.h"
// extract is declared in svp.h (provided by src/extract.cpp over biolcp).

#include <sys/stat.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <sched.h>
#include <dirent.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <memory>
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
static bool init_comp_lut() {
    for (int i = 0; i < 256; i++) COMP_LUT[i] = 'N';
    COMP_LUT[(uint8_t)'A'] = 'T'; COMP_LUT[(uint8_t)'T'] = 'A';
    COMP_LUT[(uint8_t)'C'] = 'G'; COMP_LUT[(uint8_t)'G'] = 'C';
    return true;
}
// Initialised on load, not by whichever driver remembers to call it. It was a driver's job once, and
// a driver that forgot got an all-zero table: rc_into() then produced a buffer of NULs, extract() found
// no cores in it, and the reverse-complement strand contributed NOTHING -- no error, no warning, just
// half the graph's coverage missing.
static const bool COMP_LUT_READY = init_comp_lut();

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

// ---- the build: ONE pass over the reads, exact maps -----------------------------------------
// This was two passes with a bloom cascade in front: the cascade knew a signature's MULTIPLICITY but
// not its NEIGHBOURS, so the reads had to be walked again once survivors were known. An exact map
// knows both. What that costs is holding adjacencies for error cores too, since survivorship is not
// decided until the end -- but the pruning rule is exact, not approximate: an edge is kept iff both
// endpoints survive, which is what breaking the chain at a dead core amounted to.
//
// The cascade comment said two exact hashmaps of every distinct signature would not coexist. With
// homopolymer normalisation that is no longer true -- distinct is 78.7 M per colour, and n/t share
// one record -- so the RAM the bloom was buying is no longer worth its second pass.
//
// One further difference, and it is an improvement: the cascade reports a colour's count as 0 when
// that colour stayed under the mm floor (its hits were absorbed by the filters and never recorded).
// The exact map reports what was actually seen. A node at ncov=25/tcov=2 came back as tcov=0 before.
using SigNodeMap = ankerl::unordered_dense::map<uint64_t, std::pair<uint32_t, uint32_t>>;

// ---- SVP_CASCADE: an admission gate in front of the exact maps ------------------------------
// The paragraph above holds for HiFi and stops holding for ONT. Measured on the COLO829 R10 pair:
// 96.5% of core emissions are DISTINCT (HiFi: 4.1%), because a 217 bp core survives a 0.6%/base
// error rate only 27% of the time. At 248.5 Gbp that is 5.02 G emissions -> ~3.8 G distinct sigs,
// and the two maps at their measured 29.6 and 39.0 bytes per entry want ~295 GB. The box has 157.
//
// Every one of those bytes is spent on cores that the mm floor discards at the end. So put a
// cascade of Bloom filters in front and only let a sig into the exact map once it has been seen
// min_mult times. Bloom has no false negatives, so nothing real is refused; the false positives it
// does admit are counted EXACTLY in pass 2 and then dropped by the same exact mm gate as before.
// The graph that comes out is therefore identical to the single-pass one -- verified on `hard`.
//
// This is not the old cascade. That one WAS the counter, which is why it reported a colour's count
// as 0 when that colour stayed under the floor; here it only decides who gets counted.
namespace casc {

// A blocked Bloom filter: all k bits of an item live in one 64-byte line, so an insert costs one
// cache miss instead of k. At 5 G inserts that is the whole difference.
struct BBloom {
    uint64_t* w = nullptr;
    uint64_t nblk = 0;
    static constexpr int K = 6;                     // 6 x 9 bits of position fits in one mixed word
    void init(uint64_t bits) {
        nblk = std::max<uint64_t>(1, bits >> 9);
        w = (uint64_t*)std::calloc(nblk * 8, 8);
        if (!w) { std::fprintf(stderr, "[FATAL] cascade: cannot allocate %.2f GB\n", nblk * 64 / 1e9); std::abort(); }
    }
    void release() { std::free(w); w = nullptr; nblk = 0; }
    uint64_t bytes() const { return nblk * 64; }
    static inline uint64_t mix(uint64_t x) {
        x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x;
    }
    // Sets the item's k bits and reports whether they were ALL already set (i.e. probably seen before).
    inline bool test_and_set(uint64_t h) {
        uint64_t* b = w + (((h >> 32) % nblk) << 3);
        uint64_t g = mix(h);
        bool all = true;
        for (int i = 0; i < K; i++) {
            const unsigned p = (unsigned)(g & 511); g >>= 9;
            const uint64_t m = 1ull << (p & 63);
            all &= (__atomic_fetch_or(&b[p >> 6], m, __ATOMIC_RELAXED) & m) != 0;
        }
        return all;
    }
};

constexpr int NSH = 64;                              // same sig-sharding as the node/edge globals
using SigSet = ankerl::unordered_dense::set<uint64_t>;

// min_mult-1 filters in series, then the exact set: filter i holds "seen at least i+1 times".
struct Gate {
    std::vector<BBloom> bf;
    SigSet adm[NSH];
    std::mutex mtx[NSH];
    // min_mult-2 filters, never fewer than one: admission at min_mult-1 sights, one short of the
    // floor. That slack is deliberate. Two threads can hit the same sig between one filter's load
    // and its store and both read it as unseen, which costs that sig a level; with the slack, any
    // core the exact gate would keep still gets in. Admitting a bit early is free -- pass 2 counts
    // it exactly and the mm gate throws it out -- whereas admitting late loses a real node.
    void init(int min_mult, uint64_t bits) {
        bf.resize(std::max(1, min_mult - 2));
        for (size_t i = 0; i < bf.size(); i++) bf[i].init(std::max<uint64_t>(1u << 26, bits >> (3 * i)));
    }
    void drop_filters() { for (auto& f : bf) f.release(); bf.clear(); }
    inline bool offer(uint64_t sig) {                // true => this sight admits the sig
        for (auto& f : bf) if (!f.test_and_set(sig)) return false;
        return true;
    }
    inline bool admitted(uint64_t sig) const { return adm[sig & (NSH - 1)].contains(sig); }
    uint64_t size() const { uint64_t n = 0; for (int i = 0; i < NSH; i++) n += adm[i].size(); return n; }
    uint64_t filter_bytes() const { uint64_t n = 0; for (const auto& f : bf) n += f.bytes(); return n; }
};

}  // namespace casc

// One read -> node counts AND raw consecutive adjacencies, both strands, in a single walk.
// `gate` non-null => pass 2 of the cascade: only admitted cores are counted, and the chain BREAKS at
// a non-admitted one rather than bridging over it -- which is exactly what "an edge is kept iff both
// endpoints survive" does in the single-pass path.
void process_read_single(std::string& s, int level, int margin, int col,
                         SigNodeMap& nmap, EdgeMap& emap, std::string& rc_buf, uint64_t& cpu_ext, uint64_t& ncores,
                         const casc::Gate* gate = nullptr) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    std::vector<int> keep; std::vector<Span> eff;
    auto chain = [&](const Parse& p) {
        spans_of(p, keep, eff);
        uint64_t prev = 0; bool have = false;
        ncores += eff.size();
        for (const Span& e : eff) {
            if (gate && !gate->admitted(e.sig)) { have = false; continue; }
            auto& c = nmap[e.sig]; (col == 0 ? c.first : c.second)++;
            // Every consecutive pair, survivor or not -- the mm gate is applied after the pass.
            if (have && prev != e.sig) { auto& x = emap[SEdge(prev, e.sig)]; (col == 0 ? x.first : x.second)++; }
            prev = e.sig; have = true;
        }
    };
    Parse p = parse_strand(s.data(), len, level, margin, keep, cpu_ext);
    chain(p);
    rc_into(s, rc_buf);
    p = parse_strand(rc_buf.data(), len, level, margin, keep, cpu_ext);
    chain(p);
}

// Cascade pass 1: the same walk with no maps at all -- every core is offered to the filters, and the
// sights that clear them go into a thread-local set for merging into the gate's shards.
void process_read_admit(std::string& s, int level, int margin, casc::Gate& g,
                        casc::SigSet& la, std::string& rc_buf, uint64_t& cpu_ext, uint64_t& ncores) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    std::vector<int> keep; std::vector<Span> eff;
    auto chain = [&](const Parse& p) {
        spans_of(p, keep, eff);
        ncores += eff.size();
        for (const Span& e : eff) if (g.offer(e.sig)) la.insert(e.sig);
    };
    Parse p = parse_strand(s.data(), len, level, margin, keep, cpu_ext);
    chain(p);
    rc_into(s, rc_buf);
    p = parse_strand(rc_buf.data(), len, level, margin, keep, cpu_ext);
    chain(p);
}

struct SCEdge { uint64_t a, b; uint32_t n, t; };   // sig-pair adjacency + per-colour counts (before node-id resolve)

// The single pass itself: readers stream both samples, workers accumulate exact local node/edge maps
// and periodically merge them into 64 sig-sharded globals. Same reader/worker/queue shape as the
// bloom pool, so the only variable being changed is what the workers accumulate INTO.
//
// `gate` + `admit=true` runs this same machinery as the cascade's pass 1: filters instead of maps,
// gate->adm as the only output, out_nodes/out_edges left untouched. `gate` + `admit=false` is pass 2.
void run_pool_single(const std::vector<Shard>& shards, int level, int margin, int threads, int min_mult,
                     std::vector<SigCount>& out_nodes, std::vector<SCEdge>& out_edges,
                     uint64_t& reads_n, uint64_t& reads_t,
                     casc::Gate* gate = nullptr, bool admit = false) {
    constexpr int NSHARD = 64;
    SigNodeMap gnode[NSHARD]; EdgeMap gedge[NSHARD];
    std::mutex nmtx[NSHARD], emtx[NSHARD];
    const int nshards = (int)shards.size();
    const int nrd = nshards;
    const int nw = std::max(1, threads - nrd);
    std::atomic<size_t> next{0}, done{0};
    std::atomic<int> readers_left{nrd};
    std::atomic<uint64_t> rn{0}, rt{0}, emissions{0};
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    constexpr size_t Q_CAP = 256, BATCH = 64;
    constexpr unsigned BLOCK = 1u << 20;

    auto reader = [&]() {
        std::vector<char> buf(BLOCK);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string carry; int col_local = 0; uint32_t cur_shard = 0; uint64_t shard_reads = 0;
        auto flush = [&]() {
            if (batch.empty()) return;
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(Batch{std::move(batch), col_local, cur_shard, 0}); cv_empty.notify_one();
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
            (col_local == 0 ? rn : rt) += shard_reads;
            std::fprintf(stderr, "[graph1] shard %zu/%d done (col=%d, %lu reads)\n",
                         (size_t)done.fetch_add(1) + 1, nshards, col_local, (unsigned long)shard_reads);
        }
        flush();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
    };

    auto worker = [&]() {
        std::string rc_buf; rc_buf.reserve(1u << 16);
        SigNodeMap ln; EdgeMap le; ln.reserve(1u << 20); le.reserve(1u << 20);
        casc::SigSet la;
        uint64_t ext = 0, emis = 0;
        const char* fe = std::getenv("SVP_FLUSH");   // test hook: force mid-run merges at smoke scale
        const size_t FLUSH = fe ? (size_t)std::atoll(fe) : (size_t)(1u << 22);
        auto flush_adm = [&]() {
            std::vector<std::vector<uint64_t>> ba(casc::NSH);
            for (uint64_t s : la) ba[s & (casc::NSH - 1)].push_back(s);
            for (int sh = 0; sh < casc::NSH; sh++) {
                if (ba[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(gate->mtx[sh]);
                for (uint64_t s : ba[sh]) gate->adm[sh].insert(s);
            }
            la.clear();
        };
        auto flush = [&]() {
            std::vector<std::vector<std::pair<uint64_t, std::pair<uint32_t, uint32_t>>>> bn(NSHARD);
            for (const auto& kv : ln) bn[kv.first & (NSHARD - 1)].push_back(kv);
            for (int sh = 0; sh < NSHARD; sh++) {
                if (bn[sh].empty()) continue;
                std::lock_guard<std::mutex> lk(nmtx[sh]);
                for (const auto& kv : bn[sh]) { auto& c = gnode[sh][kv.first]; c.first += kv.second.first; c.second += kv.second.second; }
            }
            ln.clear();
            std::vector<std::vector<std::pair<SEdge, std::pair<uint32_t, uint32_t>>>> be(NSHARD);
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
            if (admit) {
                for (std::string& s : b.seqs) process_read_admit(s, level, margin, *gate, la, rc_buf, ext, emis);
                if (la.size() >= FLUSH) flush_adm();
                continue;
            }
            for (std::string& s : b.seqs) process_read_single(s, level, margin, c, ln, le, rc_buf, ext, emis, gate);
            if (ln.size() >= FLUSH || le.size() >= FLUSH) flush();
        }
        admit ? flush_adm() : flush();
        emissions += emis;
    };

    { std::vector<std::thread> th;
      for (int i = 0; i < nrd; i++) th.emplace_back(reader);
      for (int i = 0; i < nw; i++) th.emplace_back(worker);
      for (auto& t : th) t.join(); }
    reads_n = rn; reads_t = rt;

    if (admit) {
        const uint64_t na = gate->size();
        std::fprintf(stderr, "[casc] pass1 emissions=%.3f M -> admitted=%.3f M (%.2f%% of emissions); "
                             "filters %.2f GB + admit set ~%.2f GB\n",
                     emissions.load() / 1e6, na / 1e6, emissions.load() ? 100.0 * na / emissions.load() : 0.0,
                     gate->filter_bytes() / 1e9, na * 21.0 / 1e9);
        return;
    }

    // ---- the mm gate, applied once, on exact counts ----
    uint64_t distinct = 0, kept = 0;
    ankerl::unordered_dense::set<uint64_t> live;
    for (int sh = 0; sh < NSHARD; sh++) distinct += gnode[sh].size();
    live.reserve(distinct);
    out_nodes.clear(); out_nodes.reserve(distinct / 2);
    for (int sh = 0; sh < NSHARD; sh++) {
        for (const auto& kv : gnode[sh]) {
            const uint32_t n = kv.second.first, t = kv.second.second;
            if (n < (uint32_t)min_mult && t < (uint32_t)min_mult) continue;   // survives in NEITHER colour
            out_nodes.push_back({kv.first, n, t}); live.insert(kv.first); kept++;
        }
        SigNodeMap().swap(gnode[sh]);                                            // free as we go
    }
    std::sort(out_nodes.begin(), out_nodes.end(), [](const SigCount& a, const SigCount& b){ return a.sig < b.sig; });

    uint64_t eraw = 0;
    for (int sh = 0; sh < NSHARD; sh++) eraw += gedge[sh].size();
    out_edges.clear(); out_edges.reserve(eraw / 3);
    for (int sh = 0; sh < NSHARD; sh++) {
        for (const auto& kv : gedge[sh]) {
            // Pass 2 breaks the chain at a dead core, so an edge exists iff BOTH ends survived.
            if (!live.count(kv.first.first) || !live.count(kv.first.second)) continue;
            out_edges.push_back({kv.first.first, kv.first.second, kv.second.first, kv.second.second});
        }
        EdgeMap().swap(gedge[sh]);
    }
    std::fprintf(stderr, "[graph1] core emissions=%.3f M\n", emissions.load() / 1e6);
    std::fprintf(stderr, "[graph1] distinct sigs=%.3f M -> kept=%.3f M (%.2f%%; discarded as read error %.3f M)\n",
                 distinct / 1e6, kept / 1e6, distinct ? 100.0 * kept / distinct : 0.0, (distinct - kept) / 1e6);
    std::fprintf(stderr, "[graph1] raw adjacencies=%.3f M -> both-ends-survive=%.3f M (%.2f%%)\n",
                 eraw / 1e6, out_edges.size() / 1e6, eraw ? 100.0 * out_edges.size() / eraw : 0.0);
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

    // ---- the build: ONE pass over the reads, exact maps, nodes and edges together ----
    // SVP_CASCADE puts an admission pass in front of it (see casc:: above). Sizing the filter needs
    // an item count before a single core has been parsed, and the input bytes give it: a FASTQ record
    // costs 2.04 bytes per base (measured on the ONT pair), each strand emits one core per `stride`
    // bases, and on a high-error read almost every emission is distinct -- so distinct ~ emissions.
    // 16 bits each is where a blocked filter's false-positive rate is a few parts in a thousand, and
    // a false positive costs nothing but an entry that pass 2 counts and the mm gate then discards.
    // Stride is level-dependent: 97.7 bases at L5, and each level down multiplies emissions by ~2.35
    // (measured 116.2 / 273.4 / 641.4 M for L5 / L4 / L3 over the same 5.75 Gbp of ONT reads).
    // SVP_CASCADE_GB caps the filter for a box with less headroom -- a smaller filter only means more
    // false admissions, never a lost node.
    std::vector<SigCount> nodes; std::vector<SCEdge> sedges; uint64_t nr_n = 0, nr_t = 0;
    std::unique_ptr<casc::Gate> gate;
    if (std::getenv("SVP_CASCADE")) {
        double bases = 0;
        for (const auto& s : shards_all) {
            struct stat st;
            if (::stat(s.path.c_str(), &st)) continue;
            bases += st.st_size * (is_svb(s.path) ? 4.0 : 1 / 2.04);   // .svb packs 4 bases a byte
        }
        const double stride = 97.7 / std::pow(2.35, 5 - level);
        uint64_t bits = std::max<uint64_t>(1u << 28, (uint64_t)(16.0 * 2.0 * bases / stride));
        if (const char* cg = std::getenv("SVP_CASCADE_GB")) bits = std::min(bits, (uint64_t)(std::atof(cg) * 8e9));
        gate = std::make_unique<casc::Gate>();
        gate->init(min_mult, bits);
        std::fprintf(stderr, "[casc] on: %d filter(s), %.2f GB, over %.1f Gbp of reads; admits at >=%d sights (mm=%d)\n",
                     (int)gate->bf.size(), gate->filter_bytes() / 1e9, bases / 1e9, std::max(1, min_mult - 2) + 1, min_mult);
        uint64_t d0 = 0, d1 = 0;
        run_pool_single(shards_all, level, margin, threads, min_mult, nodes, sedges, d0, d1, gate.get(), true);
        gate->drop_filters();                       // the filters have done their job; pass 2 only reads the set
    }
    run_pool_single(shards_all, level, margin, threads, min_mult, nodes, sedges, nr_n, nr_t, gate.get(), false);
    const uint64_t K = nodes.size();
    const Snap P1 = snap();
    tic.mark("pool_both");
    std::fprintf(stderr, "[graph]   survivors=%lu\n", (unsigned long)K);

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
    row("resolve_edges",    P2, P4);
    row("write",            P4, P5);
    row("TOTAL",            P0, P5);
    std::fprintf(stderr, "=====================================================================\n");
}
