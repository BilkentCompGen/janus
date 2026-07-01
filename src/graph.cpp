#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1   // expose popen/pclose + sched_getaffinity under strict -std=c++17
#endif
#include "graph.h"
#include "btllib/bloom_filter.hpp"   // vendored: plain (atomic, lock-free) BloomFilter
// extract is declared in svp.h (provided by src/extract.cpp over biolcp).

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
// self_cpu_s()/child_cpu_s(): kernel-accounted user+sys CPU for all in-process threads
// (SELF) and all reaped child processes i.e. rapidgzip decode (CHILDREN).
static inline uint64_t thr_ns() { struct timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec; }
static double self_cpu_s()  { struct rusage r; getrusage(RUSAGE_SELF, &r);     return r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1e6 + r.ru_stime.tv_sec + r.ru_stime.tv_usec / 1e6; }
static double child_cpu_s() { struct rusage r; getrusage(RUSAGE_CHILDREN, &r); return r.ru_utime.tv_sec + r.ru_utime.tv_usec / 1e6 + r.ru_stime.tv_sec + r.ru_stime.tv_usec / 1e6; }

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

// ---- repeat-core blacklist (opt-in SVP_BLACKLIST) ----
// Cores whose sig occurs >minfreq times in the reference assembly (built by tools/build_blacklist) are
// merged out of the build: a run of blacklisted cores is absorbed into its non-blacklisted flanks as
// two overlapping (A+X)->(X+B) cores, so A->X->B becomes A->B. The merge lives in blacklist.cpp and is
// SHARED with the anchor stage (svp.h decls), which MUST apply the identical merge to the reference.

// Node identity = core_sig (XXH3-64 over the core's bytes), computed in extract(); pass 2
// recomputes it from the cached core boundaries. The stored payload IS the core bytes — same
// region as the signature, no flank.

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

// One read -> interior signatures (NO edges; edges are emitted in pass 2 over
// the surviving-sig subset to avoid carrying 40+ GB of raw edges through the
// pool stage). DUAL extract: forward bytes via parse1 + RC bytes via the
// remap+reverse trick (≡ parse2). rc_buf is worker-local.
// cap_st/cap_w (per strand): when non-null, append each kept core's (st, width=en-st) as uint16
// for the pass-2 boundary cache. Caller guarantees len <= 65535 when capturing (so st,en fit).
void process_read(std::string& s, int level, int margin,
                  std::vector<uint64_t>& lsigs, std::string& rc_buf, uint64_t& cpu_ext,
                  std::vector<uint16_t>* cfs, std::vector<uint16_t>* cfw,
                  std::vector<uint16_t>* crs, std::vector<uint16_t>* crw) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    const uint64_t mg = (uint64_t)margin;

    std::vector<int> keep; std::vector<std::pair<uint32_t, uint32_t>> eff;   // blacklist (A+X)->(X+B) merge scratch
    auto run_pass = [&](const char* data, std::vector<uint16_t>* cap_st, std::vector<uint16_t>* cap_w) {
        int n; const uint64_t *st, *en, *sig;
        const uint64_t te = thr_ns(); extract(data, (int)len, 0, level, &n, &st, &en, &sig); cpu_ext += thr_ns() - te;
        if (blacklist_on()) {                           // merge blacklisted runs into their flanks; count the merged sigs
            keep.clear();
            for (int i = 0; i < n; i++) if (st[i] >= mg && en[i] + mg <= len) keep.push_back(i);
            blacklist_merge(st, en, sig, keep, 65535, eff);
            for (const auto& e : eff) lsigs.push_back(core_sig(data + e.first, (int)(e.second - e.first)));
            return;
        }
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            if (cap_st) { cap_st->push_back((uint16_t)st[i]); cap_w->push_back((uint16_t)(en[i] - st[i])); }
            lsigs.push_back(sig[i]);
        }
    };
    run_pass(s.data(), cfs, cfw);
    rc_buf.resize(len);
    for (uint64_t j = 0; j < len; j++) rc_buf[j] = (char)COMP_LUT[(uint8_t)s[len - 1 - j]];
    run_pass(rc_buf.data(), crs, crw);
}

// Pass 2: byte-fill surviving sigs + emit deduped within-pass adjacent-sig
// edges. `seen` is a worker-local hash of (prev, g) edge keys — an edge is
// emitted only the FIRST time its key arrives in this worker. With 30×
// coverage the same edge would otherwise be emitted 30 times across reads
// covering its locus, blowing up raw edge volume to 65 GB (the actual measured
// peak). Per-worker dedup keeps each unique edge once per worker; final radix
// + std::unique on the (much smaller) global vector dedupes across workers.
// Edge-key hash: `a * golden + b`. Collision rate ~2^-64 — negligible.
// When use_cache: replay the pass-1 boundaries (cfs/cfw fwd, crs/crw rc; already margin-filtered)
// instead of calling extract. The per-core fill+edge body is identical either way.
void process_read_fill(std::string& s, int level, int margin,
                       const ankerl::unordered_dense::map<uint64_t, uint32_t>& idx,
                       std::vector<std::vector<uint8_t>>& slot_buf,
                       std::vector<std::atomic<uint8_t>>& claimed,
                       std::vector<SEdge>& ledges,
                       ankerl::unordered_dense::set<uint64_t>& seen,
                       std::string& rc_buf, uint64_t& cpu_ext,
                       const uint16_t* cfs, const uint16_t* cfw, uint32_t cfn,
                       const uint16_t* crs, const uint16_t* crw, uint32_t crn, bool use_cache,
                       const uint8_t* is_repeat, std::vector<std::array<uint32_t, 3>>* cap) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    const uint64_t mg = (uint64_t)margin;
    std::vector<int> keep; std::vector<std::pair<uint32_t, uint32_t>> eff;   // blacklist (A+X)->(X+B) merge scratch

    auto run_pass = [&](const char* data, const uint16_t* cst, const uint16_t* cwd, uint32_t cn) {
        uint64_t prev = 0; bool have = false;
        uint32_t before = UINT32_MAX, entry = UINT32_MAX; bool in_run = false;   // X-cut: span a flagged repeat run
        auto handle = [&](uint64_t st, uint64_t en) {           // sig + stored payload are the core bytes [st,en)
            const int Lw = (int)(en - st);
            const uint64_t g = core_sig(data + st, Lw);
            auto it = idx.find(g);
            if (it == idx.end()) { have = false;
                if (cap && in_run) { cap->push_back({{entry, before, UINT32_MAX}}); in_run = false; }
                if (cap) before = UINT32_MAX;
                return; }
            const uint32_t id = it->second;
            uint8_t exp = 0;
            if (claimed[id].compare_exchange_strong(exp, 1))
                slot_buf[id].assign((const uint8_t*)(data + st), (const uint8_t*)(data + st) + Lw);
            if (have && prev != g) {
                const uint64_t key = prev * 0x9e3779b97f4a7c15ULL + g;
                if (seen.insert(key).second) ledges.push_back(SEdge(prev, g));
            }
            prev = g; have = true;
            if (cap) {
                if (is_repeat[id]) { if (!in_run) { entry = id; in_run = true; } }
                else if (in_run) { cap->push_back({{entry, before, id}}); in_run = false; before = id; }
                else before = id;
            }
        };
        if (use_cache) {
            for (uint32_t k = 0; k < cn; k++) handle((uint64_t)cst[k], (uint64_t)cst[k] + cwd[k]);
        } else {
            int n; const uint64_t *st, *en, *sig;
            const uint64_t te = thr_ns(); extract(data, (int)len, 0, level, &n, &st, &en, &sig); cpu_ext += thr_ns() - te;
            if (blacklist_on()) {                       // (A+X)->(X+B): emit the merged effective cores (overlapping, no zero-overlap edge)
                keep.clear();
                for (int i = 0; i < n; i++) if (st[i] >= mg && en[i] + mg <= len) keep.push_back(i);
                blacklist_merge(st, en, sig, keep, 65535, eff);
                for (const auto& e : eff) handle(e.first, e.second);
            } else {
                for (int i = 0; i < n; i++) { if (st[i] < mg || en[i] + mg > len) continue; handle(st[i], en[i]); }
            }
        }
        if (cap && in_run) cap->push_back({{entry, before, UINT32_MAX}});       // run reached read end -> out=NONE
    };
    run_pass(s.data(), cfs, cfw, cfn);
    rc_buf.resize(len);
    for (uint64_t j = 0; j < len; j++) rc_buf[j] = (char)COMP_LUT[(uint8_t)s[len - 1 - j]];
    run_pass(rc_buf.data(), crs, crw, crn);
}

// A sample input is a directory of FASTQ shards (read in parallel) or a single file.
std::vector<std::string> collect_shards(const std::string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) { std::fprintf(stderr, "[FATAL] no input %s\n", path.c_str()); std::abort(); }
    if (!S_ISDIR(st.st_mode)) return {path};
    std::vector<std::string> v;
    DIR* d = opendir(path.c_str());
    if (!d) { std::fprintf(stderr, "[FATAL] opendir %s\n", path.c_str()); std::abort(); }
    for (struct dirent* e = readdir(d); e; e = readdir(d)) {
        const std::string n = e->d_name;
        if (n.find(".fastq") != std::string::npos || n.find(".fq") != std::string::npos) v.push_back(path + "/" + n);
    }
    closedir(d);
    if (v.empty()) { std::fprintf(stderr, "[FATAL] no .fastq/.fq shards in %s\n", path.c_str()); std::abort(); }
    std::sort(v.begin(), v.end());
    return v;
}

struct Shard { std::string path; int col; };
struct Batch { std::vector<std::string> seqs; int col; uint32_t shard = 0; uint32_t rid0 = 0; };   // rid0 = global read index (within shard) of seqs[0]

// Parallel single-stream gzip decode via rapidgzip: decodes ONE .gz across dec_p cores,
// in-process at run time -- NO input preprocessing (zlib gzread is ~120 MB/s/core and pins
// the reader; rapidgzip is ~1.5 GB/s at ~4 cores on the same file). The reader then reads
// plain FASTQ bytes from the pipe exactly as it did from gzread. Path overridable via env.
static FILE* open_decode(const std::string& path, int dec_p) {
    const char* bin = std::getenv("SVP_RAPIDGZIP"); if (!bin) bin = "lib/rapidgzip/venv/bin/rapidgzip";   // installed by `make install-deps`; override via SVP_RAPIDGZIP
    const std::string cmd = std::string(bin) + " -P " + std::to_string(dec_p) + " -d -c '" + path + "'";
    FILE* f = popen(cmd.c_str(), "r");
    if (!f) { std::fprintf(stderr, "[FATAL] popen rapidgzip %s\n", path.c_str()); std::abort(); }
    return f;
}

// ---- cross-pass core-boundary cache --------------------------------------------------------
// Pass 1 records every read's filtered LCP core boundaries (uint16 st + uint16 width, fwd and
// rc streams) so pass 2 replays them instead of re-running extract (~4.5k core-s, 18.6%
// of the build). Key = (shard, rid=read index within its shard): deterministic across passes
// (same .gz -> same rapidgzip bytes -> same parse order), so the out-of-order worker queue is
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

// Both samples in ONE pool: shards_all (normal+tumor) -> 2 readers decompress in parallel
// (overlapping the two gzip streams instead of decompressing the samples sequentially),
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
    const int nrd = nshards;                                       // one pipe-reader per file (cheap; decode is in rapidgzip)
    const char* dpenv = std::getenv("SVP_DECP");
    const int dec_p = dpenv ? std::max(1, std::atoi(dpenv)) : std::max(2, (threads / 4) / std::max(1, nshards));  // rapidgzip threads/file (~threads/4 total decode)
    const int nw  = std::max(1, threads - nrd);   // keep full worker count; rapidgzip decode threads ride pipe backpressure
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
        size_t si;
        while ((si = next.fetch_add(1)) < shards.size()) {
            col_local = shards[si].col; cur_shard = (uint32_t)si; shard_reads = 0;
            FILE* f = open_decode(shards[si].path, dec_p);
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
            pclose(f);
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
    std::fprintf(stderr, "[graph] pool(both) CPU-sec (on-CPU, contention-immune): reader-parse=%.1f | worker-LCP=%.1f (of which extract=%.1f, sig/xxh3=%.1f) | counting=%.1f | worker-total=%.1f\n",
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

void radix_sort_edges(std::vector<SEdge>& edges);   // fwd decl (defined below); used for fill_pool's periodic edge compaction

// Pass-2 pool: re-LCPs the same shards, fills byte slots for surviving sigs,
// AND emits within-pass adjacent-sig edges (only between surviving sigs). All
// edges land in `all_edges`; the pool stage no longer carries 40+ GB of raw
// edges through the pool peak.
void fill_pool(const std::vector<Shard>& shards, int level, int margin, int threads,
               const ankerl::unordered_dense::map<uint64_t, uint32_t>& idx,
               std::vector<std::vector<uint8_t>>& slot_buf,
               std::vector<std::atomic<uint8_t>>& claimed,
               std::vector<SEdge>& all_edges, std::mutex& emtx, const CoreCache& cache,
               const uint8_t* is_repeat, std::vector<std::array<uint32_t, 3>>& all_trip, std::mutex& tmtx) {
    const int nshards = (int)shards.size();
    const int nrd = nshards;                                       // one pipe-reader per file (decode is in rapidgzip)
    const char* dpenv = std::getenv("SVP_DECP");
    const int dec_p = dpenv ? std::max(1, std::atoi(dpenv)) : std::max(2, (threads / 4) / std::max(1, nshards));
    // Fill is feed-bound (workers are cheap once extract is cached), so DON'T oversubscribe:
    // subtract the rapidgzip decode threads from the worker budget -> readers+decode get real cores
    // and drain the pipes at full speed (decode/parse floor ~130s vs the oversubscribed ~410s).
    const char* nwe = std::getenv("SVP_FILL_NW");
    const int nw  = nwe ? std::max(1, std::atoi(nwe)) : std::max(1, threads - nrd - nrd * dec_p);
    std::atomic<size_t> next{0}, done{0};
    std::atomic<int> readers_left{nrd};
    // CPU-seconds (on-CPU only): reader parse, worker total, the extract (LCP) subset, and the
    // periodic edge insert+radix-compaction. per-read fill (hash+idx+slot+seen) = total-ext-compact.
    std::atomic<uint64_t> cpu_rd_ns{0}, cpu_w_ns{0}, cpu_ext_ns{0}, cpu_comp_ns{0};
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
            const uint32_t r0 = (uint32_t)(shard_reads - batch.size());   // rid of batch[0] within shard (matches pass 1)
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(Batch{std::move(batch), col_local, cur_shard, r0}); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        size_t si;
        while ((si = next.fetch_add(1)) < shards.size()) {
            col_local = shards[si].col; cur_shard = (uint32_t)si; shard_reads = 0;
            FILE* f = open_decode(shards[si].path, dec_p);
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
            pclose(f);
            flush();                                 // shard boundary: keep batches single-shard (rid alignment with pass 1)
            std::fprintf(stderr, "[fill] shard %zu/%d done (col=%d)\n", (size_t)done.fetch_add(1) + 1, nshards, col_local);
        }
        flush();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
        cpu_rd_ns += thr_ns() - cpu0;
    };
    auto worker = [&]() {
        const uint64_t wcpu0 = thr_ns();
        uint64_t local_ext = 0, local_comp = 0;
        std::string rc_buf; rc_buf.reserve(1u << 16);
        std::vector<SEdge> ledges; ledges.reserve(1u << 22);
        std::vector<std::array<uint32_t, 3>> ltrip;                           // X-cut threading (only if is_repeat)
        ankerl::unordered_dense::set<uint64_t> seen; seen.reserve(1u << 22);   // 4 M-entry edge-key cache (~160 MB / worker; was 16M/640MB)
        constexpr size_t E_FLUSH = 1u << 22;     // 4M edges (~64 MB) flush trigger
        constexpr size_t SEEN_CAP = 1u << 22;    // cap seen at 4M; clear on overflow (more compaction, ~7GB less RAM)
        constexpr size_t EDGE_COMPACT = 1u << 28;   // collapse all_edges (sort+unique) past ~268M to bound fill_pool RAM
        while (true) {
            Batch b;
            {
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || readers_left.load() == 0; });
                if (q.empty()) break;
                b = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            const bool have_cache = cache.enabled && b.shard < cache.loc.size();
            for (size_t j = 0; j < b.seqs.size(); j++) {
                std::string& s = b.seqs[j];
                const uint16_t *cfs = nullptr, *cfw = nullptr, *crs = nullptr, *crw = nullptr;
                uint32_t cfn = 0, crn = 0; bool use_cache = false;
                if (have_cache) {
                    const uint32_t rid = b.rid0 + (uint32_t)j;
                    if (rid < cache.loc[b.shard].size()) {
                        const CoreCache::Loc& L = cache.loc[b.shard][rid];
                        if (!L.reext) {
                            use_cache = true;
                            cfs = cache.fs[L.wid].data() + L.foff; cfw = cache.fw[L.wid].data() + L.foff; cfn = L.fn;
                            crs = cache.rs[L.wid].data() + L.roff; crw = cache.rw[L.wid].data() + L.roff; crn = L.rn;
                        }
                    }
                }
                process_read_fill(s, level, margin, idx, slot_buf, claimed, ledges, seen, rc_buf, local_ext,
                                  cfs, cfw, cfn, crs, crw, crn, use_cache, is_repeat, is_repeat ? &ltrip : nullptr);
            }
            const uint64_t cm0 = thr_ns();
            if (ledges.size() >= E_FLUSH) {
                std::lock_guard<std::mutex> lk(emtx);
                all_edges.insert(all_edges.end(), ledges.begin(), ledges.end());
                ledges.clear();
                if (all_edges.size() >= EDGE_COMPACT) {   // collapse dupes in place so the global edge vector stays bounded
                    radix_sort_edges(all_edges);
                    all_edges.erase(std::unique(all_edges.begin(), all_edges.end()), all_edges.end());
                }
            }
            if (seen.size() >= SEEN_CAP) seen.clear();   // overflow-safe; loses some cross-batch dedup but bounds memory
            if (is_repeat && ltrip.size() >= (1u << 22)) { std::lock_guard<std::mutex> lk(tmtx); all_trip.insert(all_trip.end(), ltrip.begin(), ltrip.end()); ltrip.clear(); }
            local_comp += thr_ns() - cm0;
        }
        { std::lock_guard<std::mutex> lk(emtx); all_edges.insert(all_edges.end(), ledges.begin(), ledges.end()); }
        if (is_repeat && !ltrip.empty()) { std::lock_guard<std::mutex> lk(tmtx); all_trip.insert(all_trip.end(), ltrip.begin(), ltrip.end()); }
        cpu_w_ns += thr_ns() - wcpu0; cpu_ext_ns += local_ext; cpu_comp_ns += local_comp;
    };
    std::vector<std::thread> rds, wks;
    for (int i = 0; i < nrd; i++) rds.emplace_back(reader);
    for (int i = 0; i < nw;  i++) wks.emplace_back(worker);
    for (auto& t : rds) t.join();
    for (auto& t : wks) t.join();
    const double w = cpu_w_ns.load() / 1e9, e = cpu_ext_ns.load() / 1e9, cp = cpu_comp_ns.load() / 1e9;
    std::fprintf(stderr, "[fill] pool CPU-sec (on-CPU, contention-immune): reader-parse=%.1f | worker-LCP extract=%.1f | per-read fill(hash+idx+slot+seen)=%.1f | edge insert+compact=%.1f | worker-total=%.1f\n",
                 cpu_rd_ns.load() / 1e9, e, w - e - cp, cp, w);
}

// LSB radix sort over (a, b) treated as a 16-byte key: 8 passes on .second, then 8 on .first.
void radix_sort_edges(std::vector<SEdge>& edges) {
    const size_t M = edges.size();
    if (M == 0) return;
    std::vector<SEdge> tmp(M);
    SEdge *src = edges.data(), *dst = tmp.data();
    uint32_t hist[256];
    for (int key = 0; key < 2; key++) for (int pass = 0; pass < 8; pass++) {
        std::memset(hist, 0, sizeof(hist));
        const int sh = pass * 8;
        for (size_t i = 0; i < M; i++) { const uint64_t k = key == 0 ? src[i].second : src[i].first; hist[(k >> sh) & 0xFFu]++; }
        uint32_t acc = 0;
        for (int b = 0; b < 256; b++) { uint32_t c = hist[b]; hist[b] = acc; acc += c; }
        for (size_t i = 0; i < M; i++) {
            const uint64_t k = key == 0 ? src[i].second : src[i].first;
            const uint32_t bk = (uint32_t)((k >> sh) & 0xFFu);
            dst[hist[bk]++] = src[i];
        }
        SEdge *t = src; src = dst; dst = t;
    }
    if (src != edges.data()) std::memcpy((void*)edges.data(), (void*)src, M * sizeof(SEdge));
}

}  // namespace

// ---- merged stage: reads -> bloom-count -> stream-fill -> self-contained graph.bin ----
// graph.bin v2 schema: u64 K, u8 color[K], u32 ncov[K], u32 tcov[K], u64 E, IEdge[E],
//                      u64 total_payload, u16 seqlen[K], u8 payload[total_payload].
// The per-node core-byte payload lives inline in graph.bin (compact.cpp reads it directly).
// The payload is written by streaming each slot straight to disk and freeing it -- there is
// no second contiguous copy, which keeps the write-phase RAM peak low.
void run_graph_reads(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir,
                     int level, int margin, int min_mult, int threads) {
    init_comp_lut();
    if (const char* blp = std::getenv("SVP_BLACKLIST")) load_blacklist(blp);
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
    // (SELF = all in-process threads incl. readers/workers/main; CHILDREN = reaped rapidgzip
    // decode). Snapshots taken at each phase boundary; deltas attribute cost to that phase. ----
    struct Snap { double wall, self, child; };
    const auto torig = std::chrono::steady_clock::now();
    auto snap = [&]() -> Snap { return { std::chrono::duration<double>(std::chrono::steady_clock::now() - torig).count(), self_cpu_s(), child_cpu_s() }; };
    const Snap P0 = snap();

    // ---- pass 1: bloom-cascade count BOTH samples in one parallel pool -> survivors with (n,t) ----
    // The cache records every read's core boundaries so pass 2 skips extract (SVP_NOCACHE disables).
    CoreCache cache; const bool use_cc = !std::getenv("SVP_NOCACHE") && !blacklist_on();   // blacklist merge rewrites core boundaries -> both passes must re-extract+re-merge (no cache)
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

    // ---- colour + id map (node id = sig-sorted index; sigs themselves are NOT persisted) ----
    std::vector<uint8_t> color(K); std::vector<uint32_t> ncov(K), tcov(K);
    ankerl::unordered_dense::map<uint64_t, uint32_t> idx; idx.reserve(K);
    for (uint64_t i = 0; i < K; i++) {
        const SigCount& c = nodes[i];
        idx.emplace(c.sig, (uint32_t)i);
        color[i] = c.n && c.t ? 3 : c.n ? 1 : 2; ncov[i] = c.n; tcov[i] = c.t;
    }
    std::vector<SigCount>().swap(nodes);

    // ---- pass 2: fill each node's sequence once + emit adjacency (sig pairs) ----
    std::vector<std::vector<uint8_t>> slot_buf(K);
    std::vector<std::atomic<uint8_t>> claimed(K); for (auto& c : claimed) c.store(0);
    std::vector<SEdge> edges; std::mutex emtx;
    // ---- repeat-resolution capture (opt-in SVP_SPLIT_REPEATS): flag high-coverage (collapsed-
    // repeat) nodes by coverage proxy; the fill pass then records each flagged run's in->run->out
    // span to threads.bin, consumed by the X-cut in the simplify stage. ----
    const bool split_rep = std::getenv("SVP_SPLIT_REPEATS");
    std::vector<uint8_t> is_repeat;
    std::vector<std::array<uint32_t, 3>> triples; std::mutex tmtx;
    if (split_rep) {
        std::vector<uint32_t> zn, zt;
        for (uint64_t i = 0; i < K; i++) { if (ncov[i]) zn.push_back(ncov[i]); if (tcov[i]) zt.push_back(tcov[i]); }
        auto med = [](std::vector<uint32_t>& v) -> uint32_t { if (v.empty()) return 0; std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end()); return v[v.size() / 2]; };
        const uint32_t dn = med(zn), dt = med(zt);
        const char* rme = std::getenv("SVP_REPEAT_MULT");
        const double rmult = rme ? std::atof(rme) : 1.5;        // repeat flag = coverage >= rmult x median single-copy depth
        const uint32_t tn = (uint32_t)(rmult * dn), tt = (uint32_t)(rmult * dt);
        is_repeat.assign(K, 0); uint64_t nf = 0;
        for (uint64_t i = 0; i < K; i++) if ((dn && ncov[i] >= tn) || (dt && tcov[i] >= tt)) { is_repeat[i] = 1; nf++; }
        std::fprintf(stderr, "[split] candidate-repeat nodes: %lu/%lu (single-copy depth n=%u t=%u; flag >= %.2fx -> ncov>=%u or tcov>=%u; set SVP_REPEAT_MULT)\n",
                     (unsigned long)nf, (unsigned long)K, dn, dt, rmult, tn, tt);
    }
    const Snap P2 = snap();
    fill_pool(shards_all, level, margin, threads, idx, slot_buf, claimed, edges, emtx, cache,
              split_rep ? is_repeat.data() : nullptr, triples, tmtx);
    const Snap P3 = snap();
    tic.mark("fill_pool");
    if (split_rep) {
        ankerl::unordered_dense::map<uint64_t, ankerl::unordered_dense::map<uint64_t, uint32_t>> agg;
        for (const auto& t : triples) agg[t[0]][((uint64_t)t[1] << 32) | t[2]]++;
        std::vector<RepTriple> recs;
        for (const auto& m : agg) for (const auto& p : m.second) recs.push_back({(uint32_t)m.first, (uint32_t)(p.first >> 32), (uint32_t)p.first, p.second});
        FILE* ft = xopen(dir + "/threads.bin.tmp", "wb");
        uint64_t nr = recs.size(); wr(ft, &nr, 8); wr(ft, recs.data(), nr * sizeof(RepTriple));
        atomic_close(ft, dir + "/threads.bin.tmp", dir + "/threads.bin");
        std::fprintf(stderr, "[split] wrote %lu distinct in->run->out triples -> threads.bin\n", (unsigned long)nr);
    }

    // ---- resolve sig-edges -> node-id edges, sort/dedup ----
    std::vector<IEdge> ie; ie.reserve(edges.size());
    for (const SEdge& e : edges) {
        auto a = idx.find(e.first), b = idx.find(e.second);
        if (a != idx.end() && b != idx.end() && a->second != b->second) ie.push_back(IEdge(a->second, b->second));
    }
    std::vector<SEdge>().swap(edges);
    std::sort(ie.begin(), ie.end()); ie.erase(std::unique(ie.begin(), ie.end()), ie.end());
    const Snap P4 = snap();
    tic.mark("resolve_edges");

    // ---- seqlen + total; slots stay scattered until the streaming write ----
    std::vector<uint16_t> seqlen(K); uint64_t total = 0;
    for (uint64_t i = 0; i < K; i++) {
        if (slot_buf[i].empty()) { std::fprintf(stderr, "[FATAL] empty slot %lu (pass-2 missed a kept sig)\n", (unsigned long)i); std::abort(); }
        seqlen[i] = (uint16_t)slot_buf[i].size(); total += slot_buf[i].size();
    }

    // ---- write self-contained graph.bin v2; stream each slot to disk and free it (no 2nd payload copy) ----
    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen(gp + ".tmp", "wb");
    uint64_t Kw = K; wr(fg, &Kw, 8);
    wr(fg, color.data(), K); wr(fg, ncov.data(), K * 4); wr(fg, tcov.data(), K * 4);
    uint64_t E = ie.size(); wr(fg, &E, 8); wr(fg, ie.data(), E * sizeof(IEdge));
    wr(fg, &total, 8); wr(fg, seqlen.data(), K * 2);
    for (uint64_t i = 0; i < K; i++) { wr(fg, slot_buf[i].data(), slot_buf[i].size()); std::vector<uint8_t>().swap(slot_buf[i]); }
    atomic_close(fg, gp + ".tmp", gp);
    const Snap P5 = snap();
    tic.mark("write");
    std::fprintf(stderr, "[graph] normal_reads=%lu tumor_reads=%lu nodes=%lu edges=%lu payload=%.1f GB -> %s/graph.bin (self-contained v2)\n",
                 (unsigned long)nr_n, (unsigned long)nr_t, (unsigned long)K, (unsigned long)E, total / 1e9, dir.c_str());
    if (blacklist_on()) std::fprintf(stderr, "[graph] blacklist (A+X)->(X+B) merge: %lu over-long merged spans dropped (>65535 bp)\n", (unsigned long)blacklist_dropped());

    // ---- per-component CPU-second + execution-efficiency table -------------------------------
    // CPU = in-process (Δself, all threads) + decode (Δchild, rapidgzip). wall = phase wall.
    // efficiency = CPU / (ncpu * wall) = fraction of the box's cores the phase actually kept busy
    // (1.0 = every core saturated for the whole phase; <1.0 = idle cores / serial sections).
    cpu_set_t cs; CPU_ZERO(&cs);   // efficiency denominator = cores ACTUALLY allotted (taskset/affinity), not the box's total
    const double ncpu = (sched_getaffinity(0, sizeof(cs), &cs) == 0 && CPU_COUNT(&cs) > 0) ? (double)CPU_COUNT(&cs) : (double)sysconf(_SC_NPROCESSORS_ONLN);
    auto row = [&](const char* name, const Snap& lo, const Snap& hi) {
        const double wall = hi.wall - lo.wall, inproc = hi.self - lo.self, dec = hi.child - lo.child;
        const double cpu = inproc + dec;
        std::fprintf(stderr, "[eff] %-16s wall=%8.1fs  CPU=%9.1f  (in-proc=%8.1f  decode=%7.1f)  eff=%5.1f%%\n",
                     name, wall, cpu, inproc, dec, 100.0 * cpu / (ncpu * wall + 1e-9));
    };
    std::fprintf(stderr, "\n===== graph build CPU-seconds & execution efficiency (ncpu=%.0f) =====\n", ncpu);
    std::fprintf(stderr, "[eff] %-16s %9s  %9s  %9s %9s  %6s\n", "phase", "wall_s", "CPU_s", "in-proc", "decode", "eff");
    row("count(pool_both)", P0, P1);
    row("colour+idx",       P1, P2);
    row("fill(fill_pool)",  P2, P3);
    row("resolve_edges",    P3, P4);
    row("seqlen+write",     P4, P5);
    row("TOTAL",            P0, P5);
    std::fprintf(stderr, "=====================================================================\n");
}
