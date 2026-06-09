#include "graph.h"
#include "MurmurHash3.h"

extern "C" void extract_fast(const char *str, int len, uint64_t offset, int target_level,
                             int *n_out, const uint64_t **starts_out,
                             const uint64_t **ends_out, const uint32_t **labels_out);

#include <zlib.h>
#include <sys/stat.h>
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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// ---- fail-fast binary checkpoint IO (declared in svp.h; reused by modea.cpp) ----
FILE* xopen(const std::string& p, const char* m) { FILE* f = std::fopen(p.c_str(), m); if (!f) { std::fprintf(stderr, "[FATAL] fopen %s\n", p.c_str()); std::abort(); } return f; }
void wr(FILE* f, const void* p, size_t n) { if (n && std::fwrite(p, 1, n, f) != n) { std::fprintf(stderr, "[FATAL] fwrite\n"); std::abort(); } }
void rd(FILE* f, void* p, size_t n) { if (n && std::fread(p, 1, n, f) != n) { std::fprintf(stderr, "[FATAL] truncated checkpoint\n"); std::abort(); } }
void atomic_close(FILE* f, const std::string& tmp, const std::string& dst) { std::fclose(f); if (std::rename(tmp.c_str(), dst.c_str())) { std::fprintf(stderr, "[FATAL] rename %s\n", dst.c_str()); std::abort(); } }   // never leave a half-written checkpoint resume would trust

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

// Per-pass sig: murmur64(W, 0). seed=0 (label seed redundant for flank≥0, per
// memory project_label_seed_redundant). The pass's byte buffer (R for pass1,
// rc_buf for pass2) IS the canonical representation for the pass's perspective;
// strand pooling happens ACROSS READS (a +strand read's pass1 sig at locus L
// equals a -strand read's pass2 sig at locus L because the bytes coincide).
// We deliberately do NOT canonicalize across {W, RC(W)} at the per-emission
// level because storing canonical bytes breaks Mode A's KMP suffix/prefix
// splice — adjacent chain nodes can end up stored in opposite orientations
// and lose their byte overlap. Strand-pair pooling at a higher level (super-
// node / bubble byte-RC matching) is a follow-up.
static inline uint64_t emit_sig(const char* W, int Lw) {
    return murmur64(W, Lw, 0);
}

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
void process_read(std::string& s, int level, int flank, int margin,
                  std::vector<uint64_t>& lsigs, std::string& rc_buf) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    const uint64_t mu = (uint64_t)flank, mg = (uint64_t)margin;

    auto run_pass = [&](const char* data) {
        int n; const uint64_t *st, *en; const uint32_t *lb;
        extract_fast(data, (int)len, 0, level, &n, &st, &en, &lb);
        (void)lb;
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            const uint64_t L = st[i] > mu ? st[i] - mu : 0;
            const uint64_t R = en[i] + mu < len ? en[i] + mu : len;
            const uint64_t g = emit_sig(data + L, (int)(R - L));
            lsigs.push_back(g);
        }
    };
    run_pass(s.data());
    rc_buf.resize(len);
    for (uint64_t j = 0; j < len; j++) rc_buf[j] = (char)COMP_LUT[(uint8_t)s[len - 1 - j]];
    run_pass(rc_buf.data());
}

// Pass 2: byte-fill surviving sigs + emit deduped within-pass adjacent-sig
// edges. `seen` is a worker-local hash of (prev, g) edge keys — an edge is
// emitted only the FIRST time its key arrives in this worker. With 30×
// coverage the same edge would otherwise be emitted 30 times across reads
// covering its locus, blowing up raw edge volume to 65 GB (the actual measured
// peak). Per-worker dedup keeps each unique edge once per worker; final radix
// + std::unique on the (much smaller) global vector dedupes across workers.
// Edge-key hash: `a * golden + b`. Collision rate ~2^-64 — negligible.
void process_read_fill(std::string& s, int level, int flank, int margin,
                       const std::unordered_map<uint64_t, uint32_t>& idx,
                       std::vector<std::vector<uint8_t>>& slot_buf,
                       std::vector<std::atomic<uint8_t>>& claimed,
                       std::vector<SEdge>& ledges,
                       std::unordered_set<uint64_t>& seen,
                       std::string& rc_buf) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    const uint64_t mu = (uint64_t)flank, mg = (uint64_t)margin;

    auto run_pass = [&](const char* data) {
        int n; const uint64_t *st, *en; const uint32_t *lb;
        extract_fast(data, (int)len, 0, level, &n, &st, &en, &lb);
        (void)lb;
        uint64_t prev = 0; bool have = false;
        for (int i = 0; i < n; i++) {
            if (st[i] < mg || en[i] + mg > len) continue;
            const uint64_t L = st[i] > mu ? st[i] - mu : 0;
            const uint64_t R = en[i] + mu < len ? en[i] + mu : len;
            const int Lw = (int)(R - L);
            const uint64_t g = emit_sig(data + L, Lw);
            auto it = idx.find(g);
            if (it == idx.end()) { have = false; continue; }
            uint8_t exp = 0;
            if (claimed[it->second].compare_exchange_strong(exp, 1))
                slot_buf[it->second].assign((const uint8_t*)(data + L), (const uint8_t*)(data + L) + Lw);
            if (have && prev != g) {
                const uint64_t key = prev * 0x9e3779b97f4a7c15ULL + g;
                if (seen.insert(key).second) ledges.push_back(SEdge(prev, g));
            }
            prev = g; have = true;
        }
    };
    run_pass(s.data());
    rc_buf.resize(len);
    for (uint64_t j = 0; j < len; j++) rc_buf[j] = (char)COMP_LUT[(uint8_t)s[len - 1 - j]];
    run_pass(rc_buf.data());
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
struct Batch { std::vector<std::string> seqs; int col; };

// Reader/worker split with both samples in one pool: 1 reader per shard (so both .fastq.gz
// inflate in parallel CPUs) + N workers consuming colour-tagged batches. The reader is a
// block-based FASTQ parser (1 MB gzread + manual \n scan) — line-by-line gzgets was 3-5x
// slower per reader thread (zlib state crossing + string ops dominated).
void run_pool(const std::vector<Shard>& shards, int level, int flank, int margin, int threads,
              std::vector<uint64_t>& all_sigs, std::mutex& smtx,
              uint64_t& reads_n_out, uint64_t& reads_t_out) {
    const int nshards = (int)shards.size();
    const int nrd = std::min(nshards, std::max(2, threads / 4));   // ≥2 so both samples decompress in parallel
    const int nw  = std::max(1, threads - nrd);
    std::atomic<size_t> next{0}, done{0};
    std::atomic<uint64_t> nr_n{0}, nr_t{0};
    std::atomic<int> readers_left{nrd};
    // Aggregate ns spent in each phase across all threads (sum, not max). Worker idle =
    // queue empty wait; worker busy = process_read. Reader wall = gzread + parse + flush.
    std::atomic<uint64_t> tot_rd_ns{0}, tot_w_busy_ns{0}, tot_w_idle_ns{0}, tot_q_full_ns{0};
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    constexpr size_t Q_CAP = 256, BATCH = 64;
    constexpr unsigned BLOCK = 1u << 20;
    auto reader = [&]() {
        const auto t_start = std::chrono::steady_clock::now();
        uint64_t local_q_full_ns = 0;
        std::vector<char> buf(BLOCK);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string carry; int col_local = 0;
        auto flush = [&]() {
            if (batch.empty()) return;
            const auto t_w0 = std::chrono::steady_clock::now();
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            local_q_full_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_w0).count();
            q.push(Batch{std::move(batch), col_local}); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        size_t si;
        while ((si = next.fetch_add(1)) < shards.size()) {
            col_local = shards[si].col;
            gzFile f = gzopen(shards[si].path.c_str(), "rb");
            if (!f) { std::fprintf(stderr, "[FATAL] gzopen %s\n", shards[si].path.c_str()); std::abort(); }
            gzbuffer(f, 1u << 20);
            int line_no = 0; carry.clear();
            int got;
            while ((got = gzread(f, buf.data(), BLOCK)) > 0) {
                const char* base = buf.data(); const char* end = base + got; const char* start = base;
                for (const char* p = base; p < end; p++) if (*p == '\n') {
                    const size_t n = (size_t)(p - start);
                    if (line_no == 1) {
                        if (carry.empty()) batch.emplace_back(start, n);
                        else { carry.append(start, n); batch.push_back(std::move(carry)); carry.clear(); }
                        if (batch.size() >= BATCH) flush();
                    } else if (!carry.empty()) carry.clear();
                    line_no = (line_no + 1) & 3; start = p + 1;
                }
                if (start < end && line_no == 1) carry.append(start, end - start);
            }
            gzclose(f);
            std::fprintf(stderr, "[lcp] shard %zu/%d done (col=%d)\n", (size_t)done.fetch_add(1) + 1, nshards, col_local);
        }
        flush();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
        tot_rd_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_start).count();
        tot_q_full_ns += local_q_full_ns;
    };
    auto worker = [&]() {
        uint64_t local_busy_ns = 0, local_idle_ns = 0;
        std::vector<uint64_t> lsigs; lsigs.reserve(1u << 21);
        std::string rc_buf; rc_buf.reserve(1u << 16);   // worker-local RC buffer reused across reads
        uint64_t lnr_n = 0, lnr_t = 0;
        constexpr size_t FLUSH_THRESH = 1u << 21;       // ~2M entries (~24 MB sigs)
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
            const auto t_b0 = std::chrono::steady_clock::now();
            for (auto& s : b.seqs) { process_read(s, level, flank, margin, lsigs, rc_buf); (b.col == 0 ? lnr_n : lnr_t)++; }
            local_busy_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_b0).count();
            if (lsigs.size() >= FLUSH_THRESH) {
                std::lock_guard<std::mutex> lk(smtx);
                all_sigs.insert(all_sigs.end(), lsigs.begin(), lsigs.end());
                lsigs.clear();
            }
        }
        { std::lock_guard<std::mutex> lk(smtx); all_sigs.insert(all_sigs.end(), lsigs.begin(), lsigs.end()); }
        nr_n += lnr_n; nr_t += lnr_t;
        tot_w_busy_ns += local_busy_ns;
        tot_w_idle_ns += local_idle_ns;
    };
    std::vector<std::thread> rds, wks;
    for (int i = 0; i < nrd; i++) rds.emplace_back(reader);
    for (int i = 0; i < nw;  i++) wks.emplace_back(worker);
    for (auto& t : rds) t.join();
    for (auto& t : wks) t.join();
    const double rd_s    = tot_rd_ns.load()     / 1e9;
    const double wb_s    = tot_w_busy_ns.load() / 1e9;
    const double wi_s    = tot_w_idle_ns.load() / 1e9;
    const double qf_s    = tot_q_full_ns.load() / 1e9;
    std::fprintf(stderr, "[lcp] pool: %d reader(s) + %d worker(s) on %d shard(s)\n", nrd, nw, nshards);
    std::fprintf(stderr, "[lcp]    reader-wall sum = %.1fs (avg/reader = %.1fs); queue-full wait = %.1fs\n",
                 rd_s, rd_s / nrd, qf_s);
    std::fprintf(stderr, "[lcp]    worker-busy sum = %.1fs (process_read total); worker-idle sum = %.1fs (queue-empty wait)\n",
                 wb_s, wi_s);
    std::fprintf(stderr, "[lcp]    worker utilization = %.1f%% (busy / (busy+idle)); reader vs workers: %s-bound\n",
                 100.0 * wb_s / (wb_s + wi_s + 1e-9),
                 wi_s > qf_s * 2 ? "reader(inflate)" : qf_s > wi_s * 2 ? "worker(LCP)" : "balanced");
    reads_n_out = nr_n.load(); reads_t_out = nr_t.load();
}

// Pass-2 pool: re-LCPs the same shards, fills byte slots for surviving sigs,
// AND emits within-pass adjacent-sig edges (only between surviving sigs). All
// edges land in `all_edges`; the pool stage no longer carries 40+ GB of raw
// edges through the pool peak.
void fill_pool(const std::vector<Shard>& shards, int level, int flank, int margin, int threads,
               const std::unordered_map<uint64_t, uint32_t>& idx,
               std::vector<std::vector<uint8_t>>& slot_buf,
               std::vector<std::atomic<uint8_t>>& claimed,
               std::vector<SEdge>& all_edges, std::mutex& emtx) {
    const int nshards = (int)shards.size();
    const int nrd = std::min(nshards, std::max(2, threads / 4));
    const int nw  = std::max(1, threads - nrd);
    std::atomic<size_t> next{0}, done{0};
    std::atomic<int> readers_left{nrd};
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<Batch> q;
    constexpr size_t Q_CAP = 256, BATCH = 64;
    constexpr unsigned BLOCK = 1u << 20;
    auto reader = [&]() {
        std::vector<char> buf(BLOCK);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string carry; int col_local = 0;
        auto flush = [&]() {
            if (batch.empty()) return;
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(Batch{std::move(batch), col_local}); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        size_t si;
        while ((si = next.fetch_add(1)) < shards.size()) {
            col_local = shards[si].col;
            gzFile f = gzopen(shards[si].path.c_str(), "rb");
            if (!f) { std::fprintf(stderr, "[FATAL] gzopen %s\n", shards[si].path.c_str()); std::abort(); }
            gzbuffer(f, 1u << 20);
            int line_no = 0; carry.clear();
            int got;
            while ((got = gzread(f, buf.data(), BLOCK)) > 0) {
                const char* base = buf.data(); const char* end = base + got; const char* start = base;
                for (const char* p = base; p < end; p++) if (*p == '\n') {
                    const size_t n = (size_t)(p - start);
                    if (line_no == 1) {
                        if (carry.empty()) batch.emplace_back(start, n);
                        else { carry.append(start, n); batch.push_back(std::move(carry)); carry.clear(); }
                        if (batch.size() >= BATCH) flush();
                    } else if (!carry.empty()) carry.clear();
                    line_no = (line_no + 1) & 3; start = p + 1;
                }
                if (start < end && line_no == 1) carry.append(start, end - start);
            }
            gzclose(f);
            std::fprintf(stderr, "[fill] shard %zu/%d done (col=%d)\n", (size_t)done.fetch_add(1) + 1, nshards, col_local);
        }
        flush();
        if (readers_left.fetch_sub(1) == 1) { std::lock_guard<std::mutex> lk(qmtx); cv_empty.notify_all(); }
    };
    auto worker = [&]() {
        std::string rc_buf; rc_buf.reserve(1u << 16);
        std::vector<SEdge> ledges; ledges.reserve(1u << 22);
        std::unordered_set<uint64_t> seen; seen.reserve(1u << 24);   // 16 M-entry edge-key cache (~640 MB / worker peak)
        constexpr size_t E_FLUSH = 1u << 22;     // 4M edges (~64 MB) flush trigger
        constexpr size_t SEEN_CAP = 1u << 24;    // cap seen at 16M; clear on overflow
        while (true) {
            Batch b;
            {
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || readers_left.load() == 0; });
                if (q.empty()) break;
                b = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            for (auto& s : b.seqs) process_read_fill(s, level, flank, margin, idx, slot_buf, claimed, ledges, seen, rc_buf);
            if (ledges.size() >= E_FLUSH) {
                std::lock_guard<std::mutex> lk(emtx);
                all_edges.insert(all_edges.end(), ledges.begin(), ledges.end());
                ledges.clear();
            }
            if (seen.size() >= SEEN_CAP) seen.clear();   // overflow-safe; loses some cross-batch dedup but bounds memory
        }
        { std::lock_guard<std::mutex> lk(emtx); all_edges.insert(all_edges.end(), ledges.begin(), ledges.end()); }
    };
    std::vector<std::thread> rds, wks;
    for (int i = 0; i < nrd; i++) rds.emplace_back(reader);
    for (int i = 0; i < nw;  i++) wks.emplace_back(worker);
    for (auto& t : rds) t.join();
    for (auto& t : wks) t.join();
}

// LSB radix sort over bare u64 sigs: 11-bit digit × 6 passes (66 ≥ 64). size_t hist +
// offsets — N exceeds 2^32 at WG scale. Half the element bytes + 2 fewer passes than the
// old 16-byte/8-pass SigCol sort. Output overwrites src.
void radix_sort_u64(std::vector<uint64_t>& sigs) {
    const size_t N = sigs.size();
    if (N == 0) return;
    std::vector<uint64_t> tmp(N);
    uint64_t *src = sigs.data(), *dst = tmp.data();
    constexpr int RADIX = 11, BUCKETS = 1 << RADIX;
    constexpr uint64_t MASK = BUCKETS - 1;
    size_t hist[BUCKETS];
    for (int pass = 0; pass < 6; pass++) {
        std::memset(hist, 0, sizeof(hist));
        const int sh = pass * RADIX;
        for (size_t i = 0; i < N; i++) hist[(src[i] >> sh) & MASK]++;
        size_t acc = 0;
        for (int b = 0; b < BUCKETS; b++) { size_t c = hist[b]; hist[b] = acc; acc += c; }
        for (size_t i = 0; i < N; i++) {
            const uint64_t bk = (src[i] >> sh) & MASK;
            dst[hist[bk]++] = src[i];
        }
        uint64_t *t = src; src = dst; dst = t;
    }
    if (src != sigs.data()) std::memcpy(sigs.data(), src, N * sizeof(uint64_t));
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

// Reduce sorted u64 sigs into per-signature SigCount via run-length count. Each pool pass
// is single-colour, so the run length is the n-count (normal) or t-count (tumor); the
// other colour stays 0 and is combined later by merge_sorted_nodes.
void reduce_sigs_u64(const std::vector<uint64_t>& sigs, std::vector<SigCount>& nodes, bool tumor) {
    nodes.reserve(sigs.size() / 2 + 1);
    for (size_t i = 0; i < sigs.size();) {
        const uint64_t g = sigs[i]; uint32_t c = 0;
        for (; i < sigs.size() && sigs[i] == g; i++) c++;
        nodes.push_back(tumor ? SigCount{g, 0, c} : SigCount{g, c, 0});
    }
}

// Permissive Cortex bubble = a maximal connected run of non-shared (colour != 3) nodes,
// bounded by shared-backbone nodes that are not traversed. One BFS per seed gives one
// bubble. type = 'D' (only colour-1 seen, normal-only/germline), 'I' (only colour-2,
// tumor-only/somatic candidate), 'C' (both → mixed). This is the 2026-05-26 high-recall
// decomposition that hit ~99.6% recall on WG reads; the chain+union-find-by-anchor merge
// that displaced it (2026-06-01) over-merged adjacent loci and dropped recall ~3-6 pts.
void decompose(const std::vector<std::vector<int32_t>>& adj, const std::vector<uint8_t>& color,
               const std::vector<uint32_t>& ncov, const std::vector<uint32_t>& tcov,
               std::vector<CortexBubble>& out) {
    const uint32_t K = (uint32_t)color.size();
    std::vector<char> vis(K, 0); out.clear();
    std::queue<uint32_t> q;
    for (uint32_t s = 0; s < K; s++) {
        if (color[s] == 3 || vis[s]) continue;
        q.push(s); vis[s] = 1;
        uint8_t seen = 0; uint32_t nn = 0; uint64_t nc = 0, tc = 0;
        while (!q.empty()) {
            const uint32_t u = q.front(); q.pop(); seen |= color[u]; nn++; nc += ncov[u]; tc += tcov[u];
            for (int32_t v : adj[u]) if (!vis[v] && color[v] != 3) { vis[v] = 1; q.push((uint32_t)v); }
        }
        out.push_back({(char)(seen == 1 ? 'D' : seen == 2 ? 'I' : 'C'), nn, nc, tc});
    }
}
// Prune (sig, n, t) records by min_mult on max(n,t). In-place: nodes is replaced
// with the kept prefix.
void prune_by_mult(std::vector<SigCount>& nodes, int min_mult) {
    size_t k = 0;
    for (const SigCount& c : nodes) if ((c.n > c.t ? c.n : c.t) >= (uint32_t)min_mult) nodes[k++] = c;
    nodes.resize(k);
}

// Concatenate slot_buf[i] into payload + record seqlen[i]. Frees each slot as it goes.
// Aborts on empty slot — pass 2 reads the same FASTQs as pass 1 so every kept sig must
// be re-seen at least once.
void pack_payload(std::vector<std::vector<uint8_t>>& slot_buf,
                  std::vector<uint16_t>& seqlen, std::vector<uint8_t>& payload) {
    const size_t K = slot_buf.size();
    seqlen.resize(K);
    uint64_t total = 0;
    for (const auto& v : slot_buf) total += v.size();
    payload.reserve(total);
    for (size_t i = 0; i < K; i++) {
        if (slot_buf[i].empty()) { std::fprintf(stderr, "[FATAL] empty slot %zu (pass-2 missed a kept sig)\n", i); std::abort(); }
        seqlen[i] = (uint16_t)slot_buf[i].size();
        payload.insert(payload.end(), slot_buf[i].begin(), slot_buf[i].end());
        std::vector<uint8_t>().swap(slot_buf[i]);
    }
}
}  // namespace

// Two-pointer merge of two sig-sorted SigCount vectors. Combines counts at the
// same sig — `a` has tumor counts in .t implicitly zero (normal pass), `b` has
// normal counts in .n implicitly zero (tumor pass). Caller frees a/b.
static std::vector<SigCount> merge_sorted_nodes(const std::vector<SigCount>& a, const std::vector<SigCount>& b) {
    std::vector<SigCount> out; out.reserve(a.size() + b.size());
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i].sig < b[j].sig)      { out.push_back(a[i++]); }
        else if (a[i].sig > b[j].sig) { out.push_back(b[j++]); }
        else { out.push_back({a[i].sig, a[i].n, b[j].t}); i++; j++; }
    }
    while (i < a.size()) out.push_back(a[i++]);
    while (j < b.size()) out.push_back(b[j++]);
    return out;
}

// ---- stage 1: LCP + count + prune + byte fill ----
// Per-sample two-pass: pool normal alone, sort+reduce+dedup-edges, free raw
// sig vector (bare u64), THEN pool tumor alone. Peak in-pool RAM = one sample's raw
// emissions (~50 % of the dual-extract total) instead of both samples combined.
// pre_prune.bin schema: u64 N + N × SigCount (sig, n, t). Holds the merged-but-
// unpruned sig table — lets a subsequent --min-mult sweep skip the ~22-min pool
// stage and jump straight to prune + fill_pool. If --level / --flank changed, the
// caller must delete pre_prune.bin first (no header guard yet).
static bool file_exists(const std::string& p) { FILE* f = std::fopen(p.c_str(), "rb"); if (f) { std::fclose(f); return true; } return false; }

void run_lcp(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int flank, int margin, int min_mult, int threads) {
    init_comp_lut();
    mkdir(dir.c_str(), 0755);
    const std::vector<std::string> ns = collect_shards(normal_fq), ts = collect_shards(tumor_fq);
    std::vector<Shard> shards_n, shards_t, shards_all;
    for (const auto& p : ns) shards_n.push_back({p, 0});
    for (const auto& p : ts) shards_t.push_back({p, 1});
    shards_all = shards_n; shards_all.insert(shards_all.end(), shards_t.begin(), shards_t.end());
    std::fprintf(stderr, "[lcp] normal=%zu shards, tumor=%zu shards, threads=%d min_mult=%d (per-sample two-pass)\n",
                 ns.size(), ts.size(), threads, min_mult);
    Tic tic("lcp");

    std::vector<SigCount> nodes;
    uint64_t nr_n = 0, nr_t = 0;
    const std::string pp_path = dir + "/pre_prune.bin";

    if (file_exists(pp_path)) {
        // ---- RESUME: load pre_prune checkpoint, skip pool/merge ----
        FILE* fp = xopen(pp_path, "rb");
        uint64_t N; rd(fp, &N, 8);
        nodes.resize(N);
        rd(fp, nodes.data(), N * sizeof(SigCount));
        std::fclose(fp);
        tic.mark("load_pre_prune");
        std::fprintf(stderr, "[lcp]   loaded pre_prune.bin: %lu nodes (skipped pool/reduce/merge)\n", (unsigned long)N);
        // Reconstruct nr_n / nr_t from sums (best effort; not strictly accurate but informative).
        for (const SigCount& c : nodes) { nr_n += (c.n > 0); nr_t += (c.t > 0); }
    } else {
        std::vector<SigCount> nodes_n, nodes_t;
        // ---- PASS A: normal (sigs only, no edges) ----
        {
            std::vector<uint64_t> sigs; std::mutex smtx;
            uint64_t dummy_t = 0;
            run_pool(shards_n, level, flank, margin, threads, sigs, smtx, nr_n, dummy_t);
            tic.mark("pool_n");
            std::fprintf(stderr, "[lcp/n] sigs raw=%lu (%.1f GB)\n",
                         (unsigned long)sigs.size(), sigs.size() * sizeof(uint64_t) / 1e9);
            radix_sort_u64(sigs);
            tic.mark("radix_n");
            reduce_sigs_u64(sigs, nodes_n, false);
            std::vector<uint64_t>().swap(sigs);
            tic.mark("reduce_n");
            std::fprintf(stderr, "[lcp/n] nodes=%lu (%.1f GB)\n",
                         (unsigned long)nodes_n.size(), nodes_n.size() * sizeof(SigCount) / 1e9);
        }

        // ---- PASS B: tumor (sigs only) ----
        {
            std::vector<uint64_t> sigs; std::mutex smtx;
            uint64_t dummy_n = 0;
            run_pool(shards_t, level, flank, margin, threads, sigs, smtx, dummy_n, nr_t);
            tic.mark("pool_t");
            std::fprintf(stderr, "[lcp/t] sigs raw=%lu (%.1f GB)\n",
                         (unsigned long)sigs.size(), sigs.size() * sizeof(uint64_t) / 1e9);
            radix_sort_u64(sigs);
            tic.mark("radix_t");
            reduce_sigs_u64(sigs, nodes_t, true);
            std::vector<uint64_t>().swap(sigs);
            tic.mark("reduce_t");
            std::fprintf(stderr, "[lcp/t] nodes=%lu (%.1f GB)\n",
                         (unsigned long)nodes_t.size(), nodes_t.size() * sizeof(SigCount) / 1e9);
        }

        nodes = merge_sorted_nodes(nodes_n, nodes_t);
        std::vector<SigCount>().swap(nodes_n);
        std::vector<SigCount>().swap(nodes_t);
        tic.mark("merge_nodes");
        std::fprintf(stderr, "[lcp]   merged nodes=%lu (%.1f GB)\n",
                     (unsigned long)nodes.size(), nodes.size() * sizeof(SigCount) / 1e9);

        // ---- CHECKPOINT: save pre_prune.bin so future min_mult sweeps skip everything above ----
        FILE* fp = xopen(pp_path + ".tmp", "wb");
        uint64_t Npp = nodes.size(); wr(fp, &Npp, 8);
        wr(fp, nodes.data(), Npp * sizeof(SigCount));
        atomic_close(fp, pp_path + ".tmp", pp_path);
        tic.mark("save_pre_prune");
        std::fprintf(stderr, "[lcp]   wrote pre_prune.bin: %lu nodes (%.1f GB) -> %s\n",
                     (unsigned long)Npp, Npp * sizeof(SigCount) / 1e9, pp_path.c_str());
    }

    const uint64_t pre_prune = nodes.size();
    prune_by_mult(nodes, min_mult);
    tic.mark("prune");
    std::fprintf(stderr, "[lcp]   prune min_mult=%d kept=%lu of %lu (%.2f%%)\n",
                 min_mult, (unsigned long)nodes.size(), (unsigned long)pre_prune, 100.0 * nodes.size() / pre_prune);
    std::vector<SigCount>(nodes).swap(nodes);   // prune drops ~95% of nodes but resize keeps the capacity;
                                                // release the dead allocation (tens of GB at WG scale) so it
                                                // isn't carried through fill_pool's idx/slot_buf/edges peak.
    std::unordered_map<uint64_t, uint32_t> idx; idx.reserve(nodes.size());
    for (uint32_t i = 0; i < nodes.size(); i++) idx.emplace(nodes[i].sig, i);
    std::vector<std::vector<uint8_t>> slot_buf(nodes.size());
    std::vector<std::atomic<uint8_t>> claimed(nodes.size());
    for (auto& c : claimed) c.store(0);

    std::vector<SEdge> edges; std::mutex emtx;
    fill_pool(shards_all, level, flank, margin, threads, idx, slot_buf, claimed, edges, emtx);
    tic.mark("fill_pool");
    std::fprintf(stderr, "[lcp]   fill: edges raw=%lu (%.1f GB)\n",
                 (unsigned long)edges.size(), edges.size() * sizeof(SEdge) / 1e9);
    std::vector<uint16_t> seqlen; std::vector<uint8_t> payload;
    pack_payload(slot_buf, seqlen, payload);
    tic.mark("pack_payload");
    radix_sort_edges(edges);
    tic.mark("radix_edges");
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    tic.mark("dedup_edges");

    const std::string cp = dir + "/counts.bin", ep = dir + "/edges.bin";
    FILE* fc = xopen(cp + ".tmp", "wb");
    uint64_t N = nodes.size(); wr(fc, &N, 8); wr(fc, nodes.data(), N * sizeof(SigCount));
    uint64_t total = payload.size(); wr(fc, &total, 8); wr(fc, seqlen.data(), N * 2); wr(fc, payload.data(), total);
    atomic_close(fc, cp + ".tmp", cp);
    FILE* fe = xopen(ep + ".tmp", "wb"); uint64_t M = edges.size(); wr(fe, &M, 8); wr(fe, edges.data(), M * sizeof(SEdge)); atomic_close(fe, ep + ".tmp", ep);
    tic.mark("write");
    std::fprintf(stderr, "[lcp] normal_reads=%lu tumor_reads=%lu kept_sigs=%lu (payload=%.1f GB) edges=%lu -> %s/{counts,edges}.bin\n",
                 (unsigned long)nr_n, (unsigned long)nr_t, (unsigned long)N, total / 1e9, (unsigned long)M, dir.c_str());
}

// ---- stage 2: colour + resolve edges (counts.bin is already pruned in lcp stage) ----
void run_graph(const std::string& dir) {
    Tic tic("graph");
    FILE* fc = xopen(dir + "/counts.bin", "rb"); uint64_t N; rd(fc, &N, 8); std::vector<SigCount> nodes(N); rd(fc, nodes.data(), N * sizeof(SigCount)); std::fclose(fc);
    tic.mark("read_counts");
    std::unordered_map<uint64_t, uint32_t> id; id.reserve(N);
    std::vector<uint8_t> color(N); std::vector<uint32_t> ncov(N), tcov(N);
    for (uint64_t i = 0; i < N; i++) {
        const SigCount& c = nodes[i];
        id.emplace(c.sig, (uint32_t)i);
        color[i] = c.n && c.t ? 3 : c.n ? 1 : 2; ncov[i] = c.n; tcov[i] = c.t;
    }
    tic.mark("colour");
    const uint64_t kept = N;
    std::vector<SigCount>().swap(nodes);
    FILE* fe = xopen(dir + "/edges.bin", "rb"); uint64_t M; rd(fe, &M, 8); std::vector<SEdge> se(M); rd(fe, se.data(), M * sizeof(SEdge)); std::fclose(fe);
    tic.mark("read_edges");
    std::vector<IEdge> ie; ie.reserve(M);
    for (const SEdge& e : se) {
        auto a = id.find(e.first), b = id.find(e.second);
        if (a != id.end() && b != id.end() && a->second != b->second)
            ie.push_back(IEdge(a->second, b->second));
    }
    tic.mark("resolve_edges");
    std::sort(ie.begin(), ie.end()); ie.erase(std::unique(ie.begin(), ie.end()), ie.end());
    tic.mark("sort+dedup");
    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen(gp + ".tmp", "wb"); uint64_t K = kept; wr(fg, &K, 8);
    wr(fg, color.data(), K); wr(fg, ncov.data(), K * 4); wr(fg, tcov.data(), K * 4);
    uint64_t E = ie.size(); wr(fg, &E, 8); wr(fg, ie.data(), E * sizeof(IEdge)); atomic_close(fg, gp + ".tmp", gp);
    tic.mark("write");
    std::fprintf(stderr, "[graph] nodes=%lu edges=%lu -> %s/graph.bin\n",
                 (unsigned long)kept, (unsigned long)E, dir.c_str());
}

// ---- stage 3: decompose into bubbles ----
void run_bubbles(const std::string& dir, const std::string& out_tsv) {
    Tic tic("bubbles");
    FILE* fg = xopen(dir + "/graph.bin", "rb"); uint64_t K; rd(fg, &K, 8);
    std::vector<uint8_t> color(K); rd(fg, color.data(), K);
    std::vector<uint32_t> ncov(K), tcov(K); rd(fg, ncov.data(), K * 4); rd(fg, tcov.data(), K * 4);
    uint64_t E; rd(fg, &E, 8); std::vector<IEdge> ie(E); rd(fg, ie.data(), E * sizeof(IEdge)); std::fclose(fg);
    tic.mark("read_graph");
    std::vector<std::vector<int32_t>> adj(K);
    for (const IEdge& e : ie) { adj[e.first].push_back((int32_t)e.second); adj[e.second].push_back((int32_t)e.first); }
    tic.mark("build_adj");
    std::vector<CortexBubble> out; decompose(adj, color, ncov, tcov, out);
    tic.mark("bfs_decomp");
    FILE* f = xopen(out_tsv, "w"); std::fputs("type\tn_nodes\tnormal_cov\ttumor_cov\n", f);
    int nd = 0, ni = 0, nc = 0;
    for (const CortexBubble& b : out) {
        const char* t = b.type == 'D' ? "DEL" : b.type == 'I' ? "INS" : "CPLX";
        nd += b.type == 'D'; ni += b.type == 'I'; nc += b.type == 'C';
        std::fprintf(f, "%s\t%u\t%lu\t%lu\n", t, b.n_nodes, (unsigned long)b.normal_cov, (unsigned long)b.tumor_cov);
    }
    std::fclose(f);
    tic.mark("write_tsv");
    std::fprintf(stderr, "[bubbles] %zu bubbles  DEL(normal-only)=%d INS(somatic)=%d CPLX=%d -> %s\n", out.size(), nd, ni, nc, out_tsv.c_str());
}
