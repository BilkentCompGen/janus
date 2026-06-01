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
#include <utility>
#include <vector>

namespace {
struct SigCount { uint64_t sig; uint32_t n, t; };   // counts.bin record (16 B, packed)
using SEdge = std::pair<uint64_t, uint64_t>;        // signature edge (edges.bin record)
using IEdge = std::pair<uint32_t, uint32_t>;        // node-id edge (graph.bin record)
struct SigCol { uint64_t sig; uint32_t col; };      // per-emission record collected by workers
                                                     // (col 0=normal, 1=tumor); reduced by sort.

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

// ---- fail-fast binary checkpoint IO ----
FILE* xopen(const std::string& p, const char* m) { FILE* f = std::fopen(p.c_str(), m); if (!f) { std::fprintf(stderr, "[FATAL] fopen %s\n", p.c_str()); std::abort(); } return f; }
void wr(FILE* f, const void* p, size_t n) { if (n && std::fwrite(p, 1, n, f) != n) { std::fprintf(stderr, "[FATAL] fwrite\n"); std::abort(); } }
void rd(FILE* f, void* p, size_t n) { if (n && std::fread(p, 1, n, f) != n) { std::fprintf(stderr, "[FATAL] truncated checkpoint\n"); std::abort(); } }
void atomic_close(FILE* f, const std::string& tmp, const std::string& dst) { std::fclose(f); if (std::rename(tmp.c_str(), dst.c_str())) { std::fprintf(stderr, "[FATAL] rename %s\n", dst.c_str()); std::abort(); } }   // never leave a half-written checkpoint resume would trust

// One read -> interior signatures + within-read adjacent-signature edges.
// Bypasses extract_cores/compute_sigs intermediate vectors; reads the SoA
// arrays exposed by extract_fast and inlines the flank-window murmur64.
void process_read(std::string& s, int col, int level, int flank, int margin,
                  std::vector<SigCol>& lsigs, std::vector<SEdge>& ledges) {
    for (char& ch : s) if (ch >= 'a' && ch <= 'z') ch -= 32;
    const uint64_t len = s.size();
    if (len < (uint64_t)(2 * margin + 4)) return;
    int n; const uint64_t *st, *en; const uint32_t *lb;
    extract_fast(s.data(), (int)len, 0, level, &n, &st, &en, &lb);
    const uint64_t mu = (uint64_t)flank, mg = (uint64_t)margin;
    uint64_t prev = 0; bool have = false;
    for (int i = 0; i < n; i++) {
        if (st[i] < mg || en[i] + mg > len) continue;
        const uint64_t L = st[i] > mu ? st[i] - mu : 0;
        const uint64_t R = en[i] + mu < len ? en[i] + mu : len;
        const uint64_t g = murmur64(s.data() + L, (int)(R - L), lb[i]);
        lsigs.push_back({g, (uint32_t)col});
        if (have && prev != g) ledges.push_back(prev < g ? SEdge(prev, g) : SEdge(g, prev));
        prev = g; have = true;
    }
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
              std::vector<SigCol>& all_sigs, std::mutex& smtx,
              std::vector<SEdge>& all_edges, std::mutex& emtx,
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
        std::vector<SEdge> ledges; std::vector<SigCol> lsigs;
        lsigs.reserve(1u << 22); ledges.reserve(1u << 22);
        uint64_t lnr_n = 0, lnr_t = 0;
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
            for (auto& s : b.seqs) { process_read(s, b.col, level, flank, margin, lsigs, ledges); (b.col == 0 ? lnr_n : lnr_t)++; }
            local_busy_ns += (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_b0).count();
        }
        { std::lock_guard<std::mutex> lk(emtx); all_edges.insert(all_edges.end(), ledges.begin(), ledges.end()); }
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

// LSB radix sort over the .sig field (8 passes × 256 buckets). Output overwrites src.
void radix_sort_sigs(std::vector<SigCol>& sigs) {
    const size_t N = sigs.size();
    if (N == 0) return;
    std::vector<SigCol> tmp(N);
    SigCol *src = sigs.data(), *dst = tmp.data();
    uint32_t hist[256];
    for (int pass = 0; pass < 8; pass++) {
        std::memset(hist, 0, sizeof(hist));
        const int sh = pass * 8;
        for (size_t i = 0; i < N; i++) hist[(src[i].sig >> sh) & 0xFFu]++;
        uint32_t acc = 0;
        for (int b = 0; b < 256; b++) { uint32_t c = hist[b]; hist[b] = acc; acc += c; }
        for (size_t i = 0; i < N; i++) {
            const uint32_t bk = (uint32_t)((src[i].sig >> sh) & 0xFFu);
            dst[hist[bk]++] = src[i];
        }
        SigCol *t = src; src = dst; dst = t;
    }
    if (src != sigs.data()) std::memcpy(sigs.data(), src, N * sizeof(SigCol));
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

// Reduce sorted SigCols into per-signature {n,t} count records (run-length merge).
void reduce_sigs(const std::vector<SigCol>& sigs, std::vector<SigCount>& nodes) {
    nodes.reserve(sigs.size() / 2 + 1);
    for (size_t i = 0; i < sigs.size();) {
        uint64_t g = sigs[i].sig; uint32_t n = 0, t = 0;
        for (; i < sigs.size() && sigs[i].sig == g; i++) (sigs[i].col == 0 ? n : t)++;
        nodes.push_back({g, n, t});
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
}  // namespace

// ---- stage 1: LCP + count (the expensive, gz-bound pass) ----
void run_lcp(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int flank, int margin, int threads) {
    mkdir(dir.c_str(), 0755);
    const std::vector<std::string> ns = collect_shards(normal_fq), ts = collect_shards(tumor_fq);
    std::vector<Shard> shards;
    for (const auto& p : ns) shards.push_back({p, 0});
    for (const auto& p : ts) shards.push_back({p, 1});
    std::fprintf(stderr, "[lcp] normal=%zu shards, tumor=%zu shards, threads=%d\n", ns.size(), ts.size(), threads);
    Tic tic("lcp");
    std::vector<SigCol> sigs; std::mutex smtx;
    std::vector<SEdge> edges; std::mutex emtx;
    uint64_t nr_n = 0, nr_t = 0;
    run_pool(shards, level, flank, margin, threads, sigs, smtx, edges, emtx, nr_n, nr_t);
    tic.mark("pool");
    std::fprintf(stderr, "[lcp]    sigs raw=%lu (%.1f GB)  edges raw=%lu (%.1f GB)\n",
                 (unsigned long)sigs.size(), sigs.size() * sizeof(SigCol) / 1e9,
                 (unsigned long)edges.size(), edges.size() * sizeof(SEdge) / 1e9);
    radix_sort_sigs(sigs);
    tic.mark("radix_sigs");
    std::vector<SigCount> nodes; reduce_sigs(sigs, nodes);
    tic.mark("reduce_sigs");
    std::vector<SigCol>().swap(sigs);
    radix_sort_edges(edges);
    tic.mark("radix_edges");
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    tic.mark("dedup_edges");
    const std::string cp = dir + "/counts.bin", ep = dir + "/edges.bin";
    FILE* fc = xopen(cp + ".tmp", "wb"); uint64_t N = nodes.size(); wr(fc, &N, 8); wr(fc, nodes.data(), N * sizeof(SigCount)); atomic_close(fc, cp + ".tmp", cp);
    FILE* fe = xopen(ep + ".tmp", "wb"); uint64_t M = edges.size(); wr(fe, &M, 8); wr(fe, edges.data(), M * sizeof(SEdge)); atomic_close(fe, ep + ".tmp", ep);
    tic.mark("write");
    std::fprintf(stderr, "[lcp] normal_reads=%lu tumor_reads=%lu distinct_sigs=%lu edges=%lu -> %s/{counts,edges}.bin\n",
                 (unsigned long)nr_n, (unsigned long)nr_t, (unsigned long)N, (unsigned long)M, dir.c_str());
}

// ---- stage 2: prune + colour + resolve edges ----
void run_graph(const std::string& dir, int min_mult) {
    Tic tic("graph");
    FILE* fc = xopen(dir + "/counts.bin", "rb"); uint64_t N; rd(fc, &N, 8); std::vector<SigCount> nodes(N); rd(fc, nodes.data(), N * sizeof(SigCount)); std::fclose(fc);
    tic.mark("read_counts");
    std::unordered_map<uint64_t, uint32_t> id; id.reserve(N);
    std::vector<uint8_t> color; std::vector<uint32_t> ncov, tcov;
    for (const SigCount& c : nodes) {
        if ((c.n > c.t ? c.n : c.t) < (uint32_t)min_mult) continue;
        id.emplace(c.sig, (uint32_t)color.size());
        color.push_back(c.n && c.t ? 3 : c.n ? 1 : 2); ncov.push_back(c.n); tcov.push_back(c.t);
    }
    tic.mark("prune+colour");
    const uint64_t kept = color.size();
    std::vector<SigCount>().swap(nodes);
    FILE* fe = xopen(dir + "/edges.bin", "rb"); uint64_t M; rd(fe, &M, 8); std::vector<SEdge> se(M); rd(fe, se.data(), M * sizeof(SEdge)); std::fclose(fe);
    tic.mark("read_edges");
    std::vector<IEdge> ie; ie.reserve(M);
    for (const SEdge& e : se) {
        auto a = id.find(e.first), b = id.find(e.second);
        if (a != id.end() && b != id.end() && a->second != b->second)
            ie.push_back(a->second < b->second ? IEdge(a->second, b->second) : IEdge(b->second, a->second));
    }
    tic.mark("resolve_edges");
    std::sort(ie.begin(), ie.end()); ie.erase(std::unique(ie.begin(), ie.end()), ie.end());
    tic.mark("sort+dedup");
    const std::string gp = dir + "/graph.bin";
    FILE* fg = xopen(gp + ".tmp", "wb"); uint64_t K = kept; wr(fg, &K, 8);
    wr(fg, color.data(), K); wr(fg, ncov.data(), K * 4); wr(fg, tcov.data(), K * 4);
    uint64_t E = ie.size(); wr(fg, &E, 8); wr(fg, ie.data(), E * sizeof(IEdge)); atomic_close(fg, gp + ".tmp", gp);
    tic.mark("write");
    std::fprintf(stderr, "[graph] min_mult=%d kept_nodes=%lu of %lu edges=%lu -> %s/graph.bin\n",
                 min_mult, (unsigned long)kept, (unsigned long)N, (unsigned long)E, dir.c_str());
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
