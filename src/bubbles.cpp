// snape bubbles <graph.bin> <normal> <tumour> [threads=16] [k=21] [w=22] [maxlen=40] [max_dlen=12] [max_unmatched=20]
//
// Most of the somatic material in the graph is not a junction. It is SNVs and sub-SV edits, and each
// one still becomes a chunk, a contig, and a minute of alignment. This stage names those cores before
// any contig exists, so the aligner only ever sees candidates that could be an SV.
//
// A BUBBLE is a source S with out-degree >= 2, two node-disjoint linear arms, and a sink T where they
// converge -- one arm somatic (ncov == 0), the other germline. Anything else is not a bubble and is
// left alone.
//
//         S --[ somatic arm: c1 c2 ... ]--> T
//          \--[ germline arm: g1 g2 ... ]--^
//
// An L5 core spans ~217 bp, so a core signature cannot see a 50 bp change: two cores match or they do
// not. Minimizers (k=21 w=22, ~13 per core) open the arms up, and the comparison is ORDERED -- a set
// intersection measures nothing here, because a junction core is built from S's tail and T's head and
// both of those are present in the germline arm anyway. What separates the two cases is colinearity:
//
//   SNV / small edit : same sequence with a local change. The minimizer sequences agree except for a
//                      short run in the middle -> the arms are the same length and the longest common
//                      subsequence covers nearly all of the somatic arm.
//   junction         : different sequence. Minimizers may still be shared, but they do not stay in
//                      order, and the LCS collapses.
//
// Hence the two thresholds: `max_dlen` bounds the arm length difference in minimizers and
// `max_unmatched` how many somatic minimizers may fall outside the LCS. Measured on `hard`: 12/20
// drops 83% of the somatic cores and 84% of the contigs for 2 junctions of recall.
#include "bubbles.h"

#include "snape.h"
#include "unordered_dense.h"

extern "C" {
#include "encoding.h"
}

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

namespace {

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[bubbles] FATAL: %s\n", m); std::exit(1); }

int LEVEL = 5;                                 // runtime, from SNAPE_LEVEL
constexpr int MARGIN = 10;                     // MUST match the graph build

inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return x;
}

inline int b2(char c) {
    switch (c) { case 'A': case 'a': return 0; case 'C': case 'c': return 1;
                 case 'G': case 'g': return 2; case 'T': case 't': return 3; }
    return -1;
}

// Canonical (strand-independent) minimizers, one per w-window, duplicates of the previous pick
// dropped so the result is the window's sequence of choices rather than one entry per position.
void minimizers(const char* s, int len, int k, int w, std::vector<uint64_t>& out) {
    out.clear();
    if (len < k) return;
    const uint64_t mask = (k >= 32) ? ~0ull : ((1ull << (2 * k)) - 1);
    uint64_t fw = 0, rv = 0; int have = 0;
    static thread_local std::vector<uint64_t> h;
    h.clear(); h.reserve((size_t)len);
    for (int i = 0; i < len; i++) {
        const int c = b2(s[i]);
        if (c < 0) { have = 0; fw = rv = 0; h.push_back(~0ull); continue; }
        fw = ((fw << 2) | (uint64_t)c) & mask;
        rv = (rv >> 2) | ((uint64_t)(3 - c) << (2 * (k - 1)));
        if (++have >= k) h.push_back(mix64(fw < rv ? fw : rv)); else h.push_back(~0ull);
    }
    uint64_t last = ~0ull;
    for (size_t i = 0; i + (size_t)w <= h.size(); i++) {
        uint64_t best = ~0ull;
        for (int j = 0; j < w; j++) best = std::min(best, h[i + j]);
        if (best != ~0ull && best != last) { out.push_back(best); last = best; }
    }
}

// Homopolymer-compressed, so a HiFi run-length wobble does not read as a sequence difference.
void hpc_local(const char* s, int len, std::string& out) {
    out.clear();
    for (int i = 0; i < len; i++) if (i == 0 || s[i] != s[i - 1]) out.push_back(s[i]);
}

// The cores this stage needs sequence for, each given a dense id so a hit costs one lookup and no
// allocation. Membership is read-only once built, so every worker can probe it at once.
struct Want {
    ankerl::unordered_dense::map<uint64_t, uint32_t> id;
    uint32_t add(uint64_t sig) {
        const uint32_t n = (uint32_t)id.size();
        return id.emplace(sig, n).first->second;
    }
    uint32_t find(uint64_t sig) const {
        const auto it = id.find(sig);
        return it == id.end() ? UINT32_MAX : it->second;
    }
    size_t size() const { return id.size(); }
};

// One pass over one sample's reads: any read covering a wanted core gives that core its minimizers.
// WHICH read cannot matter -- equal signatures mean equal homopolymer-compressed bytes, and the
// minimizers are computed from those bytes.
//
// One reader thread fills batches of reads; the workers take a whole batch at a time and never touch
// the input again. A hit needs no lock either: `seen` is one atomic flag per wanted core and the
// worker that claims a flag writes that core's minimizers.
// [Before: every worker took the input mutex per read AND copied the read while holding it, then
//  took a global map mutex twice per hit. Measured on hcc1954: 21.9 M voluntary context switches
//  and 732% of the 1600% the run was given.]
void collect(const std::string& path, const Want& want, int NT,
             std::vector<std::vector<uint64_t>>& mins_of, int k, int w, const char* tag) {
    const size_t W = want.size();
    mins_of.assign(W, {});
    std::unique_ptr<std::atomic<uint8_t>[]> seen(new std::atomic<uint8_t>[W ? W : 1]);
    for (size_t i = 0; i < W; i++) seen[i].store(0, std::memory_order_relaxed);

    constexpr size_t Q_CAP = 128, BATCH = 32;
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<std::vector<std::string>> q;
    bool eof = false;
    std::atomic<uint64_t> nread{0}, ncov{0};

    auto reader = [&]() {
        Reads in = reads_open(path);
        std::vector<std::string> batch; batch.reserve(BATCH);
        std::string seq;
        auto push = [&]() {
            if (batch.empty()) return;
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(std::move(batch)); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        while (reads_next(in, seq, nullptr)) {
            batch.push_back(std::move(seq));
            if (batch.size() >= BATCH) push();
        }
        push();
        reads_close(in);
        { std::lock_guard<std::mutex> lk(qmtx); eof = true; }
        cv_empty.notify_all();
    };

    auto worker = [&]() {
        std::string hs; std::vector<uint64_t> mins;
        for (;;) {
            std::vector<std::string> batch;
            {
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || eof; });
                if (q.empty()) return;
                batch = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            for (const std::string& seq : batch) {
                nread++;
                int n; const uint64_t *st, *en, *sg;
                extract(seq.data(), (int)seq.size(), 0, LEVEL, &n, &st, &en, &sg);
                for (int i = 0; i < n; i++) {
                    if (st[i] < (uint64_t)MARGIN || en[i] + MARGIN > seq.size()) continue;
                    const uint32_t id = want.find(sg[i]);
                    if (id == UINT32_MAX || seen[id].load(std::memory_order_relaxed)) continue;
                    hpc_local(seq.data() + st[i], (int)(en[i] - st[i]), hs);
                    minimizers(hs.data(), (int)hs.size(), k, w, mins);
                    if (seen[id].exchange(1, std::memory_order_acq_rel) == 0) { mins_of[id] = mins; ncov++; }
                }
            }
        }
    };

    std::vector<std::thread> th;
    th.emplace_back(reader);
    for (int i = 1; i < std::max(2, NT); i++) th.emplace_back(worker);
    for (auto& t : th) t.join();
    std::fprintf(stderr, "[bubbles] %s: %llu okuma -> %zu/%zu cekirdek kapsandi\n",
                 tag, (unsigned long long)nread.load(), (size_t)ncov.load(), W);
}

size_t lcs(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    if (a.empty() || b.empty()) return 0;
    std::vector<uint16_t> prev(b.size() + 1, 0), cur(b.size() + 1, 0);
    for (size_t i = 1; i <= a.size(); i++) {
        for (size_t j = 1; j <= b.size(); j++)
            cur[j] = (a[i - 1] == b[j - 1]) ? (uint16_t)(prev[j - 1] + 1) : std::max(prev[j], cur[j - 1]);
        std::swap(prev, cur);
    }
    return prev[b.size()];
}

// ---- graph ----
struct Graph {
    uint64_t K = 0;
    std::vector<uint64_t> sig;
    std::vector<uint8_t> som;
    std::vector<uint64_t> oo, io;               // directed CSR: a bubble is a DIRECTED structure
    std::vector<uint32_t> on, in;
};

// graph.bin, streamed: the edge section is read twice in 1 M-edge chunks (degrees, then fill)
// instead of being materialised as one CEdge array. That array was 1.2 GB on hcc1954 and every byte
// of it was thrown away once the CSR existed. The per-colour counts go the same way -- they exist
// only to decide `som`, so they are consumed a chunk at a time and never all held at once.
Graph load_graph(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb"); if (!f) die("open graph.bin");
    Graph g;
    if (std::fread(&g.K, 8, 1, f) != 1) die("read K");
    g.sig.resize(g.K);
    if (std::fread(g.sig.data(), 8, g.K, f) != g.K) die("read sigs");

    g.som.assign(g.K, 0);
    constexpr size_t CH = 1u << 20;
    std::vector<uint32_t> buf(CH);
    for (uint64_t done = 0; done < g.K; ) {                       // ncov: som candidate iff 0
        const size_t m = (size_t)std::min<uint64_t>(CH, g.K - done);
        if (std::fread(buf.data(), 4, m, f) != m) die("read ncov");
        for (size_t i = 0; i < m; i++) g.som[done + i] = (buf[i] == 0);
        done += m;
    }
    for (uint64_t done = 0; done < g.K; ) {                       // tcov: and only if tumour saw it
        const size_t m = (size_t)std::min<uint64_t>(CH, g.K - done);
        if (std::fread(buf.data(), 4, m, f) != m) die("read tcov");
        for (size_t i = 0; i < m; i++) g.som[done + i] = (uint8_t)(g.som[done + i] & (buf[i] > 0));
        done += m;
    }

    uint64_t E = 0;
    if (std::fread(&E, 8, 1, f) != 1) die("read E");
    // The fill cursors below are 32-bit, which halves what they cost on a 75 M-edge graph. Node ids
    // are 32-bit anyway (CEdge), so a graph this loader could not take is one nothing else could
    // read either -- but say so rather than truncate.
    if (E > 0xffffffffull) die("graph.bin has more than 2^32 edges");
    const long epos = std::ftell(f);
    std::vector<uint32_t> od(g.K + 1, 0), idg(g.K + 1, 0);
    std::vector<CEdge> ce(CH);
    for (uint64_t done = 0; done < E; ) {                         // pass 1: degrees
        const size_t m = (size_t)std::min<uint64_t>(CH, E - done);
        if (std::fread(ce.data(), sizeof(CEdge), m, f) != m) die("read edges");
        for (size_t i = 0; i < m; i++) { od[ce[i].src]++; idg[ce[i].dst]++; }
        done += m;
    }
    g.oo.assign(g.K + 1, 0); g.io.assign(g.K + 1, 0);
    for (uint64_t i = 0; i < g.K; i++) { g.oo[i + 1] = g.oo[i] + od[i]; g.io[i + 1] = g.io[i] + idg[i]; }
    std::vector<uint32_t>().swap(od); std::vector<uint32_t>().swap(idg);
    g.on.resize(g.oo[g.K]); g.in.resize(g.io[g.K]);
    std::vector<uint32_t> c1(g.oo.begin(), g.oo.end()), c2(g.io.begin(), g.io.end());
    if (std::fseek(f, epos, SEEK_SET) != 0) die("seek edges");
    for (uint64_t done = 0; done < E; ) {                         // pass 2: fill
        const size_t m = (size_t)std::min<uint64_t>(CH, E - done);
        if (std::fread(ce.data(), sizeof(CEdge), m, f) != m) die("read edges");
        for (size_t i = 0; i < m; i++) { g.on[c1[ce[i].src]++] = ce[i].dst; g.in[c2[ce[i].dst]++] = ce[i].src; }
        done += m;
    }
    std::fclose(f);
    return g;
}

// Follow one linear arm until it converges on a node with in-degree >= 2 (the sink) or stops being
// linear. Returns the sink, or UINT32_MAX when this is not an arm.
uint32_t trace_arm(const Graph& g, uint32_t start, int maxlen, std::vector<uint32_t>& arm) {
    arm.clear();
    uint32_t cur = start;
    for (int s = 0; s <= maxlen; s++) {
        arm.push_back(cur);
        if (g.oo[cur + 1] - g.oo[cur] != 1) return UINT32_MAX;
        const uint32_t nxt = g.on[g.oo[cur]];
        const uint64_t din = g.io[nxt + 1] - g.io[nxt];
        if (din >= 2) return nxt;
        if (din != 1) return UINT32_MAX;
        cur = nxt;
    }
    return UINT32_MAX;
}

bool uniform(const Graph& g, const std::vector<uint32_t>& arm, uint8_t want) {
    for (uint32_t v : arm) if (g.som[v] != want) return false;
    return true;
}

struct Bub { std::vector<uint32_t> som_arm, germ_arm; };

// Every branching node is independent of every other one, so the sweep is split into contiguous
// ranges and the ranges are concatenated in range order -- the same bubbles in the same order as the
// single-threaded sweep produced, which is what keeps drop.sigs byte-for-byte the same.
std::vector<Bub> find_bubbles(const Graph& g, int maxlen, int NT) {
    const int nt = std::max(1, NT);
    std::vector<std::vector<Bub>> part(nt);
    std::vector<uint64_t> nsrc(nt, 0);
    std::vector<std::thread> th;
    for (int t = 0; t < nt; t++) th.emplace_back([&, t]{
        const uint64_t lo = g.K * t / nt, hi = g.K * (t + 1) / nt;
        std::vector<uint32_t> a1, a2;
        for (uint64_t s = lo; s < hi; s++) {
            if (g.oo[s + 1] - g.oo[s] < 2) continue;
            nsrc[t]++;
            for (uint64_t i = g.oo[s]; i < g.oo[s + 1]; i++) {
                const uint32_t T1 = trace_arm(g, g.on[i], maxlen, a1);
                if (T1 == UINT32_MAX) continue;
                for (uint64_t j = i + 1; j < g.oo[s + 1]; j++) {
                    if (trace_arm(g, g.on[j], maxlen, a2) != T1) continue;
                    if (g.som[a1[0]] == g.som[a2[0]]) continue;      // need exactly one somatic arm
                    if (!uniform(g, a1, g.som[a1[0]]) || !uniform(g, a2, g.som[a2[0]])) continue;
                    part[t].push_back(g.som[a1[0]] ? Bub{a1, a2} : Bub{a2, a1});
                }
            }
        }
    });
    for (auto& t : th) t.join();
    std::vector<Bub> out;
    uint64_t total = 0, srcs = 0;
    for (int t = 0; t < nt; t++) { total += part[t].size(); srcs += nsrc[t]; }
    out.reserve(total);
    for (int t = 0; t < nt; t++) {
        for (Bub& b : part[t]) out.push_back(std::move(b));
        std::vector<Bub>().swap(part[t]);
    }
    std::fprintf(stderr, "[bubbles] %llu dugum, %llu dallanma -> %zu temiz bubble (maxlen=%d)\n",
                 (unsigned long long)g.K, (unsigned long long)srcs, out.size(), maxlen);
    return out;
}

void arm_mins(const Graph& g, const std::vector<uint32_t>& arm, const Want& want,
              const std::vector<std::vector<uint64_t>>& mins, std::vector<uint64_t>& out) {
    out.clear();
    for (uint32_t v : arm) {
        const uint32_t id = want.find(g.sig[v]);
        if (id != UINT32_MAX) out.insert(out.end(), mins[id].begin(), mins[id].end());
    }
}

}  // namespace

int cmd_bubbles(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: snape bubbles <graph.bin> <normal> <tumour> [threads=16] [k=21] "
                             "[w=22] [maxlen=40] [max_dlen=12] [max_unmatched=20]\n");
        return 2;
    }
    const std::string gpath = argv[1], npath = argv[2], tpath = argv[3];
    const int NT = argc > 4 ? std::atoi(argv[4]) : 16;
    const int KMER = argc > 5 ? std::atoi(argv[5]) : 21;
    const int W = argc > 6 ? std::atoi(argv[6]) : 22;
    const int MAXLEN = argc > 7 ? std::atoi(argv[7]) : 40;
    const int MAX_DLEN = argc > 8 ? std::atoi(argv[8]) : 12;
    const int MAX_UNM = argc > 9 ? std::atoi(argv[9]) : 20;
    if (const char* lv = std::getenv("SNAPE_LEVEL")) LEVEL = std::atoi(lv);

    const Graph g = load_graph(gpath);
    const std::vector<Bub> bubs = find_bubbles(g, MAXLEN, NT);
    if (bubs.empty()) { std::fprintf(stderr, "[bubbles] temiz bubble yok\n"); return 0; }

    Want want_s, want_g;
    for (const Bub& b : bubs) {
        for (uint32_t v : b.som_arm) want_s.add(g.sig[v]);
        for (uint32_t v : b.germ_arm) want_g.add(g.sig[v]);
    }
    std::vector<std::vector<uint64_t>> mm_s, mm_g;
    collect(tpath, want_s, NT, mm_s, KMER, W, "tumor");
    collect(npath, want_g, NT, mm_g, KMER, W, "normal");

    // Both arms the same length and the somatic arm almost entirely colinear with the germline one:
    // a local edit, not a novel adjacency. Its cores never become a contig.
    ankerl::unordered_dense::set<uint64_t> drop;
    uint64_t nbub = 0;
    std::vector<uint64_t> ms, mg;
    for (const Bub& b : bubs) {
        arm_mins(g, b.som_arm, want_s, mm_s, ms);
        arm_mins(g, b.germ_arm, want_g, mm_g, mg);
        if (ms.empty() || mg.empty() || ms.size() > 400 || mg.size() > 400) continue;
        if (std::llabs((long long)ms.size() - (long long)mg.size()) > MAX_DLEN) continue;
        if ((long long)(ms.size() - lcs(ms, mg)) > MAX_UNM) continue;
        for (uint32_t v : b.som_arm) drop.insert(g.sig[v]);
        nbub++;
    }

    const size_t cut = gpath.find_last_of('/');
    const std::string out = (cut == std::string::npos ? std::string(".") : gpath.substr(0, cut)) + "/drop.sigs";
    FILE* f = std::fopen(out.c_str(), "w"); if (!f) die("open drop.sigs");
    for (uint64_t s : drop) std::fprintf(f, "%llu\n", (unsigned long long)s);
    std::fclose(f);
    std::fprintf(stderr, "[bubbles] |dlen|<=%d & eslesmeyen<=%d -> %llu/%zu bubble, %zu cekirdek -> %s\n",
                 MAX_DLEN, MAX_UNM, (unsigned long long)nbub, bubs.size(), drop.size(), out.c_str());
    return 0;
}
