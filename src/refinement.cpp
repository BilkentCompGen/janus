// refinement: topology + counts error-correction over the (already mm-floored) core graph.
//   graph.bin / refined input  ->  refined.bin  (IDENTICAL v2 format; simplify reads refined.bin)
//
// Three passes, all judging a structure by its read count RELATIVE TO LOCAL COVERAGE (never a
// global gate — that is what over-aggressive mm did, distorting real bubbles):
//   1. TIP removal          - short dead-end core-chains whose count << the junction they hang off,
//                             and tiny low-count isolated contigs (pure error islands).
//   2. LOW-COV side-branches - short side-sprouts whose count is a small fraction of a high-count
//                             neighbour path (recurrent-error cores clinging to a real path).
//   3. BUBBLE popping        - the low, error-like arm of a short bubble. CONSERVATIVE: a BALANCED
//                             bubble is a germline het (keep both); a low TUMOR-ONLY arm above a VAF
//                             floor is a candidate somatic (keep). Only a clearly imbalanced,
//                             non-somatic minority arm is popped.
//
// Removes nodes + incident edges, compacts ids, rewrites the v2 graph. Iterates tips (removal
// exposes new tips). Counts are the per-node ncov/tcov already in graph.bin.
#include "graph.h"
#include "svp.h"
#include "unordered_dense.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

struct Graph {
    uint64_t K = 0;
    std::vector<uint8_t>  color;
    std::vector<uint32_t> ncov, tcov, cnt;            // cnt = ncov + tcov
    std::vector<uint16_t> seqlen;
    std::vector<uint8_t>  payload;
    std::vector<uint64_t> seq_off;
    std::vector<IEdge>    edges;
    // CSR adjacency (static, built once) + live degrees (decremented as nodes die)
    std::vector<uint32_t> out_off, out_nbr, in_off, in_nbr;
    std::vector<uint32_t> live_out, live_in;
    std::vector<uint8_t>  dead;
};

double envd(const char* k, double d) { const char* v = std::getenv(k); return v ? std::atof(v) : d; }
int    envi(const char* k, int d)    { const char* v = std::getenv(k); return v ? std::atoi(v) : d; }

void load(const std::string& path, Graph& g) {
    FILE* f = xopen(path, "rb");
    rd(f, &g.K, 8); const uint64_t K = g.K;
    g.color.resize(K);  rd(f, g.color.data(), K);
    g.ncov.resize(K);   rd(f, g.ncov.data(), K * 4);
    g.tcov.resize(K);   rd(f, g.tcov.data(), K * 4);
    uint64_t E; rd(f, &E, 8);
    g.edges.resize(E);  rd(f, g.edges.data(), E * sizeof(IEdge));
    uint64_t total; rd(f, &total, 8);
    g.seqlen.resize(K); rd(f, g.seqlen.data(), K * 2);
    g.payload.resize(total); rd(f, g.payload.data(), total);
    std::fclose(f);
    g.seq_off.assign(K + 1, 0);
    for (uint64_t i = 0; i < K; i++) g.seq_off[i + 1] = g.seq_off[i] + g.seqlen[i];
    g.cnt.resize(K);
    for (uint64_t i = 0; i < K; i++) g.cnt[i] = g.ncov[i] + g.tcov[i];
    std::fprintf(stderr, "[refine] loaded K=%lu E=%lu payload=%.1fGB\n", (unsigned long)K, (unsigned long)E, total / 1e9);
}

// CSR for out- and in-adjacency from the directed edge list.
void build_adj(Graph& g) {
    const uint64_t K = g.K, E = g.edges.size();
    g.out_off.assign(K + 1, 0); g.in_off.assign(K + 1, 0);
    for (const IEdge& e : g.edges) { g.out_off[e.first + 1]++; g.in_off[e.second + 1]++; }
    for (uint64_t i = 0; i < K; i++) { g.out_off[i + 1] += g.out_off[i]; g.in_off[i + 1] += g.in_off[i]; }
    g.out_nbr.resize(E); g.in_nbr.resize(E);
    std::vector<uint32_t> oc(g.out_off.begin(), g.out_off.end() - 1), ic(g.in_off.begin(), g.in_off.end() - 1);
    for (const IEdge& e : g.edges) { g.out_nbr[oc[e.first]++] = e.second; g.in_nbr[ic[e.second]++] = e.first; }
    g.live_out.resize(K); g.live_in.resize(K);
    for (uint64_t i = 0; i < K; i++) { g.live_out[i] = g.out_off[i + 1] - g.out_off[i]; g.live_in[i] = g.in_off[i + 1] - g.in_off[i]; }
    g.dead.assign(K, 0);
}

inline void kill(Graph& g, uint32_t v) {
    if (g.dead[v]) return;
    g.dead[v] = 1;
    for (uint32_t k = g.out_off[v]; k < g.out_off[v + 1]; k++) { uint32_t w = g.out_nbr[k]; if (!g.dead[w] && g.live_in[w]) g.live_in[w]--; }
    for (uint32_t k = g.in_off[v];  k < g.in_off[v + 1];  k++) { uint32_t w = g.in_nbr[k];  if (!g.dead[w] && g.live_out[w]) g.live_out[w]--; }
}
// the single live out/in neighbour (caller guarantees live degree == 1)
inline uint32_t one_out(const Graph& g, uint32_t v) { for (uint32_t k = g.out_off[v]; k < g.out_off[v+1]; k++) if (!g.dead[g.out_nbr[k]]) return g.out_nbr[k]; return UINT32_MAX; }
inline uint32_t one_in (const Graph& g, uint32_t v) { for (uint32_t k = g.in_off[v];  k < g.in_off[v+1];  k++) if (!g.dead[g.in_nbr[k]])  return g.in_nbr[k];  return UINT32_MAX; }
inline uint32_t max_out_cnt_except(const Graph& g, uint32_t j, uint32_t except) {   // local coverage of j's MAIN arm
    uint32_t m = 0; for (uint32_t k = g.out_off[j]; k < g.out_off[j+1]; k++) { uint32_t w = g.out_nbr[k]; if (!g.dead[w] && w != except) m = std::max(m, g.cnt[w]); } return m;
}
inline uint32_t max_in_cnt_except(const Graph& g, uint32_t j, uint32_t except) {
    uint32_t m = 0; for (uint32_t k = g.in_off[j]; k < g.in_off[j+1]; k++) { uint32_t w = g.in_nbr[k]; if (!g.dead[w] && w != except) m = std::max(m, g.cnt[w]); } return m;
}

// ---- Pass 1: tips (dead-end short chains, low relative to the junction) -------------------------
// Trace the linear chain back from a forward dead-end (live_out==0) to the branch it hangs off.
// Returns the chain and the junction (UINT32_MAX if none / too long / isolated). dir: +1 backward
// from forward-dead-end via in-edges, -1 forward from backward-dead-end via out-edges.
uint32_t trace_tip(const Graph& g, uint32_t v, int dir, int maxlen, std::vector<uint32_t>& chain) {
    chain.clear(); uint32_t u = v;
    for (int s = 0; s <= maxlen; s++) {
        chain.push_back(u);
        const uint32_t deg = dir > 0 ? g.live_in[u] : g.live_out[u];
        if (deg == 0) return UINT32_MAX;                            // reached a free end -> isolated contig
        if (deg != 1) return u;                                     // u is a convergence node -> tip attaches here (ambiguous; treat u as anchor)
        const uint32_t p = dir > 0 ? one_in(g, u) : one_out(g, u);
        const uint32_t pbranch = dir > 0 ? g.live_out[p] : g.live_in[p];
        if (pbranch >= 2) return p;                                 // p is the junction; chain (v..u) is the tip
        u = p;                                                      // linear interior, keep tracing
    }
    return UINT32_MAX;                                              // too long -> not a short tip
}

// nuke: aggressive mode — ignore the abs_keep floor and remove EVERY tip <= maxlen unconditionally
// (REF_TIP_NUKE). Measured cost vs the conservative default: it clears ~116k extra tip-error but
// erodes only ~0.03% of somatic cores (short SV-end trims, not lost calls) — tips carry ~no variant
// signal (28 of 4.1M removable tips before refine). Long real tips (> maxlen) are untouched either way.
uint64_t pass_tips(Graph& g, int maxlen, double frac, uint32_t island_cnt, uint32_t abs_keep, int rounds, bool nuke) {
    uint64_t removed = 0;
    std::vector<uint32_t> chain; std::vector<uint32_t> tokill;
    for (int r = 0; r < rounds; r++) {
        tokill.clear();
        for (uint64_t v = 0; v < g.K; v++) {
            if (g.dead[v]) continue;
            for (int dir = +1; dir >= -1; dir -= 2) {               // +1: forward dead-end, -1: backward dead-end
                const uint32_t end_deg = dir > 0 ? g.live_out[v] : g.live_in[v];
                if (end_deg != 0) continue;
                const uint32_t j = trace_tip(g, (uint32_t)v, dir, maxlen, chain);
                if (chain.size() > (size_t)maxlen) continue;
                uint32_t tipmax = 0; for (uint32_t c : chain) tipmax = std::max(tipmax, g.cnt[c]);
                if (!nuke && tipmax >= abs_keep) continue;          // a well-covered tip is ~99.6% real (count>=10) -> never drop (off in nuke)
                bool drop = nuke;                                   // nuke: remove every tip<=maxlen unconditionally
                if (j == UINT32_MAX) {                              // isolated short contig -> error island if low count
                    if (tipmax < island_cnt) drop = true;
                } else {                                            // tip hanging off junction j; base = chain.back()
                    const uint32_t base = chain.back();
                    const uint32_t main = dir > 0 ? max_out_cnt_except(g, j, base) : max_in_cnt_except(g, j, base);
                    if (main && tipmax < (uint32_t)(frac * main)) drop = true;
                }
                if (drop) for (uint32_t c : chain) tokill.push_back(c);
            }
        }
        if (tokill.empty()) break;
        for (uint32_t c : tokill) if (!g.dead[c]) { kill(g, c); removed++; }
    }
    std::fprintf(stderr, "[refine] pass 1 (tips): removed %lu nodes\n", (unsigned long)removed);
    return removed;
}

// ---- Pass 2: low-cov side-branches (a node far below its highest-count neighbour) ----------------
// A core whose count is a tiny fraction of its strongest neighbour AND that is a leaf-ish branch
// (in- or out-degree <= 1, so it is a sprout off a real path, not a hub on it). Local-relative,
// so genuine low-coverage regions (neighbours also low) are spared.
uint64_t pass_lowcov(Graph& g, double frac, uint32_t abs_keep) {
    uint64_t removed = 0;
    std::vector<uint32_t> tokill;
    for (uint64_t v = 0; v < g.K; v++) {
        if (g.dead[v]) continue;
        if (g.cnt[v] >= abs_keep) continue;                        // well-covered -> ~99.6% real -> keep
        if (g.live_in[v] > 1 && g.live_out[v] > 1) continue;       // hub on a path -> don't touch
        uint32_t nmax = 0;
        for (uint32_t k = g.out_off[v]; k < g.out_off[v+1]; k++) { uint32_t w = g.out_nbr[k]; if (!g.dead[w]) nmax = std::max(nmax, g.cnt[w]); }
        for (uint32_t k = g.in_off[v];  k < g.in_off[v+1];  k++) { uint32_t w = g.in_nbr[k];  if (!g.dead[w]) nmax = std::max(nmax, g.cnt[w]); }
        if (nmax && g.cnt[v] < (uint32_t)(frac * nmax)) tokill.push_back((uint32_t)v);
    }
    for (uint32_t c : tokill) if (!g.dead[c]) { kill(g, c); removed++; }
    std::fprintf(stderr, "[refine] pass 2 (low-cov branches): removed %lu nodes\n", (unsigned long)removed);
    return removed;
}

// ---- Pass 3: bubble popping (conservative; spares germline-het and candidate-somatic arms) -------
// Trace the linear arm forward from a successor of branch S until it converges (live_in>=2) or ends.
uint32_t trace_arm(const Graph& g, uint32_t start, int maxlen, std::vector<uint32_t>& arm) {
    arm.clear(); uint32_t cur = start;
    for (int s = 0; s <= maxlen; s++) {
        arm.push_back(cur);
        if (g.live_out[cur] != 1) return UINT32_MAX;
        const uint32_t nxt = one_out(g, cur);
        if (g.live_in[nxt] >= 2) return nxt;                       // convergence node T
        if (g.live_in[nxt] != 1) return UINT32_MAX;
        cur = nxt;
    }
    return UINT32_MAX;
}

uint64_t pass_bubbles(Graph& g, int maxlen, double het_frac, double som_frac, uint32_t abs_keep) {
    uint64_t removed = 0;
    std::vector<uint32_t> a1, a2, tokill;
    for (uint64_t s = 0; s < g.K; s++) {
        if (g.dead[s] || g.live_out[s] < 2) continue;
        // collect (T, arm, mean-count, tumor-only?) for each live successor
        std::vector<uint32_t> succ;
        for (uint32_t k = g.out_off[s]; k < g.out_off[s+1]; k++) { uint32_t w = g.out_nbr[k]; if (!g.dead[w]) succ.push_back(w); }
        for (size_t i = 0; i < succ.size(); i++) for (size_t jx = i + 1; jx < succ.size(); jx++) {
            const uint32_t T1 = trace_arm(g, succ[i], maxlen, a1); if (T1 == UINT32_MAX) continue;
            const uint32_t T2 = trace_arm(g, succ[jx], maxlen, a2); if (T2 == UINT32_MAX || T2 != T1) continue;
            // two arms s -> ... -> T : a bubble. mean count per arm, and is the LOW arm tumor-only?
            auto mean = [&](std::vector<uint32_t>& a){ uint64_t s2 = 0; for (uint32_t c : a) s2 += g.cnt[c]; return (uint32_t)(s2 / std::max<size_t>(1, a.size())); };
            uint32_t m1 = mean(a1), m2 = mean(a2);
            std::vector<uint32_t>& lo = (m1 <= m2) ? a1 : a2;       // the minority (lower-coverage) arm
            if (std::max(m1, m2) == 0) continue;
            // Per-colour, FLANK-RELATIVE classification (Cortex genotyping + Severus somatic logic): judge
            // the minority arm against the bubble's anchor (flank) coverage in EACH sample, not an absolute
            // count. flank = the coverage the arm SHOULD have if it were a real allele there.
            const double fn = std::min(g.ncov[s], g.ncov[T1]);      // normal flank (anchor) coverage
            const double ft = std::min(g.tcov[s], g.tcov[T1]);      // tumor  flank coverage
            double ln = 0, lt = 0; for (uint32_t c : lo) { ln += g.ncov[c]; lt += g.tcov[c]; }
            ln /= lo.size(); lt /= lo.size();
            const bool germline = fn > 0 && ln >= het_frac * fn;    // germline-level in normal -> real allele
            const bool somatic  = ft > 0 && lt >= som_frac * ft && (fn == 0 || ln < som_frac * fn);  // present in tumor AND absent in normal
            if (germline || somatic) continue;                      // real (germline-het or low-VAF somatic) -> KEEP
            uint32_t lo_max = 0; for (uint32_t c : lo) lo_max = std::max(lo_max, g.cnt[c]);
            if (lo_max >= abs_keep) continue;                       // well-covered arm: indistinguishable real-vs-error by coverage -> KEEP (safe)
            for (uint32_t c : lo) tokill.push_back(c);              // low-count, sub-het, non-somatic minority arm = error -> POP
        }
    }
    for (uint32_t c : tokill) if (!g.dead[c]) { kill(g, c); removed++; }
    std::fprintf(stderr, "[refine] pass 3 (bubbles): removed %lu nodes\n", (unsigned long)removed);
    return removed;
}

void write_refined(const std::string& path, const Graph& g) {
    const uint64_t K = g.K;
    std::vector<uint32_t> newid(K, UINT32_MAX); uint64_t nk = 0;
    for (uint64_t i = 0; i < K; i++) if (!g.dead[i]) newid[i] = (uint32_t)nk++;
    std::vector<uint8_t> color; std::vector<uint32_t> ncov, tcov; std::vector<uint16_t> seqlen;
    color.reserve(nk); ncov.reserve(nk); tcov.reserve(nk); seqlen.reserve(nk);
    uint64_t total = 0;
    for (uint64_t i = 0; i < K; i++) if (!g.dead[i]) { color.push_back(g.color[i]); ncov.push_back(g.ncov[i]); tcov.push_back(g.tcov[i]); seqlen.push_back(g.seqlen[i]); total += g.seqlen[i]; }
    std::vector<IEdge> ie; ie.reserve(g.edges.size());
    for (const IEdge& e : g.edges) if (!g.dead[e.first] && !g.dead[e.second]) ie.push_back(IEdge(newid[e.first], newid[e.second]));
    std::sort(ie.begin(), ie.end()); ie.erase(std::unique(ie.begin(), ie.end()), ie.end());

    FILE* f = xopen(path + ".tmp", "wb");
    static char wbuf[1u << 24]; std::setvbuf(f, wbuf, _IOFBF, sizeof(wbuf));
    uint64_t Kw = nk; wr(f, &Kw, 8);
    wr(f, color.data(), nk); wr(f, ncov.data(), nk * 4); wr(f, tcov.data(), nk * 4);
    uint64_t E = ie.size(); wr(f, &E, 8); wr(f, ie.data(), E * sizeof(IEdge));
    wr(f, &total, 8); wr(f, seqlen.data(), nk * 2);
    for (uint64_t i = 0; i < K; i++) if (!g.dead[i]) wr(f, g.payload.data() + g.seq_off[i], g.seqlen[i]);
    atomic_close(f, path + ".tmp", path);
    std::fprintf(stderr, "[refine] wrote %s: nodes %lu -> %lu  edges %lu -> %lu  payload %.1fGB\n",
                 path.c_str(), (unsigned long)K, (unsigned long)nk, (unsigned long)g.edges.size(), (unsigned long)E, total / 1e9);
}

}  // namespace

// graph.bin -> refined.bin (same v2 format). Tunable via env: REF_TIP_MAX, REF_TIP_FRAC,
// REF_ISLAND_CNT, REF_LOWCOV_FRAC, REF_BUB_MAX, REF_HET_FRAC, REF_SOM_FRAC, REF_ABS_KEEP, REF_ROUNDS.
// REF_TIP_NUKE=1 -> aggressive tip removal (see pass_tips); writes a DISTINCT dir/refined_nuketip.bin
// so the conservative dir/refined.bin is never clobbered.
void run_refine(const std::string& dir) {
    const bool   nuke        = envi("REF_TIP_NUKE", 0);
    const std::string in = dir + "/graph.bin";
    const std::string out = dir + (nuke ? "/refined_nuketip.bin" : "/refined.bin");
    const int    tip_max     = envi("REF_TIP_MAX", 12);
    const double tip_frac    = envd("REF_TIP_FRAC", 0.34);
    const int    island_cnt  = envi("REF_ISLAND_CNT", 8);
    const double lowcov_frac = envd("REF_LOWCOV_FRAC", 0.12);
    const int    bub_max     = envi("REF_BUB_MAX", 16);
    const double het_frac    = envd("REF_HET_FRAC", 0.25);          // bubble arm >= this x normal-flank -> germline allele -> keep
    const double som_frac    = envd("REF_SOM_FRAC", 0.05);          // arm >= this x tumor-flank AND normal-absent -> somatic -> keep
    const uint32_t abs_keep  = (uint32_t)envi("REF_ABS_KEEP", 10);   // tips/low-cov only: count>=this never removed (bubbles use flank-relative)
    const int    rounds      = envi("REF_ROUNDS", 4);
    std::fprintf(stderr, "[refine] params: tip_max=%d tip_frac=%.2f island_cnt=%d lowcov_frac=%.2f bub_max=%d het_frac=%.2f som_frac=%.2f abs_keep=%u rounds=%d nuke=%d -> %s\n",
                 tip_max, tip_frac, island_cnt, lowcov_frac, bub_max, het_frac, som_frac, abs_keep, rounds, (int)nuke, out.c_str());
    Graph g; load(in, g); build_adj(g);
    uint64_t r1 = pass_tips(g, tip_max, tip_frac, (uint32_t)island_cnt, abs_keep, rounds, nuke);
    uint64_t r2 = pass_lowcov(g, lowcov_frac, abs_keep);
    uint64_t r3 = pass_bubbles(g, bub_max, het_frac, som_frac, abs_keep);
    pass_tips(g, tip_max, tip_frac, (uint32_t)island_cnt, abs_keep, rounds, nuke);   // re-run tips: bubble/branch removal exposes new tips
    std::fprintf(stderr, "[refine] total removed: tips=%lu lowcov=%lu bubbles=%lu (+ tip re-pass)\n",
                 (unsigned long)r1, (unsigned long)r2, (unsigned long)r3);
    write_refined(out, g);
}
