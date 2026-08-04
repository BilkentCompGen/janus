// Topology + counts error-correction: <dir>/graph.bin -> <dir>/refined.bin, same v3 format.
// What each pass does and why it is wired in: src/refinement.h.
#include "refinement.h"

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
    std::vector<uint64_t> sig;                        // murmur node identity: carried through unchanged, never read by a pass
    std::vector<uint32_t> ncov, tcov, cnt;            // cnt = ncov + tcov
    std::vector<CEdge>    edges;                       // directed src->dst with per-colour traversal counts
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
    g.sig.resize(K);    rd(f, g.sig.data(), K * 8);
    g.ncov.resize(K);   rd(f, g.ncov.data(), K * 4);
    g.tcov.resize(K);   rd(f, g.tcov.data(), K * 4);
    uint64_t E; rd(f, &E, 8);
    g.edges.resize(E);  rd(f, g.edges.data(), E * sizeof(CEdge));
    std::fclose(f);
    g.cnt.resize(K);
    for (uint64_t i = 0; i < K; i++) g.cnt[i] = g.ncov[i] + g.tcov[i];
    std::fprintf(stderr, "[refine] loaded K=%lu E=%lu (v3: sig + edge counts, no payload)\n", (unsigned long)K, (unsigned long)E);
}

// CSR for out- and in-adjacency from the directed edge list (endpoints only; per-edge counts stay
// in g.edges and are carried straight through to the output -- the passes judge nodes, not edges).
void build_adj(Graph& g) {
    const uint64_t K = g.K, E = g.edges.size();
    g.out_off.assign(K + 1, 0); g.in_off.assign(K + 1, 0);
    for (const CEdge& e : g.edges) { g.out_off[e.src + 1]++; g.in_off[e.dst + 1]++; }
    for (uint64_t i = 0; i < K; i++) { g.out_off[i + 1] += g.out_off[i]; g.in_off[i + 1] += g.in_off[i]; }
    g.out_nbr.resize(E); g.in_nbr.resize(E);
    std::vector<uint32_t> oc(g.out_off.begin(), g.out_off.end() - 1), ic(g.in_off.begin(), g.in_off.end() - 1);
    for (const CEdge& e : g.edges) { g.out_nbr[oc[e.src]++] = e.dst; g.in_nbr[ic[e.dst]++] = e.src; }
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

// Rebuild CSR + live degrees from the CURRENT edge list, honoring dead[] (does NOT reset it, unlike
// build_adj). Edges with a dead endpoint are dropped from g.edges first. Used between iterative-peel
// rounds so the re-run node passes see the post-deflation topology. cnt is refreshed by the caller.
void rebuild_live_adj(Graph& g) {
    std::vector<CEdge> live; live.reserve(g.edges.size());
    for (const CEdge& e : g.edges) if (!g.dead[e.src] && !g.dead[e.dst]) live.push_back(e);
    g.edges.swap(live);
    const uint64_t K = g.K, E = g.edges.size();
    g.out_off.assign(K + 1, 0); g.in_off.assign(K + 1, 0);
    for (const CEdge& e : g.edges) { g.out_off[e.src + 1]++; g.in_off[e.dst + 1]++; }
    for (uint64_t i = 0; i < K; i++) { g.out_off[i + 1] += g.out_off[i]; g.in_off[i + 1] += g.in_off[i]; }
    g.out_nbr.resize(E); g.in_nbr.resize(E);
    std::vector<uint32_t> oc(g.out_off.begin(), g.out_off.end() - 1), ic(g.in_off.begin(), g.in_off.end() - 1);
    for (const CEdge& e : g.edges) { g.out_nbr[oc[e.src]++] = e.dst; g.in_nbr[ic[e.dst]++] = e.src; }
    g.live_out.resize(K); g.live_in.resize(K);
    for (uint64_t i = 0; i < K; i++) { g.live_out[i] = g.out_off[i+1]-g.out_off[i]; g.live_in[i] = g.in_off[i+1]-g.in_off[i]; }
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

// ---- Pass 4: edge filtering (the node passes judge nodes; edges are gated here) ------------------
// Node counts are floored to 0-or->=mm by the build-time gate, but edge traversal counts are raw
// (min 1), so a "tumor-only" edge can rest on a single tumor read -- the dominant somatic-junction
// false positive. Two gates, both per-colour and local (never a global cut):
//   (1) HARD  : drop an edge unless it is crossed >= mm times in at least one colour (the edge
//               analogue of the node mm floor). max(ncov,tcov) >= mm.
//   (2) RELATIVE: even above the floor, drop an edge that carries only a small fraction (< frac) of
//               its endpoints' coverage in EVERY colour it has -- a weak spur between two cores that
//               are themselves much better covered (reads mostly take a different path). frac=0 off.
// Only edges between two live nodes are gated; edges with a killed endpoint are dropped regardless.
uint64_t pass_edges(Graph& g, uint32_t mm, double frac) {
    std::vector<CEdge> keep; keep.reserve(g.edges.size());
    uint64_t d_hard = 0, d_rel = 0, d_dead = 0;
    for (const CEdge& e : g.edges) {
        if (g.dead[e.src] || g.dead[e.dst]) { d_dead++; continue; }        // endpoint killed by a node pass
        if (e.ncov < mm && e.tcov < mm) { d_hard++; continue; }            // below the mm floor in both colours
        if (frac > 0.0) {
            bool strong = false;                                           // strong in >= 1 colour it carries -> keep
            if (e.tcov) { uint32_t m = std::min(g.tcov[e.src], g.tcov[e.dst]); if (!m || (double)e.tcov >= frac * m) strong = true; }
            if (e.ncov) { uint32_t m = std::min(g.ncov[e.src], g.ncov[e.dst]); if (!m || (double)e.ncov >= frac * m) strong = true; }
            if (!strong) { d_rel++; continue; }
        }
        keep.push_back(e);
    }
    g.edges.swap(keep);
    std::fprintf(stderr, "[refine] pass 4 (edges): dropped hard(mm=%u)=%lu relative(frac=%.2f)=%lu (dead-endpoint=%lu) -> %lu edges\n",
                 mm, (unsigned long)d_hard, frac, (unsigned long)d_rel, (unsigned long)d_dead, (unsigned long)g.edges.size());
    return d_hard + d_rel;
}

// ---- iterative-peel edge pass: gate edges AND deflate node counts by every removed edge ---------
// The idea: an error edge's traversal reads also inflated its endpoint nodes' coverage, so when we
// drop the edge we subtract its per-colour counts from the endpoint node(s). Gate-dropped edges (hard
// mm + relative) deflate BOTH endpoints; edges whose endpoint already died in a node pass deflate the
// lone survivor. Deflation is accumulated and applied AFTER the pass, so the relative gate reads
// consistent start-of-round counts (order-independent). Returns gate-dropped edges (drives convergence).
// CAVEAT: a read crosses a node once but up to two edges, so a node with one kept + one dropped edge
// is over-deflated -- an intentional approximation of the requested rule; the measurement shows its cost.
uint64_t pass_edges_deflate(Graph& g, uint32_t mm, double frac, bool clamp) {
    const uint64_t K = g.K;
    std::vector<uint32_t> dn(K, 0), dt(K, 0);
    std::vector<CEdge> keep; keep.reserve(g.edges.size());
    uint64_t d_gate = 0, d_dead = 0;
    for (const CEdge& e : g.edges) {
        const bool sd = g.dead[e.src], dd = g.dead[e.dst];
        if (sd && dd) { d_dead++; continue; }                              // both gone: nothing left to deflate
        if (sd || dd) { const uint32_t s = sd ? e.dst : e.src; dn[s] += e.ncov; dt[s] += e.tcov; d_dead++; continue; }
        bool drop = (e.ncov < mm && e.tcov < mm);                          // hard floor
        if (!drop && frac > 0.0) {
            bool strong = false;
            if (e.tcov) { uint32_t m = std::min(g.tcov[e.src], g.tcov[e.dst]); if (!m || (double)e.tcov >= frac * m) strong = true; }
            if (e.ncov) { uint32_t m = std::min(g.ncov[e.src], g.ncov[e.dst]); if (!m || (double)e.ncov >= frac * m) strong = true; }
            drop = !strong;
        }
        if (drop) { dn[e.src] += e.ncov; dt[e.src] += e.tcov; dn[e.dst] += e.ncov; dt[e.dst] += e.tcov; d_gate++; }
        else keep.push_back(e);
    }
    g.edges.swap(keep);
    for (uint64_t i = 0; i < K; i++) { g.ncov[i] = dn[i] < g.ncov[i] ? g.ncov[i] - dn[i] : 0; g.tcov[i] = dt[i] < g.tcov[i] ? g.tcov[i] - dt[i] : 0; }
    // Reverse coupling (opt-in, REF_EDGE_CLAMP): an edge can't carry more reads than either endpoint,
    // so clamp each surviving edge to min(endpoint counts). MEASURED TO HURT edges: when deflation zeros
    // a node's ncov, the clamp zeros its edges' ncov too -> shared edges recolour to tumor-only -> ~55k
    // false-somatic edges. Off by default; the flow-conservation filter is the safe alternative (it drops
    // over-flow edges instead of recolouring them).
    uint64_t clamped = 0;
    if (clamp) for (CEdge& e : g.edges) {
        const uint32_t cn = std::min(g.ncov[e.src], g.ncov[e.dst]);
        const uint32_t ct = std::min(g.tcov[e.src], g.tcov[e.dst]);
        if (e.ncov > cn) { e.ncov = cn; clamped++; }
        if (e.tcov > ct) { e.tcov = ct; clamped++; }
    }
    std::fprintf(stderr, "[refine]   peel edges: gate-dropped=%lu node-death-dropped=%lu edge-clamped=%lu -> %lu edges kept\n",
                 (unsigned long)d_gate, (unsigned long)d_dead, (unsigned long)clamped, (unsigned long)g.edges.size());
    return d_gate + clamped;   // total edge activity (gate drops + clamp reductions) -> drives convergence
}

// ---- flow-conservation node filter -------------------------------------------------------------
// Every read through a node is partitioned among its in-edges (and among its out-edges), so in EACH
// colour Sum(in-edge counts) <= node.count and Sum(out-edge counts) <= node.count. A node that
// violates this claims more edge-traversals than it has reads -- an inconsistency that is empirically
// ~95% a read-error core (real cores violate ~0.01%). Remove violators. MUST run on CONSISTENT counts
// (before the peel deflates node counts below their still-original edge counts), so it is applied once
// after the initial node refine, ahead of the peel loop.
uint64_t flow_filter(Graph& g) {
    const uint64_t K = g.K;
    std::vector<uint32_t> in_n(K, 0), in_t(K, 0), out_n(K, 0), out_t(K, 0);
    for (const CEdge& e : g.edges) {
        if (g.dead[e.src] || g.dead[e.dst]) continue;
        out_n[e.src] += e.ncov; out_t[e.src] += e.tcov;
        in_n[e.dst]  += e.ncov; in_t[e.dst]  += e.tcov;
    }
    uint64_t removed = 0;
    for (uint64_t i = 0; i < K; i++)
        if (!g.dead[i] && (in_n[i] > g.ncov[i] || out_n[i] > g.ncov[i] || in_t[i] > g.tcov[i] || out_t[i] > g.tcov[i])) { g.dead[i] = 1; removed++; }
    std::fprintf(stderr, "[refine] flow-conservation filter: removed %lu nodes (Sum incident edges > node count)\n", (unsigned long)removed);
    return removed;
}

void write_refined(const std::string& path, const Graph& g) {
    const uint64_t K = g.K;
    std::vector<uint32_t> newid(K, UINT32_MAX); uint64_t nk = 0;
    for (uint64_t i = 0; i < K; i++) if (!g.dead[i]) newid[i] = (uint32_t)nk++;
    std::vector<uint64_t> sig; std::vector<uint32_t> ncov, tcov;
    sig.reserve(nk); ncov.reserve(nk); tcov.reserve(nk);
    for (uint64_t i = 0; i < K; i++) if (!g.dead[i]) { sig.push_back(g.sig[i]); ncov.push_back(g.ncov[i]); tcov.push_back(g.tcov[i]); }
    std::vector<CEdge> ce; ce.reserve(g.edges.size());
    for (const CEdge& e : g.edges) if (!g.dead[e.src] && !g.dead[e.dst]) ce.push_back({newid[e.src], newid[e.dst], e.ncov, e.tcov});
    std::sort(ce.begin(), ce.end(), [](const CEdge& x, const CEdge& y){ return x.src != y.src ? x.src < y.src : x.dst < y.dst; });

    FILE* f = xopen(path + ".tmp", "wb");
    static char wbuf[1u << 24]; std::setvbuf(f, wbuf, _IOFBF, sizeof(wbuf));
    uint64_t Kw = nk; wr(f, &Kw, 8);
    wr(f, sig.data(), nk * 8); wr(f, ncov.data(), nk * 4); wr(f, tcov.data(), nk * 4);
    uint64_t E = ce.size(); wr(f, &E, 8); wr(f, ce.data(), E * sizeof(CEdge));
    atomic_close(f, path + ".tmp", path);
    std::fprintf(stderr, "[refine] wrote %s: nodes %lu -> %lu  edges %lu -> %lu (v3)\n",
                 path.c_str(), (unsigned long)K, (unsigned long)nk, (unsigned long)g.edges.size(), (unsigned long)E);
}

}  // namespace

// graph.bin -> refined.bin (same v3 format). Tunable via env: REF_TIP_MAX, REF_TIP_FRAC,
// REF_ISLAND_CNT, REF_LOWCOV_FRAC, REF_BUB_MAX, REF_HET_FRAC, REF_SOM_FRAC, REF_ABS_KEEP, REF_ROUNDS,
// REF_EDGE_MM (edge hard mm floor; 0 disables), REF_EDGE_FRAC (edge relative gate; 0 disables).
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
    const uint32_t edge_mm   = (uint32_t)envi("REF_EDGE_MM", 3);     // edge hard floor: keep only if max(ncov,tcov) >= this
    const double edge_frac   = envd("REF_EDGE_FRAC", 0.10);          // edge relative gate: drop if < this x min-endpoint coverage in every colour
    const int    iter_n      = envi("REF_ITER", 0);                  // >0: iterative peel (edge-drop -> deflate node counts -> re-refine); value = MAX rounds, stops early at convergence; writes refined_iter.bin
    const bool   edge_clamp  = envi("REF_EDGE_CLAMP", 0);            // clamp edges to min(endpoint counts) each peel round -- OFF: measured to recolour shared edges -> +FP
    const bool   flow_filt   = envi("REF_FLOW", 1);                  // flow-conservation filter before the peel: drop nodes whose summed incident edges exceed the node count (~95% error)
    std::fprintf(stderr, "[refine] params: tip_max=%d tip_frac=%.2f island_cnt=%d lowcov_frac=%.2f bub_max=%d het_frac=%.2f som_frac=%.2f abs_keep=%u rounds=%d edge_mm=%u edge_frac=%.2f nuke=%d -> %s\n",
                 tip_max, tip_frac, island_cnt, lowcov_frac, bub_max, het_frac, som_frac, abs_keep, rounds, edge_mm, edge_frac, (int)nuke, out.c_str());
    Graph g; load(in, g); build_adj(g);
    uint64_t r1 = pass_tips(g, tip_max, tip_frac, (uint32_t)island_cnt, abs_keep, rounds, nuke);
    uint64_t r2 = pass_lowcov(g, lowcov_frac, abs_keep);
    uint64_t r3 = pass_bubbles(g, bub_max, het_frac, som_frac, abs_keep);
    pass_tips(g, tip_max, tip_frac, (uint32_t)island_cnt, abs_keep, rounds, nuke);   // re-run tips: bubble/branch removal exposes new tips
    std::fprintf(stderr, "[refine] total removed: tips=%lu lowcov=%lu bubbles=%lu (+ tip re-pass)\n",
                 (unsigned long)r1, (unsigned long)r2, (unsigned long)r3);

    if (iter_n > 0) {
        // ---- iterative peel: drop error edges, subtract their reads from the endpoint nodes, then
        // re-refine the (now deflated) nodes; removing nodes exposes more dead edges next round. The
        // peel is monotonic (dead[] only grows, edges only shrink) so it always converges; REF_ITER is
        // just the safety bound on rounds -- the ed==0 && nd==0 check below stops at true convergence. ----
        const int N = iter_n;
        const std::string iout = dir + "/refined_iter.bin";
        std::fprintf(stderr, "[refine] ITERATIVE PEEL: edge-drop -> deflate node counts -> re-refine, up to %d rounds (edge_mm=%u edge_frac=%.2f clamp=%d flow=%d) -> %s\n",
                     N, edge_mm, edge_frac, (int)edge_clamp, (int)flow_filt, iout.c_str());
        if (flow_filt) flow_filter(g);   // one-shot error removal on consistent counts, before deflation distorts them
        for (int it = 1; it <= N; it++) {
            std::fprintf(stderr, "[refine] === peel round %d/%d ===\n", it, N);
            const uint64_t ea = pass_edges_deflate(g, edge_mm, edge_frac, edge_clamp);   // gate + deflate nodes (+ clamp edges if enabled)
            // A node with ZERO coverage in BOTH samples has no read support left -> remove it outright
            // (deflation only zeroed the depleted colour; single-colour-zero nodes stay, just recoloured).
            uint64_t z = 0;
            for (uint64_t i = 0; i < g.K; i++) if (!g.dead[i] && g.ncov[i] == 0 && g.tcov[i] == 0) { g.dead[i] = 1; z++; }
            for (uint64_t i = 0; i < g.K; i++) g.cnt[i] = g.ncov[i] + g.tcov[i];    // refresh cnt after deflation
            rebuild_live_adj(g);
            uint64_t nd = pass_tips(g, tip_max, tip_frac, (uint32_t)island_cnt, abs_keep, rounds, nuke);
            nd += pass_lowcov(g, lowcov_frac, abs_keep);
            nd += pass_bubbles(g, bub_max, het_frac, som_frac, abs_keep);
            nd += pass_tips(g, tip_max, tip_frac, (uint32_t)island_cnt, abs_keep, rounds, nuke);
            std::fprintf(stderr, "[refine] round %d done: edge-activity(gate+clamp)=%lu zero-cov-removed=%lu node-pass-removed=%lu\n",
                         it, (unsigned long)ea, (unsigned long)z, (unsigned long)nd);
            if (ea == 0 && z == 0 && nd == 0) { std::fprintf(stderr, "[refine] converged at round %d\n", it); break; }
        }
        write_refined(iout, g);
        return;
    }

    if (edge_mm > 1 || edge_frac > 0.0) pass_edges(g, edge_mm, edge_frac);
    write_refined(out, g);
}
