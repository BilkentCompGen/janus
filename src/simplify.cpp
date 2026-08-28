// simplify: the single graph-simplification stage. It builds the compacted super-graph in RAM
// (unitig compaction; see build_super_graph in compact.cpp), then collapses it one step
// further and emits the WHOLE graph — not extracted parts. Every compacted super-node passes
// through; for each std_both (S1,S2) bubble the prefix/suffix shared across ALL of an
// anchor's bubble arms is absorbed into the anchor (grown), leaving each arm as just its
// divergent middle at base-pair resolution. Coverage is the per-base (n,t) track,
// max-merged on shared regions.
//
// Output schema = the compacted node schema (SuperNode{block_len,runs_len,colour} + blocks
// + runs) PLUS edges that carry an explicit overlap: reconstruct a path by dropping
// `overlap` bytes from the head of each `to` node before concatenating. Bubble-seam edges
// (grown-anchor -> residual -> anchor) have overlap 0 — the arms diverge at exactly one
// base, so no overlap can span it; pass-through edges keep the compaction suffix/prefix
// overlap. An anchor is grown only when every edge on the grown side is a zipped arm
// (otherwise growing its tail would corrupt a surviving backbone edge's overlap).
#include "graph.h"
#include "svp.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

// KMP-based longest k such that A[la-k:la] == B[0:k]. The fall-back at i==la-1 only
// triggers if B fully matched as a suffix of A — that case must return lb, not fail[lb-1].
uint32_t suffix_prefix(const uint8_t* A, uint32_t la, const uint8_t* B, uint32_t lb) {
    if (la == 0 || lb == 0) return 0;
    std::vector<uint32_t> fail(lb, 0);
    for (uint32_t i = 1, j = 0; i < lb; i++) {
        while (j > 0 && B[i] != B[j]) j = fail[j - 1];
        if (B[i] == B[j]) j++;
        fail[i] = j;
    }
    uint32_t j = 0;
    for (uint32_t i = 0; i < la; i++) {
        while (j > 0 && A[i] != B[j]) j = fail[j - 1];
        if (A[i] == B[j]) j++;
        if (j == lb) { if (i + 1 == la) return lb; j = fail[j - 1]; }
    }
    return j;
}

namespace {

void walk_chain(uint32_t u, const std::vector<uint32_t>& succ,
                const std::vector<uint32_t>& in_deg, const std::vector<uint32_t>& out_deg,
                const std::vector<uint8_t>& color, std::vector<char>& visited,
                std::vector<uint32_t>& chain) {
    chain.clear();
    chain.push_back(u);
    visited[u] = 1;
    uint32_t curr = u;
    while (out_deg[curr] == 1) {
        uint32_t next = succ[curr];
        if (in_deg[next] != 1 || color[next] != color[curr] || visited[next]) break;
        chain.push_back(next);
        visited[next] = 1;
        curr = next;
    }
}

// k == 0 is FATAL — a chain edge with zero byte overlap means upstream data is broken
// (LCP-derived adjacent windows always share ~85 bp at L5; see [[feedback_lcp_core_overlap]]).
void splice_chain(const std::vector<uint32_t>& chain,
                  const std::vector<uint16_t>& seqlen, const uint8_t* payload,
                  const std::vector<uint64_t>& seq_off,
                  std::vector<uint8_t>& block, std::vector<uint16_t>& overlap) {
    block.clear(); overlap.clear();
    const uint32_t l0 = seqlen[chain[0]];
    block.insert(block.end(), payload + seq_off[chain[0]], payload + seq_off[chain[0]] + l0);
    for (size_t i = 1; i < chain.size(); i++) {
        const uint16_t li = seqlen[chain[i]];
        const uint8_t* B = payload + seq_off[chain[i]];
        const uint16_t lt = seqlen[chain[i - 1]];
        const uint8_t* A = block.data() + (block.size() - lt);
        const uint32_t k = suffix_prefix(A, lt, B, li);
        if (k == 0) { std::fprintf(stderr, "[FATAL] compact: zero-overlap chain edge node=%u->%u\n", chain[i-1], chain[i]); std::abort(); }
        overlap.push_back((uint16_t)k);
        block.insert(block.end(), B + k, B + li);
    }
}

// Sweep p over [0, block_len). [l, r) tracks the active chain range whose [seg_start,
// seg_end) interval contains p. Emit a CovRun whenever max(n_cov, t_cov) changes — RLE.
void compute_runs(const std::vector<uint32_t>& chain,
                  const std::vector<uint16_t>& seqlen,
                  const std::vector<uint32_t>& ncov, const std::vector<uint32_t>& tcov,
                  const std::vector<uint16_t>& overlap, uint64_t block_len,
                  std::vector<CovRun>& runs) {
    runs.clear();
    const size_t m = chain.size();
    std::vector<uint32_t> seg_start(m, 0);
    for (size_t i = 1; i < m; i++) seg_start[i] = seg_start[i - 1] + seqlen[chain[i - 1]] - overlap[i - 1];
    size_t l = 0, r = 0;
    uint32_t cur_n = 0, cur_t = 0;
    bool started = false;
    for (uint32_t p = 0; p < (uint32_t)block_len; p++) {
        while (r < m && seg_start[r] == p) r++;
        while (l < r && seg_start[l] + seqlen[chain[l]] == p) l++;
        uint32_t n_max = 0, t_max = 0;
        for (size_t i = l; i < r; i++) {
            const uint32_t nc = ncov[chain[i]], tc = tcov[chain[i]];
            if (nc > n_max) n_max = nc;
            if (tc > t_max) t_max = tc;
        }
        if (!started) { cur_n = n_max; cur_t = t_max; started = true; continue; }
        if (n_max != cur_n || t_max != cur_t) {
            runs.push_back({p, cur_n, cur_t});
            cur_n = n_max; cur_t = t_max;
        }
    }
    runs.push_back({(uint32_t)block_len, cur_n, cur_t});
}

void build_super_edges(const std::vector<IEdge>& edges,
                       const std::vector<uint32_t>& node_to_super,
                       std::vector<IEdge>& super_edges) {
    super_edges.clear();
    super_edges.reserve(edges.size());
    for (const IEdge& e : edges) {
        const uint32_t su = node_to_super[e.first], sv = node_to_super[e.second];
        if (su != sv) super_edges.push_back({su, sv});
    }
    std::sort(super_edges.begin(), super_edges.end());
    super_edges.erase(std::unique(super_edges.begin(), super_edges.end()), super_edges.end());
}

void emit_super_node(const std::vector<uint32_t>& chain, const std::vector<uint16_t>& seqlen,
                     const uint8_t* payload, const std::vector<uint64_t>& seq_off,
                     const std::vector<uint32_t>& ncov, const std::vector<uint32_t>& tcov, const std::vector<uint8_t>& color,
                     std::vector<SuperNode>& snodes,
                     std::vector<uint8_t>& blocks, std::vector<CovRun>& runs, std::vector<uint32_t>& node_to_super,
                     std::vector<uint8_t>& block_buf, std::vector<uint16_t>& overlap_buf, std::vector<CovRun>& runs_buf) {
    const uint32_t sid = (uint32_t)snodes.size();
    for (uint32_t nid : chain) node_to_super[nid] = sid;
    splice_chain(chain, seqlen, payload, seq_off, block_buf, overlap_buf);
    compute_runs(chain, seqlen, ncov, tcov, overlap_buf, block_buf.size(), runs_buf);
    snodes.push_back({(uint32_t)block_buf.size(), (uint32_t)runs_buf.size(), color[chain[0]]});
    blocks.insert(blocks.end(), block_buf.begin(), block_buf.end());
    runs.insert(runs.end(), runs_buf.begin(), runs_buf.end());
}

}  // namespace

// ---- unitig compaction over the directed colored graph, IN MEMORY (phase 1 of `simplify`) ----
// Builds the compacted super-graph (unitig-compacted, per-base coverage) directly into the
// caller's vectors — never written to disk. run_simplify consumes these straight away.
void build_super_graph(const std::string& dir, std::vector<SuperNode>& snodes,
                       std::vector<uint8_t>& blocks, std::vector<CovRun>& runs, std::vector<IEdge>& super_edges,
                       std::vector<uint32_t>& node_to_super) {
    FILE* fg = xopen(dir + "/graph.bin", "rb");
    uint64_t K; rd(fg, &K, 8);
    std::vector<uint8_t> color(K);                rd(fg, color.data(), K);
    std::vector<uint32_t> ncov(K), tcov(K);       rd(fg, ncov.data(), K * 4); rd(fg, tcov.data(), K * 4);
    uint64_t E; rd(fg, &E, 8);
    std::vector<IEdge> edges(E);                  rd(fg, edges.data(), E * sizeof(IEdge));
    // self-contained graph.bin (2-pass `graph`): the per-node sequence payload follows the
    // edges inline (total, seqlen[K], payload) -- no separate counts.bin round-trip.
    uint64_t total_payload; rd(fg, &total_payload, 8);
    std::vector<uint16_t> seqlen(K);              rd(fg, seqlen.data(), K * 2);
    std::vector<uint8_t> payload(total_payload);  rd(fg, payload.data(), total_payload);
    std::fclose(fg);
    std::vector<uint64_t> seq_off(K + 1, 0);
    for (uint64_t i = 0; i < K; i++) seq_off[i + 1] = seq_off[i] + seqlen[i];
    std::vector<uint32_t> in_deg(K, 0), out_deg(K, 0), succ(K, UINT32_MAX), pred(K, UINT32_MAX);
    for (const IEdge& e : edges) { out_deg[e.first]++; in_deg[e.second]++; succ[e.first] = e.second; pred[e.second] = e.first; }
    std::vector<char> visited(K, 0);
    snodes.clear(); snodes.reserve(K / 2);
    node_to_super.assign(K, UINT32_MAX);
    std::vector<uint32_t> chain;
    std::vector<uint16_t> overlap_buf;
    std::vector<uint8_t> block_buf;
    std::vector<CovRun> runs_buf;
    for (uint32_t u = 0; u < (uint32_t)K; u++) {
        if (visited[u]) continue;
        const bool left_end = (in_deg[u] != 1) || (out_deg[pred[u]] != 1) || (color[pred[u]] != color[u]);
        if (!left_end) continue;
        walk_chain(u, succ, in_deg, out_deg, color, visited, chain);
        emit_super_node(chain, seqlen, payload.data(), seq_off, ncov, tcov, color,
                        snodes, blocks, runs, node_to_super,
                        block_buf, overlap_buf, runs_buf);
    }
    uint64_t cycle_singletons = 0;
    for (uint32_t u = 0; u < (uint32_t)K; u++) {
        if (visited[u]) continue;
        chain.assign(1, u); visited[u] = 1; cycle_singletons++;
        emit_super_node(chain, seqlen, payload.data(), seq_off, ncov, tcov, color,
                        snodes, blocks, runs, node_to_super,
                        block_buf, overlap_buf, runs_buf);
    }
    build_super_edges(edges, node_to_super, super_edges);
    std::fprintf(stderr, "[compact] K=%lu -> S=%lu (cycle-singletons=%lu) edges=%lu->%lu block=%.2fGB runs=%lu\n",
                 (unsigned long)K, (unsigned long)snodes.size(), (unsigned long)cycle_singletons,
                 (unsigned long)E, (unsigned long)super_edges.size(), blocks.size() / 1e9, (unsigned long)runs.size());
}


namespace {

struct GEdge { uint32_t from, to, overlap; };

constexpr uint32_t NODE_BUDGET = 500;   // per-branch DFS budget (mirrors the audit)

struct PairInfo { uint32_t c1_branches = 0, c2_branches = 0; bool direct = false, cycle = false; };
struct BubblePair { uint32_t s1, s2; uint8_t bucket; };   // 0=std_both, 1=T_only_no_N, 2=N_only_no_T

// CSR adjacency: neighbours of u are nbr[off[u] .. off[u+1]). One flat allocation instead
// of S per-node vectors; counting-sort fill preserves edge-insertion order per node.
struct CSR { std::vector<uint32_t> off, nbr; };

void build_csr(const std::vector<IEdge>& super_edges, uint64_t S, bool fwd, CSR& g) {
    g.off.assign(S + 1, 0);
    for (const IEdge& e : super_edges) g.off[(fwd ? e.first : e.second) + 1]++;
    for (uint64_t i = 0; i < S; i++) g.off[i + 1] += g.off[i];
    g.nbr.resize(super_edges.size());
    std::vector<uint32_t> w(g.off.begin(), g.off.begin() + S);   // per-node write cursor
    for (const IEdge& e : super_edges) {
        const uint32_t u = fwd ? e.first : e.second, v = fwd ? e.second : e.first;
        g.nbr[w[u]++] = v;
    }
}

struct OffsetTable { std::vector<uint64_t> blk_off, runs_off; };

void build_offset_table(const std::vector<SuperNode>& snodes, OffsetTable& t) {
    const uint64_t S = snodes.size();
    t.blk_off.assign(S + 1, 0); t.runs_off.assign(S + 1, 0);
    for (uint64_t i = 0; i < S; i++) {
        t.blk_off[i + 1]  = t.blk_off[i]  + snodes[i].block_len;
        t.runs_off[i + 1] = t.runs_off[i] + snodes[i].runs_len;
    }
}

// BFS from s1 to s2; interior nodes restricted by a colour bitmask, start out-neighbour by
// start_mask. Deterministic (layer order, lowest-ID parent tiebreak). path = [s1, ..., s2].
// seen/parent/dist are caller-owned scratch reused across calls; a node belongs to THIS
// call iff seen[v]==gen, so there is no per-call O(S) allocation or reset.
bool bfs_path(uint32_t s1, uint32_t s2, const CSR& og, const std::vector<uint8_t>& color,
              uint8_t start_mask, uint8_t interior_mask, uint32_t max_hops,
              std::vector<uint32_t>& seen, std::vector<uint32_t>& parent, std::vector<uint32_t>& dist,
              uint32_t gen, std::vector<uint32_t>& path) {
    path.clear();
    seen[s1] = gen; parent[s1] = s1; dist[s1] = 0;
    std::queue<uint32_t> q; q.push(s1);
    while (!q.empty()) {
        const uint32_t u = q.front(); q.pop();
        if (dist[u] >= max_hops) continue;
        const uint8_t mask = (u == s1) ? start_mask : interior_mask;
        for (uint32_t i = og.off[u]; i < og.off[u + 1]; i++) {
            const uint32_t v = og.nbr[i];
            if (v == s2 && color[v] == 3) {
                if (seen[v] != gen || u < parent[v]) { seen[v] = gen; parent[v] = u; dist[v] = dist[u] + 1; }
                continue;
            }
            if (!(mask & (1u << (color[v] - 1)))) continue;
            if (seen[v] == gen) continue;
            seen[v] = gen; parent[v] = u; dist[v] = dist[u] + 1; q.push(v);
        }
    }
    if (seen[s2] != gen) return false;
    for (uint32_t cur = s2; cur != s1; cur = parent[cur]) { path.push_back(cur); if (parent[cur] == cur) return false; }
    path.push_back(s1); std::reverse(path.begin(), path.end()); return true;
}

// Splice a path of super-nodes into one byte sequence via KMP suffix-prefix overlap.
void splice_arm(const std::vector<uint32_t>& path, const std::vector<SuperNode>& snodes,
                const std::vector<uint8_t>& blocks, const OffsetTable& ot,
                std::vector<uint8_t>& arm, std::vector<uint64_t>& seg_start) {
    arm.clear(); seg_start.assign(path.size(), 0);
    if (path.empty()) return;
    const uint64_t b0 = ot.blk_off[path[0]]; const uint32_t l0 = snodes[path[0]].block_len;
    arm.insert(arm.end(), blocks.begin() + b0, blocks.begin() + b0 + l0);
    for (size_t i = 1; i < path.size(); i++) {
        const uint32_t li = snodes[path[i]].block_len; const uint8_t* B = blocks.data() + ot.blk_off[path[i]];
        const uint32_t lp = snodes[path[i - 1]].block_len; const uint8_t* A = arm.data() + (arm.size() - lp);
        const uint32_t k = suffix_prefix(A, lp, B, li);
        if (k == 0) { std::fprintf(stderr, "[FATAL] simplify: zero-overlap arm splice path[%zu]->path[%zu]\n", i - 1, i); std::abort(); }
        seg_start[i] = seg_start[i - 1] + lp - k;
        arm.insert(arm.end(), B + k, B + li);
    }
}

// Expand a path's per-base RLE coverage into arm-coordinate cov_n/cov_t, MAX over overlaps.
void expand_arm_coverage(const std::vector<uint32_t>& path, const std::vector<SuperNode>& snodes,
                         const std::vector<CovRun>& runs, const OffsetTable& ot,
                         const std::vector<uint64_t>& seg_start, uint64_t arm_len,
                         std::vector<uint32_t>& cov_n, std::vector<uint32_t>& cov_t) {
    cov_n.assign(arm_len, 0); cov_t.assign(arm_len, 0);
    for (size_t i = 0; i < path.size(); i++) {
        const uint64_t origin = seg_start[i]; const uint32_t bl = snodes[path[i]].block_len;
        const uint64_t r0 = ot.runs_off[path[i]]; const uint32_t rn = snodes[path[i]].runs_len;
        uint32_t prev_end = 0;
        for (uint32_t j = 0; j < rn; j++) {
            const CovRun& r = runs[r0 + j];
            for (uint32_t p = prev_end; p < r.end_pos && p < bl; p++) {
                const uint64_t a = origin + p;
                if (r.n_cov > cov_n[a]) cov_n[a] = r.n_cov;
                if (r.t_cov > cov_t[a]) cov_t[a] = r.t_cov;
            }
            prev_end = r.end_pos;
        }
    }
}

// Expand ONE super-node's RLE into flat per-base coverage.
void expand_node_cov(uint32_t nid, const std::vector<SuperNode>& snodes,
                     const std::vector<CovRun>& runs, const OffsetTable& ot,
                     std::vector<uint32_t>& cn, std::vector<uint32_t>& ct) {
    const uint32_t bl = snodes[nid].block_len;
    cn.assign(bl, 0); ct.assign(bl, 0);
    const uint64_t r0 = ot.runs_off[nid]; const uint32_t rn = snodes[nid].runs_len;
    uint32_t prev = 0;
    for (uint32_t j = 0; j < rn; j++) {
        const CovRun& r = runs[r0 + j];
        for (uint32_t p = prev; p < r.end_pos && p < bl; p++) { cn[p] = r.n_cov; ct[p] = r.t_cov; }
        prev = r.end_pos;
    }
}

// Expand an arm's RLE coverage into flat cn/ct over positions [lo, hi). Arms keep coverage
// RLE (not flat uint32 per base) — a long multi-node bubble path would otherwise carry
// 8 bytes/base resident for the whole stage; only these contiguous slices are ever read.
void cov_slice(const std::vector<CovRun>& cov, uint32_t lo, uint32_t hi,
               std::vector<uint32_t>& cn, std::vector<uint32_t>& ct) {
    cn.assign(hi - lo, 0); ct.assign(hi - lo, 0);
    uint32_t prev = 0;
    for (const CovRun& r : cov) {
        const uint32_t a = std::max(prev, lo), b = std::min(r.end_pos, hi);
        for (uint32_t p = a; p < b; p++) { cn[p - lo] = r.n_cov; ct[p - lo] = r.t_cov; }
        prev = r.end_pos;
        if (prev >= hi) break;
    }
}

// Compress flat cov_n/cov_t into CovRun RLE.
void rle_emit(const std::vector<uint32_t>& cov_n, const std::vector<uint32_t>& cov_t, std::vector<CovRun>& out) {
    out.clear(); const uint64_t L = cov_n.size(); if (L == 0) return;
    uint32_t cur_n = cov_n[0], cur_t = cov_t[0];
    for (uint64_t p = 1; p < L; p++) {
        if (cov_n[p] != cur_n || cov_t[p] != cur_t) {
            out.push_back({(uint32_t)p, cur_n, cur_t}); cur_n = cov_n[p]; cur_t = cov_t[p];
        }
    }
    out.push_back({(uint32_t)L, cur_n, cur_t});
}

// DFS from branch start x, color != 3 interior, record every c3 closure node reached.
// vis/stack are caller-owned scratch; a node is visited in THIS call iff vis[n]==gen.
void dfs_branch(uint32_t s1, uint32_t x, const std::vector<uint8_t>& color, const CSR& og,
                std::vector<uint32_t>& vis, uint32_t gen, std::vector<uint32_t>& stack,
                std::vector<uint32_t>& closed, bool& cycle_hit) {
    closed.clear(); cycle_hit = false;
    stack.clear(); stack.push_back(x); vis[x] = gen;
    uint32_t nodes = 0;
    while (!stack.empty() && nodes < NODE_BUDGET) {
        const uint32_t n = stack.back(); stack.pop_back(); nodes++;
        for (uint32_t i = og.off[n]; i < og.off[n + 1]; i++) {
            const uint32_t nx = og.nbr[i];
            if (nx == s1) { cycle_hit = true; continue; }
            if (color[nx] == 3) { closed.push_back(nx); continue; }
            if (vis[nx] == gen) { cycle_hit = true; continue; }
            vis[nx] = gen; stack.push_back(nx);
        }
    }
}

// Iterate every c3 S1, enumerate (S1, S2) closures by branch colour, classify into buckets.
void enumerate_bubbles(const std::vector<uint8_t>& color, const CSR& og,
                       std::vector<BubblePair>& out, uint64_t& direct_only, uint64_t& cycle_skipped, uint64_t& multi_path) {
    out.clear(); direct_only = cycle_skipped = multi_path = 0;
    const uint64_t S = color.size();
    std::vector<uint32_t> vis(S, 0), stack, closed; uint32_t gen = 0;
    std::unordered_map<uint64_t, PairInfo> pair_info;
    for (uint32_t s1 = 0; s1 < (uint32_t)S; s1++) {
        if (color[s1] != 3) continue;
        pair_info.clear();
        for (uint32_t i = og.off[s1]; i < og.off[s1 + 1]; i++) { const uint32_t x = og.nbr[i]; if (color[x] == 3) pair_info[((uint64_t)s1 << 32) | x].direct = true; }
        for (uint32_t i = og.off[s1]; i < og.off[s1 + 1]; i++) {
            const uint32_t x = og.nbr[i];
            if (color[x] == 3) continue;
            const uint8_t branch_c = color[x]; bool cy = false;
            dfs_branch(s1, x, color, og, vis, ++gen, stack, closed, cy);
            for (uint32_t s2 : closed) {
                PairInfo& p = pair_info[((uint64_t)s1 << 32) | s2];
                if (branch_c == 1) p.c1_branches++; else p.c2_branches++;
                if (cy) p.cycle = true;
            }
        }
        for (const auto& kv : pair_info) {
            const uint32_t s2 = (uint32_t)kv.first; const PairInfo& v = kv.second;
            const bool has1 = v.c1_branches > 0, has2 = v.c2_branches > 0;
            if (!has1 && !has2) { direct_only++; continue; }
            if (v.cycle) { cycle_skipped++; continue; }
            if (v.c1_branches > 1 || v.c2_branches > 1) { multi_path++; continue; }   // v1: leave tangles intact
            if (has1 && has2)        out.push_back({s1, s2, 0});
            else if (has2 && !has1)  out.push_back({s1, s2, 1});
            else                     out.push_back({s1, s2, 2});
        }
    }
}

// ---- output graph accumulator ----
struct OutGraph { std::vector<SuperNode> snodes; std::vector<uint8_t> blocks; std::vector<CovRun> runs; std::vector<GEdge> edges; };

uint32_t out_add(OutGraph& g, const std::vector<uint8_t>& bytes,
                 const std::vector<uint32_t>& cn, const std::vector<uint32_t>& ct, uint8_t color) {
    const uint32_t id = (uint32_t)g.snodes.size();
    std::vector<CovRun> r; rle_emit(cn, ct, r);
    g.snodes.push_back({(uint32_t)bytes.size(), (uint32_t)r.size(), color});
    g.blocks.insert(g.blocks.end(), bytes.begin(), bytes.end());
    g.runs.insert(g.runs.end(), r.begin(), r.end());
    return id;
}

struct Arm { uint32_t s1, s2; std::vector<uint8_t> bytes; std::vector<CovRun> cov; uint8_t color; };
struct ZBubble { int arm_n, arm_t; };

// Longest common prefix / suffix length over a set of allele byte vectors.
uint32_t lcp_over(const std::vector<Arm>& arms, const std::vector<int>& idx) {
    if (idx.empty()) return 0;
    uint32_t L = UINT32_MAX; for (int i : idx) L = std::min(L, (uint32_t)arms[i].bytes.size());
    uint32_t p = 0;
    while (p < L) {
        const uint8_t c = arms[idx[0]].bytes[p]; bool all = true;
        for (int i : idx) if (arms[i].bytes[p] != c) { all = false; break; }
        if (!all) break;
        p++;
    }
    return p;
}
uint32_t lcs_over(const std::vector<Arm>& arms, const std::vector<int>& idx) {
    if (idx.empty()) return 0;
    uint32_t L = UINT32_MAX; for (int i : idx) L = std::min(L, (uint32_t)arms[i].bytes.size());
    uint32_t s = 0;
    while (s < L) {
        const Arm& a0 = arms[idx[0]]; const uint8_t c = a0.bytes[a0.bytes.size() - 1 - s]; bool all = true;
        for (int i : idx) { const Arm& a = arms[i]; if (a.bytes[a.bytes.size() - 1 - s] != c) { all = false; break; } }
        if (!all) break;
        s++;
    }
    return s;
}

void write_simplified_graph(const std::string& path, const OutGraph& g) {
    const std::string tmp = path + ".tmp"; FILE* f = xopen(tmp, "wb");
    uint64_t S = g.snodes.size();  wr(f, &S, 8);  wr(f, g.snodes.data(), S * sizeof(SuperNode));
    uint64_t tb = g.blocks.size(); wr(f, &tb, 8); wr(f, g.blocks.data(), tb);
    uint64_t tr = g.runs.size();   wr(f, &tr, 8); wr(f, g.runs.data(),   tr * sizeof(CovRun));
    uint64_t E = g.edges.size();   wr(f, &E,  8); wr(f, g.edges.data(),  E * sizeof(GEdge));
    atomic_close(f, tmp, path);
}

}  // namespace

void run_simplify(const std::string& dir) {
    std::vector<SuperNode> snodes;
    std::vector<uint8_t> blocks; std::vector<CovRun> runs; std::vector<IEdge> super_edges;
    std::vector<uint32_t> node_to_super;
    build_super_graph(dir, snodes, blocks, runs, super_edges, node_to_super);    // unitig compaction, straight into RAM
    resolve_repeats(dir, snodes, blocks, runs, super_edges, node_to_super);      // X-cut repeat resolution (no-op if no threads.bin)
    const uint32_t S = (uint32_t)snodes.size();
    std::vector<uint8_t> color(S); for (uint32_t i = 0; i < S; i++) color[i] = snodes[i].color;
    CSR out_adj, in_adj; build_csr(super_edges, S, true, out_adj); build_csr(super_edges, S, false, in_adj);
    OffsetTable ot; build_offset_table(snodes, ot);

    std::vector<BubblePair> bubbles; uint64_t direct_only = 0, cycle_skipped = 0, multi_path = 0;
    enumerate_bubbles(color, out_adj, bubbles, direct_only, cycle_skipped, multi_path);

    // ---- Phase 1: reconstruct std_both arms; mark consumed interiors ----
    constexpr uint8_t M1 = 0x1, M2 = 0x2, M12 = M1 | M2;
    std::vector<Arm> arms; std::vector<ZBubble> zb;
    std::vector<char> consumed(S, 0), is_anchor(S, 0);
    std::vector<uint32_t> bseen(S, 0), bpar(S, 0), bdist(S, 0); uint32_t bgen = 0;   // bfs scratch
    std::vector<uint32_t> fcn, fct;                                                  // reused flat cov buffers
    uint64_t zipped = 0, passed = 0;
    for (const BubblePair& bp : bubbles) {
        if (bp.bucket != 0) { passed++; continue; }            // v1: only std_both is zipped
        std::vector<uint32_t> pn, pt;
        if (!bfs_path(bp.s1, bp.s2, out_adj, color, M1, M12, NODE_BUDGET, bseen, bpar, bdist, ++bgen, pn) ||
            !bfs_path(bp.s1, bp.s2, out_adj, color, M2, M12, NODE_BUDGET, bseen, bpar, bdist, ++bgen, pt)) { passed++; continue; }
        bool conflict = false;                                  // don't double-consume a shared interior
        for (size_t i = 1; i + 1 < pn.size() && !conflict; i++) conflict = consumed[pn[i]];
        for (size_t i = 1; i + 1 < pt.size() && !conflict; i++) conflict = consumed[pt[i]];
        if (conflict) { passed++; continue; }
        for (size_t i = 1; i + 1 < pn.size(); i++) consumed[pn[i]] = 1;
        for (size_t i = 1; i + 1 < pt.size(); i++) consumed[pt[i]] = 1;
        is_anchor[bp.s1] = 1; is_anchor[bp.s2] = 1;
        Arm an, at; an.s1 = at.s1 = bp.s1; an.s2 = at.s2 = bp.s2; an.color = 1; at.color = 2;
        std::vector<uint64_t> sgn, sgt;
        splice_arm(pn, snodes, blocks, ot, an.bytes, sgn);
        splice_arm(pt, snodes, blocks, ot, at.bytes, sgt);
        expand_arm_coverage(pn, snodes, runs, ot, sgn, an.bytes.size(), fcn, fct); rle_emit(fcn, fct, an.cov);
        expand_arm_coverage(pt, snodes, runs, ot, sgt, at.bytes.size(), fcn, fct); rle_emit(fcn, fct, at.cov);
        const int in = (int)arms.size(); arms.push_back(std::move(an));
        const int it = (int)arms.size(); arms.push_back(std::move(at));
        zb.push_back({in, it}); zipped++;
    }

    // ---- Phase 2: per-anchor grow amounts (cross-bubble common prefix/suffix), only when "private" ----
    std::unordered_map<uint32_t, std::vector<int>> s1_arms, s2_arms;
    for (int i = 0; i < (int)arms.size(); i++) { s1_arms[arms[i].s1].push_back(i); s2_arms[arms[i].s2].push_back(i); }
    std::unordered_map<uint32_t, uint32_t> grow_pre, grow_suf;
    for (auto& kv : s1_arms) {
        bool priv = true; for (uint32_t i = out_adj.off[kv.first]; i < out_adj.off[kv.first + 1]; i++) if (!consumed[out_adj.nbr[i]]) { priv = false; break; }
        if (priv) grow_pre[kv.first] = lcp_over(arms, kv.second);
    }
    for (auto& kv : s2_arms) {
        bool priv = true; for (uint32_t i = in_adj.off[kv.first]; i < in_adj.off[kv.first + 1]; i++) if (!consumed[in_adj.nbr[i]]) { priv = false; break; }
        if (priv) grow_suf[kv.first] = lcs_over(arms, kv.second);
    }

    // ---- Phase 3a: emit kept nodes (interiors dropped; anchors grown) ----
    OutGraph g; std::vector<uint32_t> new_id(S, UINT32_MAX);
    std::vector<uint32_t> own_cn, own_ct, scn, sct;
    for (uint32_t nid = 0; nid < S; nid++) {
        if (consumed[nid]) continue;
        const uint32_t bl = snodes[nid].block_len;
        const uint8_t* blk = blocks.data() + ot.blk_off[nid];
        expand_node_cov(nid, snodes, runs, ot, own_cn, own_ct);
        if (!is_anchor[nid]) { std::vector<uint8_t> b(blk, blk + bl); new_id[nid] = out_add(g, b, own_cn, own_ct, color[nid]); continue; }
        const uint32_t gp = grow_pre.count(nid) ? grow_pre[nid] : bl;     // >= bl (lcp includes the block)
        const uint32_t gs = grow_suf.count(nid) ? grow_suf[nid] : bl;
        const uint32_t add_pre = gp > bl ? gp - bl : 0;                   // appended toward the downstream bubble
        const uint32_t add_suf = gs > bl ? gs - bl : 0;                   // prepended toward the upstream bubble
        std::vector<uint8_t> b; std::vector<uint32_t> cn, ct;
        b.reserve(add_suf + bl + add_pre); cn.reserve(add_suf + bl + add_pre); ct.reserve(add_suf + bl + add_pre);
        if (add_suf) { const std::vector<int>& idx = s2_arms[nid]; const Arm& a0 = arms[idx[0]]; const uint32_t L0 = (uint32_t)a0.bytes.size();
            for (uint32_t p = 0; p < add_suf; p++) b.push_back(a0.bytes[L0 - gs + p]);
            std::vector<uint32_t> mn(add_suf, 0), mt(add_suf, 0);
            for (int i : idx) { const Arm& a = arms[i]; const uint32_t L = (uint32_t)a.bytes.size();
                cov_slice(a.cov, L - gs, L - gs + add_suf, scn, sct);
                for (uint32_t p = 0; p < add_suf; p++) { mn[p] = std::max(mn[p], scn[p]); mt[p] = std::max(mt[p], sct[p]); } }
            for (uint32_t p = 0; p < add_suf; p++) { cn.push_back(mn[p]); ct.push_back(mt[p]); } }
        for (uint32_t p = 0; p < bl; p++) b.push_back(blk[p]);
        for (uint32_t p = 0; p < bl; p++) { cn.push_back(own_cn[p]); ct.push_back(own_ct[p]); }
        if (add_pre) { const std::vector<int>& idx = s1_arms[nid]; const Arm& a0 = arms[idx[0]];
            for (uint32_t p = 0; p < add_pre; p++) b.push_back(a0.bytes[bl + p]);
            std::vector<uint32_t> mn(add_pre, 0), mt(add_pre, 0);
            for (int i : idx) { const Arm& a = arms[i];
                cov_slice(a.cov, bl, bl + add_pre, scn, sct);
                for (uint32_t p = 0; p < add_pre; p++) { mn[p] = std::max(mn[p], scn[p]); mt[p] = std::max(mt[p], sct[p]); } }
            for (uint32_t p = 0; p < add_pre; p++) { cn.push_back(mn[p]); ct.push_back(mt[p]); } }
        new_id[nid] = out_add(g, b, cn, ct, 3);
    }

    // ---- Phase 3b: residual arm nodes + abutting bubble-seam edges ----
    uint64_t residuals = 0, collapsed = 0;
    for (const ZBubble& z : zb) {
        for (int ai : {z.arm_n, z.arm_t}) {
            const Arm& a = arms[ai];
            const uint32_t gp = grow_pre.count(a.s1) ? grow_pre[a.s1] : snodes[a.s1].block_len;
            const uint32_t gs = grow_suf.count(a.s2) ? grow_suf[a.s2] : snodes[a.s2].block_len;
            const uint32_t lo = gp, hi = (uint32_t)a.bytes.size() - gs;
            const uint32_t n1 = new_id[a.s1], n2 = new_id[a.s2];
            if (hi <= lo) { g.edges.push_back({n1, n2, 0}); collapsed++; continue; }   // alleles identical -> direct
            std::vector<uint8_t> rb(a.bytes.begin() + lo, a.bytes.begin() + hi);
            std::vector<uint32_t> rcn, rct; cov_slice(a.cov, lo, hi, rcn, rct);
            const uint32_t rid = out_add(g, rb, rcn, rct, a.color);
            g.edges.push_back({n1, rid, 0}); g.edges.push_back({rid, n2, 0}); residuals++;
        }
    }

    // ---- Phase 3c: pass-through edges (neither endpoint consumed); recompute overlap on output bytes ----
    OffsetTable got; build_offset_table(g.snodes, got);
    for (const IEdge& e : super_edges) {
        if (consumed[e.first] || consumed[e.second]) continue;
        const uint32_t a = new_id[e.first], b = new_id[e.second];
        if (a == UINT32_MAX || b == UINT32_MAX) continue;
        const uint32_t la = g.snodes[a].block_len, lb = g.snodes[b].block_len;
        const uint32_t k = suffix_prefix(g.blocks.data() + got.blk_off[a], la, g.blocks.data() + got.blk_off[b], lb);
        g.edges.push_back({a, b, k});
    }
    std::sort(g.edges.begin(), g.edges.end(), [](const GEdge& x, const GEdge& y){ return x.from != y.from ? x.from < y.from : (x.to != y.to ? x.to < y.to : x.overlap < y.overlap); });
    g.edges.erase(std::unique(g.edges.begin(), g.edges.end(), [](const GEdge& x, const GEdge& y){ return x.from == y.from && x.to == y.to && x.overlap == y.overlap; }), g.edges.end());

    write_simplified_graph(dir + "/simplify.bin", g);
    std::fprintf(stderr, "[simplify] in: S=%u E=%lu | bubbles=%lu (direct=%lu cycle=%lu multi=%lu) zipped=%lu passed=%lu | grew_pre=%lu grew_suf=%lu residual_arms=%lu collapsed=%lu | out: S=%lu E=%lu blocks=%.2fMB -> %s/simplify.bin\n",
                 S, (unsigned long)super_edges.size(), (unsigned long)bubbles.size(),
                 (unsigned long)direct_only, (unsigned long)cycle_skipped, (unsigned long)multi_path,
                 (unsigned long)zipped, (unsigned long)passed,
                 (unsigned long)grow_pre.size(), (unsigned long)grow_suf.size(), (unsigned long)residuals, (unsigned long)collapsed,
                 (unsigned long)g.snodes.size(), (unsigned long)g.edges.size(), g.blocks.size() / 1e6, dir.c_str());
}
