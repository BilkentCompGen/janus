#include "graph.h"
#include "svp.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
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
                       std::vector<uint8_t>& blocks, std::vector<CovRun>& runs, std::vector<IEdge>& super_edges) {
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
    std::vector<uint32_t> node_to_super(K, UINT32_MAX), chain;
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
