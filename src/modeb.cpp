// modeb: per-bubble N-vs-T mutual-trim refinement on top of Mode A. For each (S1, S2)
// bubble — both color=3 anchors with a parallel allele path pair between them — find
// the shared prefix and shared suffix of the N-allele and T-allele byte sequences,
// attribute them logically to S1's tail / S2's head, and emit the residuals tagged
// EQ / D / I / S / C. modea.bin is NOT mutated; modeb.bin holds the per-bubble records.
//
// Three processable bucket classes mirror tools/modea_audit.py:
//   std_both       : parallel color=1 and color=2 paths between S1, S2
//   T_only_no_N    : color=2 path + color=3 backbone detour (≤8 hops) as the N-allele
//   N_only_no_T    : color=1 path + color=3 backbone detour (≤8 hops) as the T-allele
// direct_only / *_cycle / multi_path_both are skipped or absorbed into C.
//
// Initial stage: enumeration only — bucket counts must agree with the Python audit
// before splice/trim is added on top (Plan agent's first-commit gate).
#include "graph.h"
#include "svp.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <queue>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>

namespace {

struct SuperNode { uint32_t block_len, runs_len; uint8_t color; };
struct CovRun    { uint32_t end_pos, n_cov, t_cov; };

// Mirror tools/modea_audit.py:44 — per-branch DFS budget. Tracks compatibility.
constexpr uint32_t NODE_BUDGET = 500;
// Mirror tools/modea_audit_orphans.py:74 — c3-only detour BFS hop cap.
constexpr uint32_t MAX_DETOUR_HOPS = 8;

struct PairInfo { uint32_t c1_branches = 0, c2_branches = 0; bool direct = false, cycle = false; };
struct BubblePair { uint32_t s1, s2; uint8_t bucket; };   // 0=std_both, 1=T_only_no_N, 2=N_only_no_T

struct BubbleTag {
    uint32_t s1, s2;
    uint8_t  tag, bucket; uint16_t _pad;     // tag: 'E' | 'D' | 'I' | 'S' | 'C'
    uint32_t pre, suf;
    uint32_t n_resid_off, n_resid_len, t_resid_off, t_resid_len;
    uint32_t n_runs_off, n_runs_len, t_runs_off, t_runs_len;
};

uint64_t expected_modea_size(uint64_t S, uint64_t tb, uint64_t tr, uint64_t E) {
    return 8 + S * sizeof(SuperNode) + 8 + tb + 8 + tr * sizeof(CovRun) + 8 + E * sizeof(IEdge);
}

void read_modea(const std::string& dir,
                std::vector<SuperNode>& snodes, std::vector<uint8_t>& blocks,
                std::vector<CovRun>& runs, std::vector<IEdge>& super_edges) {
    const std::string path = dir + "/modea.bin";
    struct stat st; if (::stat(path.c_str(), &st)) { std::fprintf(stderr, "[FATAL] stat %s\n", path.c_str()); std::abort(); }
    FILE* f = xopen(path, "rb");
    uint64_t S; rd(f, &S, 8); snodes.resize(S); rd(f, snodes.data(), S * sizeof(SuperNode));
    uint64_t tb; rd(f, &tb, 8); blocks.resize(tb); rd(f, blocks.data(), tb);
    uint64_t tr; rd(f, &tr, 8); runs.resize(tr); rd(f, runs.data(), tr * sizeof(CovRun));
    uint64_t E;  rd(f, &E,  8); super_edges.resize(E); rd(f, super_edges.data(), E * sizeof(IEdge));
    std::fclose(f);
    const uint64_t expect = expected_modea_size(S, tb, tr, E);
    if ((uint64_t)st.st_size != expect) { std::fprintf(stderr, "[FATAL] modea.bin size %lu != expected %lu\n", (unsigned long)st.st_size, (unsigned long)expect); std::abort(); }
}

void build_out_adj(const std::vector<IEdge>& super_edges, uint64_t S, std::vector<std::vector<uint32_t>>& out_adj) {
    out_adj.assign(S, {});
    for (const IEdge& e : super_edges) out_adj[e.first].push_back(e.second);
}

struct OffsetTable { std::vector<uint64_t> blk_off, runs_off; };

void build_offset_table(const std::vector<SuperNode>& snodes, OffsetTable& t) {
    const uint64_t S = snodes.size();
    t.blk_off.assign(S + 1, 0); t.runs_off.assign(S + 1, 0);
    for (uint64_t i = 0; i < S; i++) {
        t.blk_off[i + 1]   = t.blk_off[i]   + snodes[i].block_len;
        t.runs_off[i + 1]  = t.runs_off[i]  + snodes[i].runs_len;
    }
}

// KMP longest exact suffix-of-A == prefix-of-B. Verbatim from src/modea.cpp:16-31.
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

// BFS from s1 to s2 with interior nodes restricted by a colour bitmask (bit 0 = c1, bit 1 = c2,
// bit 2 = c3) and starting out-neighbour restricted by start_mask. Both s1 and s2 must be c3.
// Deterministic: in BFS layer order with lowest-ID parent tiebreak. Returns false if no path.
bool bfs_path(uint32_t s1, uint32_t s2, const std::vector<std::vector<uint32_t>>& out_adj,
              const std::vector<uint8_t>& color, uint8_t start_mask, uint8_t interior_mask, uint32_t max_hops,
              std::vector<uint32_t>& path) {
    path.clear();
    const uint64_t S = out_adj.size();
    std::vector<uint32_t> parent(S, UINT32_MAX), dist(S, UINT32_MAX);
    parent[s1] = s1; dist[s1] = 0;
    std::queue<uint32_t> q; q.push(s1);
    while (!q.empty()) {
        const uint32_t u = q.front(); q.pop();
        if (dist[u] >= max_hops) continue;
        const uint8_t mask = (u == s1) ? start_mask : interior_mask;
        for (uint32_t v : out_adj[u]) {
            if (v == s2 && color[v] == 3) {
                if (parent[v] == UINT32_MAX || u < parent[v]) { parent[v] = u; dist[v] = dist[u] + 1; }
                continue;
            }
            if (!(mask & (1u << (color[v] - 1)))) continue;
            if (parent[v] != UINT32_MAX) continue;
            parent[v] = u; dist[v] = dist[u] + 1; q.push(v);
        }
    }
    if (parent[s2] == UINT32_MAX) return false;
    for (uint32_t cur = s2; cur != s1; cur = parent[cur]) { path.push_back(cur); if (parent[cur] == cur) return false; }
    path.push_back(s1); std::reverse(path.begin(), path.end()); return true;
}

// Splice a path of super-nodes into a single arm byte sequence via KMP suffix-prefix overlap.
// k_arm[i] = exact overlap between path[i] and path[i+1] (in [1, min(block_len[i], block_len[i+1])]).
// k_arm == 0 is FATAL — adjacent LCP-derived super-nodes always share ≥1 overlap byte.
void splice_arm(const std::vector<uint32_t>& path, const std::vector<SuperNode>& snodes,
                const std::vector<uint8_t>& blocks, const OffsetTable& ot,
                std::vector<uint8_t>& arm, std::vector<uint32_t>& k_arm, std::vector<uint64_t>& seg_start) {
    arm.clear(); k_arm.clear(); seg_start.assign(path.size(), 0);
    if (path.empty()) return;
    const uint64_t b0 = ot.blk_off[path[0]]; const uint32_t l0 = snodes[path[0]].block_len;
    arm.insert(arm.end(), blocks.begin() + b0, blocks.begin() + b0 + l0);
    for (size_t i = 1; i < path.size(); i++) {
        const uint32_t li = snodes[path[i]].block_len; const uint8_t* B = blocks.data() + ot.blk_off[path[i]];
        const uint32_t lp = snodes[path[i - 1]].block_len; const uint8_t* A = arm.data() + (arm.size() - lp);
        const uint32_t k = suffix_prefix(A, lp, B, li);
        if (k == 0) { std::fprintf(stderr, "[FATAL] modeb: zero-overlap arm splice path[%zu]->path[%zu]\n", i - 1, i); std::abort(); }
        k_arm.push_back(k); seg_start[i] = seg_start[i - 1] + lp - k;
        arm.insert(arm.end(), B + k, B + li);
    }
}

// Expand each super-node's per-base RLE into the arm-coordinate cov_n/cov_t arrays,
// taking MAX over overlapping windows. seg_start[i] is the arm-coord origin of path[i].
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

// Compress flat cov_n/cov_t arrays into the same CovRun RLE format Mode A uses.
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

// LCP(N, T) and reverse-LCP(N, T) with clamp pre + suf ≤ min(|N|, |T|), prefix priority.
void mutual_trim(const std::vector<uint8_t>& N, const std::vector<uint8_t>& T, uint32_t& pre, uint32_t& suf) {
    const uint32_t lim = (uint32_t)std::min(N.size(), T.size());
    pre = 0; while (pre < lim && N[pre] == T[pre]) pre++;
    suf = 0; while (suf < lim - pre && N[N.size() - 1 - suf] == T[T.size() - 1 - suf]) suf++;
}

// Slice an RLE coverage stream into the arm-coord range [lo, hi). Output positions are
// re-based so out's last CovRun.end_pos == (hi - lo).
void slice_runs(const std::vector<CovRun>& runs, uint64_t lo, uint64_t hi, std::vector<CovRun>& out) {
    out.clear(); if (lo >= hi) return; uint32_t prev_end = 0;
    for (const CovRun& r : runs) {
        const uint32_t r_lo = prev_end, r_hi = r.end_pos; prev_end = r.end_pos;
        if (r_hi <= lo) continue;
        if (r_lo >= hi) break;
        const uint32_t clipped_hi = (uint32_t)std::min<uint64_t>(r_hi, hi);
        out.push_back({(uint32_t)(clipped_hi - lo), r.n_cov, r.t_cov});
    }
}

// EQ/D/I/S/C dispatch. multi_path_both (collapsed into bucket 0 with c1>1 or c2>1) and any
// pair where the structure resisted clean classification falls into C.
char classify_tag(uint32_t n_resid_len, uint32_t t_resid_len) {
    if (n_resid_len == 0 && t_resid_len == 0) return 'E';
    if (t_resid_len == 0)                     return 'D';
    if (n_resid_len == 0)                     return 'I';
    if (n_resid_len == t_resid_len)           return 'S';
    return 'C';
}

// Reconstruct N-arm and T-arm super-node paths for one bubble, splice their bytes, expand
// per-base coverage with MAX-merge, mutual-trim, emit residuals + RLE. Path picker uses
// bfs_path: lowest-ID at parent set (deterministic).
bool process_bubble(const BubblePair& bp, const std::vector<SuperNode>& snodes,
                    const std::vector<uint8_t>& blocks, const std::vector<CovRun>& runs,
                    const OffsetTable& ot, const std::vector<uint8_t>& color,
                    const std::vector<std::vector<uint32_t>>& out_adj,
                    BubbleTag& tag_out, std::vector<uint8_t>& bytes_out, std::vector<CovRun>& runs_out) {
    std::vector<uint32_t> path_n, path_t;
    // Colour masks per bucket — start_mask = which colours at the s1 out-neighbour are valid
    // arm starts; interior_mask = the audit's permissive non-c3 walk (c1|c2) or c3-only detour.
    constexpr uint8_t M1 = 0x1, M2 = 0x2, M3 = 0x4, M12 = M1 | M2;
    uint8_t sn_mask = 0, in_n_mask = 0, st_mask = 0, in_t_mask = 0;
    uint32_t hops_n = NODE_BUDGET, hops_t = NODE_BUDGET;
    if (bp.bucket == 0)      { sn_mask = M1; in_n_mask = M12; st_mask = M2; in_t_mask = M12; }
    else if (bp.bucket == 1) { sn_mask = M3; in_n_mask = M3;  hops_n = MAX_DETOUR_HOPS; st_mask = M2; in_t_mask = M12; }
    else                     { sn_mask = M1; in_n_mask = M12; st_mask = M3; in_t_mask = M3; hops_t = MAX_DETOUR_HOPS; }
    // For orphan buckets (1, 2) the c3 detour may not exist within MAX_DETOUR_HOPS — treat
    // the missing arm as empty bytes so the residual collapses into the surviving arm,
    // tagging as I (T_only_no_N without detour) or D (N_only_no_T without detour).
    const bool ok_n = bfs_path(bp.s1, bp.s2, out_adj, color, sn_mask, in_n_mask, hops_n, path_n);
    const bool ok_t = bfs_path(bp.s1, bp.s2, out_adj, color, st_mask, in_t_mask, hops_t, path_t);
    if (bp.bucket == 0 && (!ok_n || !ok_t)) return false;
    if (bp.bucket == 1 && !ok_t) return false;
    if (bp.bucket == 2 && !ok_n) return false;
    std::vector<uint8_t> N, T; std::vector<uint32_t> kn, kt; std::vector<uint64_t> sn, st;
    if (ok_n) splice_arm(path_n, snodes, blocks, ot, N, kn, sn);
    if (ok_t) splice_arm(path_t, snodes, blocks, ot, T, kt, st);
    std::vector<uint32_t> cn_n, cn_t, ct_n, ct_t;
    if (ok_n) expand_arm_coverage(path_n, snodes, runs, ot, sn, N.size(), cn_n, cn_t);
    if (ok_t) expand_arm_coverage(path_t, snodes, runs, ot, st, T.size(), ct_n, ct_t);
    std::vector<CovRun> rn, rt; rle_emit(cn_n, cn_t, rn); rle_emit(ct_n, ct_t, rt);
    uint32_t pre = 0, suf = 0; mutual_trim(N, T, pre, suf);
    const uint64_t n_lo = pre, n_hi = N.size() - suf, t_lo = pre, t_hi = T.size() - suf;
    std::vector<CovRun> rn_resid, rt_resid; slice_runs(rn, n_lo, n_hi, rn_resid); slice_runs(rt, t_lo, t_hi, rt_resid);
    tag_out.s1 = bp.s1; tag_out.s2 = bp.s2; tag_out.bucket = bp.bucket; tag_out._pad = 0;
    tag_out.pre = pre; tag_out.suf = suf;
    tag_out.n_resid_off = (uint32_t)bytes_out.size(); tag_out.n_resid_len = (uint32_t)(n_hi - n_lo);
    bytes_out.insert(bytes_out.end(), N.begin() + n_lo, N.begin() + n_hi);
    tag_out.t_resid_off = (uint32_t)bytes_out.size(); tag_out.t_resid_len = (uint32_t)(t_hi - t_lo);
    bytes_out.insert(bytes_out.end(), T.begin() + t_lo, T.begin() + t_hi);
    tag_out.n_runs_off = (uint32_t)runs_out.size(); tag_out.n_runs_len = (uint32_t)rn_resid.size();
    runs_out.insert(runs_out.end(), rn_resid.begin(), rn_resid.end());
    tag_out.t_runs_off = (uint32_t)runs_out.size(); tag_out.t_runs_len = (uint32_t)rt_resid.size();
    runs_out.insert(runs_out.end(), rt_resid.begin(), rt_resid.end());
    tag_out.tag = classify_tag(tag_out.n_resid_len, tag_out.t_resid_len);
    return true;
}

void write_modeb_bin(const std::string& path, const std::vector<BubbleTag>& tags,
                     const std::vector<uint8_t>& bytes, const std::vector<CovRun>& runs) {
    const std::string tmp = path + ".tmp"; FILE* f = xopen(tmp, "wb");
    uint64_t B = tags.size();   wr(f, &B,  8); wr(f, tags.data(),  B  * sizeof(BubbleTag));
    uint64_t Lb = bytes.size(); wr(f, &Lb, 8); wr(f, bytes.data(), Lb);
    uint64_t Lr = runs.size();  wr(f, &Lr, 8); wr(f, runs.data(),  Lr * sizeof(CovRun));
    atomic_close(f, tmp, path);
}

// DFS from branch start x, color != 3 interior, record every c3 closure node reached.
// Mirror modea_audit.py:73-100 termination rule; cycle_hit fires on back-edge or s1 self-loop.
void dfs_branch(uint32_t s1, uint32_t x, const std::vector<uint8_t>& color,
                const std::vector<std::vector<uint32_t>>& out_adj,
                std::vector<uint32_t>& closed, bool& cycle_hit, bool& truncated) {
    closed.clear(); cycle_hit = false; truncated = false;
    std::vector<char> visited(out_adj.size(), 0);
    std::vector<uint32_t> stack; stack.push_back(x); visited[x] = 1;
    uint32_t nodes = 0;
    while (!stack.empty() && nodes < NODE_BUDGET) {
        const uint32_t n = stack.back(); stack.pop_back(); nodes++;
        for (uint32_t nx : out_adj[n]) {
            if (nx == s1) { cycle_hit = true; continue; }
            if (color[nx] == 3) { closed.push_back(nx); continue; }
            if (visited[nx]) { cycle_hit = true; continue; }
            visited[nx] = 1; stack.push_back(nx);
        }
    }
    if (nodes >= NODE_BUDGET) truncated = true;
}

// Iterate every color=3 super-node S1, enumerate (S1, S2) closures by branch colour, and
// classify into the three processable buckets. Exact mirror of tools/modea_audit.py:59-133.
void enumerate_bubbles(const std::vector<uint8_t>& color, const std::vector<std::vector<uint32_t>>& out_adj,
                       std::vector<BubblePair>& out, uint64_t& direct_only, uint64_t& cycle_skipped,
                       uint64_t& multi_path_C, uint64_t& truncated_starts) {
    out.clear(); direct_only = cycle_skipped = multi_path_C = truncated_starts = 0;
    const uint64_t S = color.size();
    std::vector<uint32_t> closed; std::unordered_map<uint64_t, PairInfo> pair_info;
    for (uint32_t s1 = 0; s1 < (uint32_t)S; s1++) {
        if (color[s1] != 3) continue;
        pair_info.clear();
        for (uint32_t x : out_adj[s1]) { if (color[x] == 3) pair_info[((uint64_t)s1 << 32) | x].direct = true; }
        for (uint32_t x : out_adj[s1]) {
            if (color[x] == 3) continue;
            const uint8_t branch_c = color[x]; bool cy = false, tr = false;
            dfs_branch(s1, x, color, out_adj, closed, cy, tr);
            if (tr) truncated_starts++;
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
            if (v.c1_branches > 1 || v.c2_branches > 1) { multi_path_C++; out.push_back({s1, s2, 0}); continue; }
            if (has1 && has2)        out.push_back({s1, s2, 0});       // std_both
            else if (has2 && !has1)  out.push_back({s1, s2, 1});       // T_only_no_N
            else                     out.push_back({s1, s2, 2});       // N_only_no_T
        }
    }
}

void summarize_bubbles(const std::vector<BubblePair>& bubbles, uint64_t& n_std, uint64_t& n_t_only, uint64_t& n_n_only) {
    n_std = n_t_only = n_n_only = 0;
    for (const BubblePair& b : bubbles) {
        if      (b.bucket == 0) n_std++;
        else if (b.bucket == 1) n_t_only++;
        else                    n_n_only++;
    }
}

}  // namespace

void run_modeb(const std::string& dir) {
    std::vector<SuperNode> snodes;
    std::vector<uint8_t> blocks; std::vector<CovRun> runs; std::vector<IEdge> super_edges;
    read_modea(dir, snodes, blocks, runs, super_edges);
    std::vector<uint8_t> color(snodes.size()); for (size_t i = 0; i < snodes.size(); i++) color[i] = snodes[i].color;
    std::vector<std::vector<uint32_t>> out_adj; build_out_adj(super_edges, snodes.size(), out_adj);
    OffsetTable ot; build_offset_table(snodes, ot);
    std::vector<BubblePair> bubbles; uint64_t direct_only = 0, cycle_skipped = 0, multi_path_C = 0, truncated_starts = 0;
    enumerate_bubbles(color, out_adj, bubbles, direct_only, cycle_skipped, multi_path_C, truncated_starts);
    uint64_t n_std = 0, n_t_only = 0, n_n_only = 0; summarize_bubbles(bubbles, n_std, n_t_only, n_n_only);
    std::vector<BubbleTag> tags; std::vector<uint8_t> bytes_out; std::vector<CovRun> runs_out;
    tags.reserve(bubbles.size()); BubbleTag bt = {};
    uint64_t no_path_skipped = 0; uint64_t tag_count[5] = {0, 0, 0, 0, 0};   // E D I S C
    for (const BubblePair& bp : bubbles) {
        if (!process_bubble(bp, snodes, blocks, runs, ot, color, out_adj, bt, bytes_out, runs_out)) { no_path_skipped++; continue; }
        tags.push_back(bt);
        switch (bt.tag) { case 'E': tag_count[0]++; break; case 'D': tag_count[1]++; break; case 'I': tag_count[2]++; break; case 'S': tag_count[3]++; break; default: tag_count[4]++; }
    }
    write_modeb_bin(dir + "/modeb.bin", tags, bytes_out, runs_out);
    std::fprintf(stderr, "[modeb] S=%lu E=%lu enumerated=%lu (std_both=%lu T_only_no_N=%lu N_only_no_T=%lu) direct=%lu cycle=%lu multi_path=%lu trunc=%lu no_path=%lu | tags: EQ=%lu D=%lu I=%lu S=%lu C=%lu | bytes=%.2fMB runs=%lu -> %s/modeb.bin\n",
                 (unsigned long)snodes.size(), (unsigned long)super_edges.size(),
                 (unsigned long)bubbles.size(), (unsigned long)n_std, (unsigned long)n_t_only, (unsigned long)n_n_only,
                 (unsigned long)direct_only, (unsigned long)cycle_skipped, (unsigned long)multi_path_C, (unsigned long)truncated_starts, (unsigned long)no_path_skipped,
                 (unsigned long)tag_count[0], (unsigned long)tag_count[1], (unsigned long)tag_count[2], (unsigned long)tag_count[3], (unsigned long)tag_count[4],
                 bytes_out.size() / 1e6, (unsigned long)runs_out.size(), dir.c_str());
}
