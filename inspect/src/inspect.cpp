// inspect — Stage 2 of inspect/.
//
// Given a target (bubble id / region / truth-SV id), extract the subgraph
// containing the target bubble(s) plus directly-adjacent color=3 anchors,
// and emit four files per subgraph:
//
//   <key>/raw.gfa       full subgraph in GFA1 (one node per LCP-core)
//   <key>/modeA.gfa     unitig-compacted (exact suffix/prefix merge)
//   <key>/modeB.gfa     allele-diff resolved (refinement on Mode A)
//   <key>/meta.json     annotations sidecar
//
// Usage:
//   inspect <dir> bubble <id>
//   inspect <dir> region <chr>:<start>-<end>
//   inspect <dir> sv     <truth_vcf> <sv_id>

#include "common.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace insp;

namespace {

constexpr int MIN_MULT = 3;

// ---------- helpers ----------
std::string sanitize(const std::string& s) {
    std::string out = s;
    for (char& c : out) if (c == ':' || c == '-' || c == '/' || c == ' ') c = '_';
    return out;
}

void mkdir_p(const std::string& p) {
    ::mkdir(p.c_str(), 0755);
}

void build_node_to_sig(const std::string& path, uint64_t K, std::vector<uint64_t>& node_to_sig) {
    FILE* f = xopen(path, "rb");
    uint64_t N; xread(f, &N, 8);
    node_to_sig.resize(K);
    uint64_t kept = 0;
    SigCount rec;
    for (uint64_t i = 0; i < N; i++) {
        xread(f, &rec, sizeof(rec));
        if (std::max(rec.n, rec.t) < (uint32_t)MIN_MULT) continue;
        if (kept >= K) die("counts.bin: more kept sigs than graph K");
        node_to_sig[kept++] = rec.sig;
    }
    std::fclose(f);
    if (kept != K) die("counts.bin kept != graph K");
}

// ---------- bubble enumeration (matches src/graph.cpp:267-283 + node_id tracking) ----------
struct Bubble {
    char type = 0;
    std::vector<uint32_t> node_ids;
    std::vector<uint32_t> anchor_ids;
};

void decompose_with_ids(const std::vector<std::vector<int32_t>>& adj,
                        const std::vector<uint8_t>& color,
                        std::vector<Bubble>& out) {
    const uint32_t K = (uint32_t)color.size();
    std::vector<char> vis(K, 0);
    std::queue<uint32_t> q;
    std::vector<uint32_t> anchor_buf;
    anchor_buf.reserve(64);
    for (uint32_t s = 0; s < K; s++) {
        if (color[s] == 3 || vis[s]) continue;
        q.push(s); vis[s] = 1;
        Bubble b; uint8_t seen = 0;
        anchor_buf.clear();
        while (!q.empty()) {
            uint32_t u = q.front(); q.pop();
            b.node_ids.push_back(u);
            seen |= color[u];
            for (int32_t v : adj[u]) {
                if (color[v] == 3) anchor_buf.push_back((uint32_t)v);
                else if (!vis[v]) { vis[v] = 1; q.push((uint32_t)v); }
            }
        }
        b.type = (seen == 1 ? 'D' : seen == 2 ? 'I' : 'C');
        std::sort(anchor_buf.begin(), anchor_buf.end());
        anchor_buf.erase(std::unique(anchor_buf.begin(), anchor_buf.end()), anchor_buf.end());
        b.anchor_ids = anchor_buf;
        out.push_back(std::move(b));
    }
}

// ---------- anchor analysis ----------
struct AnchorAnalysis {
    bool anchorable = false;
    bool mixed_strand = false;
    bool mixed_chrom = false;
    int32_t chrom_idx = -1;
    uint32_t L_pos = 0, R_pos = 0;
    uint8_t L_strand = 0, R_strand = 0;
    uint32_t L_node = UINT32_MAX, R_node = UINT32_MAX;
    std::string chrom_name;
};

void analyze_anchors(const Bubble& bub,
                     const std::vector<uint64_t>& node_to_sig,
                     const Posmap& pm,
                     AnchorAnalysis& aa) {
    struct Hit { uint16_t ci; uint32_t ps, pe; uint8_t strand; uint32_t node_id; };
    std::vector<Hit> hits;
    for (uint32_t aid : bub.anchor_ids) {
        uint64_t sig = node_to_sig[aid];
        auto it = pm.sig_pos.find(sig);
        if (it == pm.sig_pos.end()) continue;
        for (const Pos& p : it->second) hits.push_back({p.chrom_idx, p.pos_start, p.pos_end, p.strand, aid});
    }
    if (hits.empty()) return;

    std::map<uint16_t, std::vector<Hit>> by_chrom;
    for (const Hit& h : hits) by_chrom[h.ci].push_back(h);
    if (by_chrom.size() > 1) aa.mixed_chrom = true;
    auto best = by_chrom.begin();
    for (auto it = by_chrom.begin(); it != by_chrom.end(); ++it)
        if (it->second.size() > best->second.size()) best = it;
    auto& bh = best->second;
    std::map<uint8_t, int> by_strand;
    for (const Hit& h : bh) by_strand[h.strand]++;
    if (by_strand.size() > 1) aa.mixed_strand = true;
    uint8_t bs = 0; int bsc = -1;
    for (auto& kv : by_strand) if (kv.second > bsc) { bsc = kv.second; bs = kv.first; }

    uint32_t lo = UINT32_MAX, hi = 0, lo_node = UINT32_MAX, hi_node = UINT32_MAX;
    for (const Hit& h : bh) {
        if (h.strand != bs) continue;
        if (h.ps < lo) { lo = h.ps; lo_node = h.node_id; }
        if (h.pe > hi) { hi = h.pe; hi_node = h.node_id; }
    }
    if (lo == UINT32_MAX) return;
    if ((uint64_t)hi - (uint64_t)lo > 1'000'000ULL) return;
    aa.anchorable = true;
    aa.chrom_idx = (int32_t)best->first;
    aa.L_pos = lo; aa.R_pos = hi;
    aa.L_strand = bs; aa.R_strand = bs;
    aa.L_node = lo_node; aa.R_node = hi_node;
    if ((size_t)aa.chrom_idx < pm.chrom_names.size()) aa.chrom_name = pm.chrom_names[aa.chrom_idx];
}

// ---------- subgraph extraction ----------
struct SubNode {
    uint32_t node_id;            // original svp node id
    uint8_t  color;
    uint32_t ncov, tcov;
    std::string seq;             // raw sequence bytes
    // For mode A/B super-nodes:
    std::vector<uint32_t> members;  // original node ids that compacted into this super-node
    std::string label;           // GFA-friendly id
    // Posmap (for anchors only)
    bool has_pos = false;
    int32_t chrom_idx = -1;
    std::string chrom_name;
    uint32_t pos_start = 0, pos_end = 0;
    uint8_t strand = 0;
};

struct Subgraph {
    std::vector<SubNode> nodes;             // indexed 0..nodes.size()-1
    std::vector<std::pair<uint32_t, uint32_t>> edges;  // (i, j) into nodes[]
    std::unordered_map<uint32_t, uint32_t> orig_to_sub;  // orig node_id → sub index
    // Per-target-bubble info (for Mode B grouping)
    struct BubbleInfo {
        char type;
        std::vector<uint32_t> interior_sub_ids;     // sub-indexes (color != 3) inside the bubble
        std::vector<uint32_t> anchor_sub_ids;       // sub-indexes (color == 3) adjacent to interior
        uint32_t bubble_id_orig;
    };
    std::vector<BubbleInfo> bubbles;
    // Per-super-bubble Mode B classification (multiple bubbles sharing the same (L,R) anchor pair
    // are merged into one super-bubble before allele-difference resolution).
    struct SuperBubble {
        uint32_t L_sub_id, R_sub_id;
        std::vector<uint32_t> bubble_ids_orig;     // contributing bubble_id from the WG decompose
        std::vector<uint32_t> n_path;              // path sub-ids L→...→R through color=1 or 3
        std::vector<uint32_t> t_path;              // path sub-ids L→...→R through color=2 or 3
        std::string n_allele;                      // overlap-merged normal allele
        std::string t_allele;                      // overlap-merged tumor allele
        size_t shared_prefix = 0, shared_suffix = 0;
        uint32_t diff_n = 0, diff_t = 0;
        std::string sv_tag;                        // EQ / D / I / S / C / ?
        bool rc_applied = false;                   // T-allele was reverse-complemented (INV case)
    };
    std::vector<SuperBubble> super_bubbles;
};

// Position-indexed list of color=3 nodes per chrom (sorted by pos_start).
// Built once at startup; used by build_subgraph to include the dup-interior shared
// nodes whose posmap puts them between L_anchor and R_anchor.
using PosIndex = std::map<uint16_t, std::vector<std::pair<uint32_t, uint32_t>>>;

void build_pos_index(const std::vector<uint8_t>& color,
                     const std::vector<uint64_t>& node_to_sig,
                     const Posmap& pm,
                     PosIndex& idx) {
    const uint64_t K = color.size();
    for (uint64_t i = 0; i < K; i++) {
        if (color[i] != 3) continue;
        uint64_t sig = node_to_sig[i];
        auto it = pm.sig_pos.find(sig);
        if (it == pm.sig_pos.end()) continue;
        for (const Pos& p : it->second) {
            idx[p.chrom_idx].push_back({p.pos_start, (uint32_t)i});
        }
    }
    for (auto& kv : idx) std::sort(kv.second.begin(), kv.second.end());
}

void build_subgraph(const std::vector<Bubble>& bubs,
                    const std::vector<uint32_t>& target_bubble_ids,
                    const std::vector<std::vector<int32_t>>& adj,
                    const std::vector<uint8_t>& color,
                    const std::vector<uint32_t>& ncov,
                    const std::vector<uint32_t>& tcov,
                    const std::vector<uint64_t>& node_to_sig,
                    const Posmap& pm,
                    const Seqs& seqs,
                    const PosIndex& pos_index,
                    Subgraph& sg) {
    // Collect all node ids: interior of target bubbles + their color=3 neighbors
    std::unordered_set<uint32_t> nodes;
    for (uint32_t bid : target_bubble_ids) {
        const Bubble& b = bubs[bid];
        for (uint32_t n : b.node_ids) nodes.insert(n);
        for (uint32_t a : b.anchor_ids) nodes.insert(a);
    }

    // Expansion: for each target bubble, find its bounding genomic region via the
    // anchors' posmap entries, then include every color=3 node whose posmap puts
    // it inside that region. This captures the dup body's shared interior (whose
    // tcov/ncov ratio is the DUP signal) and any cycle-closing edges spanning the
    // dup. Critical guard: LCP-core sigs can hash to multiple ref positions
    // (murmur64 collisions), so we cluster the anchors' posmap entries per chrom
    // and take only the densest tight cluster (≤ 1 Mbp span). This rejects
    // LCP-collision phantom positions that would otherwise pull in the entire
    // genome between the collision and the real locus.
    constexpr uint32_t MAX_EXPANSION_SPAN = 1'000'000;
    for (uint32_t bid : target_bubble_ids) {
        const Bubble& b = bubs[bid];
        // Gather all anchor posmap entries
        std::map<uint16_t, std::vector<uint32_t>> per_chrom_positions;
        for (uint32_t a : b.anchor_ids) {
            uint64_t sig = node_to_sig[a];
            auto it = pm.sig_pos.find(sig);
            if (it == pm.sig_pos.end()) continue;
            for (const Pos& p : it->second) {
                per_chrom_positions[p.chrom_idx].push_back(p.pos_start);
            }
        }
        // Find best (chrom, lo, hi): chrom with the most positions; within that
        // chrom, the densest sliding window of width MAX_EXPANSION_SPAN.
        uint16_t best_chrom = UINT16_MAX;
        uint32_t best_count = 0;
        uint32_t best_lo = 0, best_hi = 0;
        for (auto& [ci, pos_vec] : per_chrom_positions) {
            std::sort(pos_vec.begin(), pos_vec.end());
            // Sliding window: for each i, find the largest j such that pos[j] - pos[i] <= span
            for (size_t i = 0; i < pos_vec.size(); i++) {
                size_t j = i;
                while (j + 1 < pos_vec.size() && pos_vec[j + 1] - pos_vec[i] <= MAX_EXPANSION_SPAN) j++;
                uint32_t count = (uint32_t)(j - i + 1);
                if (count > best_count) {
                    best_count = count;
                    best_chrom = ci;
                    best_lo = pos_vec[i];
                    best_hi = pos_vec[j];
                }
            }
        }
        if (best_chrom == UINT16_MAX) continue;
        auto pit = pos_index.find(best_chrom);
        if (pit == pos_index.end()) continue;
        const auto& v = pit->second;
        auto lo_it = std::lower_bound(v.begin(), v.end(),
                                      std::pair<uint32_t,uint32_t>{best_lo, 0});
        for (auto it = lo_it; it != v.end() && it->first <= best_hi; ++it) {
            nodes.insert(it->second);
        }
    }

    // Create SubNode entries
    std::vector<uint32_t> sorted_nodes(nodes.begin(), nodes.end());
    std::sort(sorted_nodes.begin(), sorted_nodes.end());
    sg.nodes.reserve(sorted_nodes.size());
    for (uint32_t nid : sorted_nodes) {
        SubNode sn;
        sn.node_id = nid;
        sn.color = color[nid];
        sn.ncov = ncov[nid];
        sn.tcov = tcov[nid];
        uint16_t L; const uint8_t* p = seqs.seq(nid, L);
        sn.seq.assign((const char*)p, (size_t)L);
        sn.members.push_back(nid);
        sn.label = "n" + std::to_string(nid);
        // Posmap for color=3 nodes
        if (color[nid] == 3) {
            uint64_t sig = node_to_sig[nid];
            auto it = pm.sig_pos.find(sig);
            if (it != pm.sig_pos.end() && !it->second.empty()) {
                const Pos& p0 = it->second[0];
                sn.has_pos = true;
                sn.chrom_idx = (int32_t)p0.chrom_idx;
                if ((size_t)p0.chrom_idx < pm.chrom_names.size())
                    sn.chrom_name = pm.chrom_names[p0.chrom_idx];
                sn.pos_start = p0.pos_start;
                sn.pos_end = p0.pos_end;
                sn.strand = p0.strand;
            }
        }
        uint32_t sub_idx = (uint32_t)sg.nodes.size();
        sg.nodes.push_back(std::move(sn));
        sg.orig_to_sub[nid] = sub_idx;
    }

    // Edges (only between nodes both in subgraph)
    std::set<std::pair<uint32_t, uint32_t>> eset;
    for (uint32_t nid : sorted_nodes) {
        uint32_t sub_a = sg.orig_to_sub[nid];
        for (int32_t v : adj[nid]) {
            auto it = sg.orig_to_sub.find((uint32_t)v);
            if (it == sg.orig_to_sub.end()) continue;
            uint32_t sub_b = it->second;
            if (sub_a < sub_b) eset.emplace(sub_a, sub_b);
            else if (sub_b < sub_a) eset.emplace(sub_b, sub_a);
        }
    }
    sg.edges.assign(eset.begin(), eset.end());

    // Bubble info per target
    sg.bubbles.reserve(target_bubble_ids.size());
    for (uint32_t bid : target_bubble_ids) {
        const Bubble& b = bubs[bid];
        Subgraph::BubbleInfo bi;
        bi.type = b.type;
        bi.bubble_id_orig = bid;
        for (uint32_t n : b.node_ids) bi.interior_sub_ids.push_back(sg.orig_to_sub[n]);
        for (uint32_t a : b.anchor_ids) bi.anchor_sub_ids.push_back(sg.orig_to_sub[a]);
        sg.bubbles.push_back(std::move(bi));
    }
}

// ---------- Mode A: unitig compaction with EXACT suffix/prefix overlap merge ----------
size_t suffix_prefix_overlap(const std::string& A, const std::string& B) {
    size_t maxk = std::min(A.size(), B.size());
    for (size_t k = maxk; k > 0; k--) {
        if (A.compare(A.size() - k, k, B, 0, k) == 0) return k;
    }
    return 0;
}

void merge_node_into_chain(std::string& chain, std::vector<uint32_t>& members,
                           uint32_t& ncov_sum, uint32_t& tcov_sum,
                           const SubNode& next, std::vector<std::string>& warnings,
                           bool& merge_ok) {
    // Bidirectional overlap: the chain walk picks a direction arbitrarily, but the
    // genome-order of the two LCP-cores can be either way. Try both A.suffix-B.prefix
    // AND B.suffix-A.prefix; take whichever gives the longer exact match.
    const size_t kf = suffix_prefix_overlap(chain, next.seq);
    const size_t kr = suffix_prefix_overlap(next.seq, chain);
    if (kf == 0 && kr == 0) {
        merge_ok = false;
        std::ostringstream os;
        os << "no_overlap between members ";
        for (uint32_t m : members) os << m << ",";
        os << " and node " << next.node_id;
        warnings.push_back(os.str());
        return;
    }
    if (kf >= kr) {
        chain.append(next.seq, kf, next.seq.size() - kf);
    } else {
        // Reverse-order: next.seq comes BEFORE the current chain in the genome.
        // Prepend next.seq to chain (using kr-byte overlap from next's suffix into chain's prefix).
        std::string new_chain = next.seq + chain.substr(kr);
        chain = std::move(new_chain);
    }
    members.push_back(next.node_id);
    ncov_sum += next.ncov;
    tcov_sum += next.tcov;
    merge_ok = true;
}

void run_mode_a(const Subgraph& raw, Subgraph& mA, std::vector<std::string>& warnings) {
    // adjacency on raw
    const size_t N = raw.nodes.size();
    std::vector<std::vector<uint32_t>> adj(N);
    for (auto& e : raw.edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    // mark which raw-node indexes belong to which super-node
    std::vector<int32_t> super_of(N, -1);
    int32_t next_super = 0;

    // Find chains: walk forward then backward from each unvisited candidate (deg=2)
    // until we hit a branching or a color-change boundary.
    auto can_extend = [&](uint32_t cur, uint32_t next) {
        // Both must be in the same chain: 1-in/1-out semantics implies deg<=2 + color match
        return raw.nodes[cur].color == raw.nodes[next].color &&
               adj[next].size() <= 2;
    };

    for (uint32_t start = 0; start < N; start++) {
        if (super_of[start] != -1) continue;
        // Initialize the chain with this node, then extend in both directions
        std::vector<uint32_t> chain_nodes; chain_nodes.push_back(start);
        super_of[start] = next_super;
        // Walk forward
        uint32_t cur = start;
        while (true) {
            if (adj[cur].size() > 2) break;          // branch (interior node must be 1-in / 1-out total deg ≤ 2)
            uint32_t prev = chain_nodes.size() >= 2 ? chain_nodes[chain_nodes.size() - 2] : UINT32_MAX;
            uint32_t next = UINT32_MAX;
            for (uint32_t n : adj[cur]) if (n != prev) { next = n; break; }
            if (next == UINT32_MAX) break;
            if (super_of[next] != -1) break;
            if (!can_extend(cur, next)) break;
            // also the predecessor side: next must have deg<=2 too (which can_extend covers)
            chain_nodes.push_back(next);
            super_of[next] = next_super;
            cur = next;
        }
        // Walk backward from start
        cur = start;
        while (true) {
            if (adj[cur].size() > 2) break;
            // Find the predecessor — the neighbor that's NOT in chain_nodes (we go backwards)
            uint32_t in_chain_succ = chain_nodes.size() >= 2 && chain_nodes[1] != UINT32_MAX
                                     ? chain_nodes[chain_nodes[0] == cur ? 1 : 0]
                                     : UINT32_MAX;
            uint32_t back = UINT32_MAX;
            for (uint32_t n : adj[cur]) if (n != in_chain_succ) { back = n; break; }
            if (back == UINT32_MAX) break;
            if (super_of[back] != -1) break;
            if (!can_extend(cur, back)) break;
            chain_nodes.insert(chain_nodes.begin(), back);
            super_of[back] = next_super;
            cur = back;
        }

        // Build super-node by merging in chain order
        SubNode sn;
        sn.color = raw.nodes[start].color;
        sn.ncov = 0; sn.tcov = 0;
        sn.seq = raw.nodes[chain_nodes[0]].seq;
        sn.members.push_back(raw.nodes[chain_nodes[0]].node_id);
        sn.ncov += raw.nodes[chain_nodes[0]].ncov;
        sn.tcov += raw.nodes[chain_nodes[0]].tcov;
        // For posmap on color=3 super-nodes: take leftmost member's start, rightmost member's end
        if (sn.color == 3 && raw.nodes[chain_nodes[0]].has_pos) {
            sn.has_pos = true;
            sn.chrom_idx = raw.nodes[chain_nodes[0]].chrom_idx;
            sn.chrom_name = raw.nodes[chain_nodes[0]].chrom_name;
            sn.pos_start = raw.nodes[chain_nodes[0]].pos_start;
            sn.pos_end = raw.nodes[chain_nodes[0]].pos_end;
            sn.strand = raw.nodes[chain_nodes[0]].strand;
        }
        for (size_t i = 1; i < chain_nodes.size(); i++) {
            const SubNode& nxt = raw.nodes[chain_nodes[i]];
            bool ok = false;
            std::vector<uint32_t> members_tmp = sn.members;
            uint32_t nc = sn.ncov, tc = sn.tcov;
            std::string chain_tmp = sn.seq;
            merge_node_into_chain(chain_tmp, members_tmp, nc, tc, nxt, warnings, ok);
            if (!ok) {
                // Break the chain here: this and subsequent nodes form a new super-node
                // (we mark super_of accordingly and start a new chain via the outer loop)
                for (size_t j = i; j < chain_nodes.size(); j++) super_of[chain_nodes[j]] = -1;
                chain_nodes.resize(i);
                break;
            }
            sn.seq = std::move(chain_tmp);
            sn.members = std::move(members_tmp);
            sn.ncov = nc;
            sn.tcov = tc;
            if (sn.has_pos && nxt.has_pos &&
                nxt.chrom_idx == sn.chrom_idx && nxt.strand == sn.strand) {
                if (nxt.pos_start < sn.pos_start) sn.pos_start = nxt.pos_start;
                if (nxt.pos_end > sn.pos_end)     sn.pos_end = nxt.pos_end;
            }
        }
        sn.label = "u" + std::to_string(next_super);
        sn.node_id = raw.nodes[chain_nodes[0]].node_id;  // representative (first member)
        mA.nodes.push_back(std::move(sn));
        next_super++;
    }

    // Build edges between super-nodes: any raw edge whose endpoints are in different supers
    std::set<std::pair<uint32_t, uint32_t>> eset;
    for (const auto& e : raw.edges) {
        int32_t sa = super_of[e.first], sb = super_of[e.second];
        if (sa < 0 || sb < 0 || sa == sb) continue;
        uint32_t a = (uint32_t)sa, b = (uint32_t)sb;
        if (a < b) eset.emplace(a, b);
        else        eset.emplace(b, a);
    }
    mA.edges.assign(eset.begin(), eset.end());

    // Map orig→super index in mA
    for (size_t super = 0; super < mA.nodes.size(); super++) {
        for (uint32_t m : mA.nodes[super].members) mA.orig_to_sub[m] = (uint32_t)super;
    }

    // Translate bubble info from raw to mA
    mA.bubbles.reserve(raw.bubbles.size());
    for (const Subgraph::BubbleInfo& bi : raw.bubbles) {
        Subgraph::BubbleInfo nbi;
        nbi.type = bi.type;
        nbi.bubble_id_orig = bi.bubble_id_orig;
        std::unordered_set<uint32_t> intset, ancset;
        for (uint32_t s : bi.interior_sub_ids) {
            // s here is a raw sub-index; translate via members
            uint32_t orig = raw.nodes[s].node_id;
            intset.insert(mA.orig_to_sub[orig]);
        }
        for (uint32_t s : bi.anchor_sub_ids) {
            uint32_t orig = raw.nodes[s].node_id;
            ancset.insert(mA.orig_to_sub[orig]);
        }
        nbi.interior_sub_ids.assign(intset.begin(), intset.end());
        nbi.anchor_sub_ids.assign(ancset.begin(), ancset.end());
        std::sort(nbi.interior_sub_ids.begin(), nbi.interior_sub_ids.end());
        std::sort(nbi.anchor_sub_ids.begin(), nbi.anchor_sub_ids.end());
        mA.bubbles.push_back(std::move(nbi));
    }
}

// ---------- Mode B: refinement on Mode A ----------
size_t lcp(const std::string& a, const std::string& b) {
    size_t k = 0, m = std::min(a.size(), b.size());
    while (k < m && a[k] == b[k]) k++;
    return k;
}
size_t lcs(const std::string& a, const std::string& b) {
    size_t k = 0, m = std::min(a.size(), b.size());
    while (k < m && a[a.size() - 1 - k] == b[b.size() - 1 - k]) k++;
    return k;
}

std::string reverse_complement(const std::string& s) {
    std::string r(s.size(), 'N');
    for (size_t i = 0; i < s.size(); i++) {
        char c = s[s.size() - 1 - i];
        r[i] = (c == 'A') ? 'T' : (c == 'T') ? 'A' : (c == 'C') ? 'G' : (c == 'G') ? 'C' : 'N';
    }
    return r;
}

// Bidirectional overlap concat: try both A.suffix-B.prefix and B.suffix-A.prefix,
// take whichever finds the longer exact byte match. If neither finds a positive
// overlap, fall back to plain concat (logged as a no-overlap warning).
std::string concat_bidir(const std::string& A, const std::string& B,
                         const std::string& ctx, std::vector<std::string>& warnings) {
    if (A.empty()) return B;
    if (B.empty()) return A;
    size_t kf = suffix_prefix_overlap(A, B);
    size_t kr = suffix_prefix_overlap(B, A);
    if (kf == 0 && kr == 0) {
        warnings.push_back("concat no-overlap (" + ctx + ")");
        return A + B;
    }
    if (kf >= kr) return A + B.substr(kf);
    // Reverse-order win: B comes before A in the genome direction
    return B + A.substr(kr);
}

// BFS L→R through allowed nodes. allowed_color: 1 (normal+anchors) or 2 (tumor+anchors).
// Anchors (color=3) are always allowed; interior nodes must match allowed_color.
std::vector<uint32_t> walk_color_path(uint32_t L_sub, uint32_t R_sub,
                                     const Subgraph& sg,
                                     const std::vector<std::vector<uint32_t>>& adj,
                                     uint8_t allowed_color) {
    if (L_sub == R_sub) return {L_sub};
    const size_t N = sg.nodes.size();
    std::vector<int32_t> parent(N, -1);
    parent[L_sub] = (int32_t)L_sub;
    std::queue<uint32_t> q; q.push(L_sub);
    while (!q.empty()) {
        uint32_t u = q.front(); q.pop();
        if (u == R_sub) break;
        for (uint32_t v : adj[u]) {
            if (parent[v] != -1) continue;
            if (v != R_sub) {
                uint8_t c = sg.nodes[v].color;
                if (c != 3 && c != allowed_color) continue;
            }
            parent[v] = (int32_t)u;
            q.push(v);
        }
    }
    std::vector<uint32_t> path;
    if (parent[R_sub] == -1) return path;  // no valid path
    for (uint32_t u = R_sub; u != L_sub; u = (uint32_t)parent[u]) path.push_back(u);
    path.push_back(L_sub);
    std::reverse(path.begin(), path.end());
    return path;
}

// Merge a path of super-node sequences in WALK ORDER (L stays at the start, R at
// the end). Forward overlap only: `A.suffix vs B.prefix`. If forward overlap is 0,
// we also try RC(B) — adjacent cores can have flipped strand-orientation bytes if
// they originated from opposite-strand reads. If RC also fails, plain concat with
// warning. Never swap A/B order (that would corrupt the L→R orientation Mode B
// relies on for the trim).
std::string merge_path(const std::vector<uint32_t>& path, const Subgraph& sg,
                       const std::string& ctx, std::vector<std::string>& warnings) {
    if (path.empty()) return "";
    std::string merged = sg.nodes[path[0]].seq;
    for (size_t i = 1; i < path.size(); i++) {
        const std::string& nxt = sg.nodes[path[i]].seq;
        if (merged.empty()) { merged = nxt; continue; }
        if (nxt.empty())    continue;
        size_t kf = suffix_prefix_overlap(merged, nxt);
        if (kf > 0) { merged.append(nxt, kf, nxt.size() - kf); continue; }
        std::string nxt_rc = reverse_complement(nxt);
        size_t krc = suffix_prefix_overlap(merged, nxt_rc);
        if (krc > 0) {
            merged.append(nxt_rc, krc, nxt_rc.size() - krc);
            warnings.push_back("merge_path RC-joined (" + ctx + ")");
            continue;
        }
        warnings.push_back("merge_path no-overlap forward (" + ctx + ")");
        merged += nxt;
    }
    return merged;
}

// Mode B v2: group bubbles sharing the same (L, R) anchor pair into super-bubbles,
// walk N-path and T-path through color-constrained interior, trim shared prefix/suffix,
// absorb into the anchor super-nodes in place, classify EQ/D/I/S/C.
void run_mode_b(Subgraph& mA, Subgraph& mB, std::vector<std::string>& warnings) {
    mB = mA;  // deep copy; mutate anchors + collapse interior in place
    const size_t N = mB.nodes.size();
    std::vector<std::vector<uint32_t>> adj(N);
    for (auto& e : mB.edges) { adj[e.first].push_back(e.second); adj[e.second].push_back(e.first); }

    // Pre-compute (L, R) for each input bubble and group by anchor pair.
    std::map<std::pair<uint32_t, uint32_t>, Subgraph::SuperBubble> groups;
    for (Subgraph::BubbleInfo& bi : mB.bubbles) {
        uint32_t L = UINT32_MAX, R = UINT32_MAX;
        uint32_t L_pos = UINT32_MAX, R_pos = 0;
        for (uint32_t a : bi.anchor_sub_ids) {
            const SubNode& sn = mB.nodes[a];
            if (!sn.has_pos) continue;
            if (sn.pos_start < L_pos) { L_pos = sn.pos_start; L = a; }
            if (sn.pos_end > R_pos)   { R_pos = sn.pos_end;   R = a; }
        }
        if (L == UINT32_MAX || R == UINT32_MAX || L == R) {
            std::ostringstream os; os << "no usable anchor pair for bubble " << bi.bubble_id_orig;
            warnings.push_back(os.str());
            continue;
        }
        auto& sb = groups[{L, R}];
        sb.L_sub_id = L;
        sb.R_sub_id = R;
        sb.bubble_ids_orig.push_back(bi.bubble_id_orig);
    }

    // Classify each super-bubble.
    for (auto& [key, sb] : groups) {
        const std::string ctx = "sb L=" + std::to_string(sb.L_sub_id) + " R=" + std::to_string(sb.R_sub_id);

        sb.n_path = walk_color_path(sb.L_sub_id, sb.R_sub_id, mB, adj, /*allowed=*/1);
        sb.t_path = walk_color_path(sb.L_sub_id, sb.R_sub_id, mB, adj, /*allowed=*/2);

        // Walk fallback: if no color=1 path exists (e.g., DUP), the normal allele is just L+R
        // (direct anchor-to-anchor) — reconstructed by overlap-concat below.
        std::string nA = sb.n_path.empty()
            ? concat_bidir(mB.nodes[sb.L_sub_id].seq, mB.nodes[sb.R_sub_id].seq, ctx + " n-fallback", warnings)
            : merge_path(sb.n_path, mB, ctx + " n", warnings);
        std::string tA = sb.t_path.empty()
            ? concat_bidir(mB.nodes[sb.L_sub_id].seq, mB.nodes[sb.R_sub_id].seq, ctx + " t-fallback", warnings)
            : merge_path(sb.t_path, mB, ctx + " t", warnings);

        // INV candidate: anchor strands disagree → RC the tumor allele before trim
        if (mB.nodes[sb.L_sub_id].has_pos && mB.nodes[sb.R_sub_id].has_pos &&
            mB.nodes[sb.L_sub_id].strand != mB.nodes[sb.R_sub_id].strand) {
            tA = reverse_complement(tA);
            sb.rc_applied = true;
        }

        size_t p = lcp(nA, tA);
        size_t s = (p < nA.size() && p < tA.size())
            ? lcs(nA.substr(p), tA.substr(p))
            : 0;
        sb.shared_prefix = p;
        sb.shared_suffix = s;
        sb.n_allele = nA;
        sb.t_allele = tA;

        std::string dN = (p + s <= nA.size()) ? nA.substr(p, nA.size() - p - s) : "";
        std::string dT = (p + s <= tA.size()) ? tA.substr(p, tA.size() - p - s) : "";
        sb.diff_n = (uint32_t)dN.size();
        sb.diff_t = (uint32_t)dT.size();

        if (dN.empty() && dT.empty())                       sb.sv_tag = "EQ";
        else if (dT.empty())                                 sb.sv_tag = "D";
        else if (dN.empty())                                 sb.sv_tag = "I";
        else if (dN.size() == 1 && dT.size() == 1)           sb.sv_tag = "S";
        else                                                 sb.sv_tag = "C";

        // Absorb shared prefix into L_anchor, shared suffix into R_anchor.
        // The new L_anchor.seq is nA[0..p) (== tA[0..p)); R_anchor.seq is nA[len-s..).
        mB.nodes[sb.L_sub_id].seq = nA.substr(0, p);
        mB.nodes[sb.R_sub_id].seq = nA.substr(nA.size() - s, s);

        // Shrink/collapse interior super-nodes that participated in this super-bubble.
        // For each interior super-node in n_path (excluding L,R), replace its seq with the
        // portion of dN it should own. v1 lumps the whole dN into the first interior
        // super-node and collapses the rest; v2 (TODO) could divide proportionally.
        auto collapse_interior = [&](std::vector<uint32_t>& path, const std::string& diff) {
            bool first = true;
            for (uint32_t sid : path) {
                if (sid == sb.L_sub_id || sid == sb.R_sub_id) continue;
                if (first && !diff.empty()) {
                    mB.nodes[sid].seq = diff;
                    if (mB.nodes[sid].label.find("_d") == std::string::npos)
                        mB.nodes[sid].label += "_d";
                    first = false;
                } else {
                    mB.nodes[sid].seq.clear();
                    if (mB.nodes[sid].label.find("_collapsed") == std::string::npos)
                        mB.nodes[sid].label += "_collapsed";
                }
            }
        };
        collapse_interior(sb.n_path, dN);
        collapse_interior(sb.t_path, dT);

        mB.super_bubbles.push_back(sb);
    }
}

// ---------- GFA emission ----------
const char* color_label(uint8_t c) {
    return c == 1 ? "NORMAL" : c == 2 ? "TUMOR" : c == 3 ? "SHARED" : "UNK";
}
const char* color_hex(uint8_t c) {
    return c == 1 ? "#5B9BD5" : c == 2 ? "#E07B5C" : c == 3 ? "#8E8E8E" : "#000000";
}

void emit_gfa(const Subgraph& sg, const std::string& path, bool skip_empty) {
    FILE* f = xopen(path, "w");
    std::fputs("H\tVN:Z:1.0\n", f);
    for (size_t i = 0; i < sg.nodes.size(); i++) {
        const SubNode& sn = sg.nodes[i];
        if (skip_empty && sn.seq.empty()) continue;  // collapsed segment in Mode B
        std::fprintf(f, "S\t%s\t%s\tLN:i:%zu\tCL:Z:%s\tCC:Z:%s\tKC:i:%u\tRC:i:%u\tNN:i:%zu\tNI:i:%u",
                     sn.label.c_str(),
                     sn.seq.empty() ? "*" : sn.seq.c_str(),
                     sn.seq.size(), color_label(sn.color), color_hex(sn.color),
                     sn.ncov, sn.tcov, sn.members.size(), sn.node_id);
        if (sn.has_pos) std::fprintf(f, "\tPS:Z:%s:%u-%u:%u", sn.chrom_name.empty() ? "?" : sn.chrom_name.c_str(), sn.pos_start, sn.pos_end, sn.strand);
        // Depth ratio (DUP signal): tumor coverage / normal coverage, averaged across members.
        // Shared nodes with DR ≥ 1.5 are dup-region candidates; DR ≥ 1.8 is strong evidence.
        if (sn.color == 3 && sn.ncov > 0) {
            double dr = (double)sn.tcov / (double)sn.ncov;
            std::fprintf(f, "\tDR:f:%.3f", dr);
        }
        std::fputc('\n', f);
    }
    for (const auto& e : sg.edges) {
        // skip edges to collapsed (empty-seq) segments under Mode B
        if (skip_empty && (sg.nodes[e.first].seq.empty() || sg.nodes[e.second].seq.empty())) continue;
        std::fprintf(f, "L\t%s\t+\t%s\t+\t0M\n",
                     sg.nodes[e.first].label.c_str(),
                     sg.nodes[e.second].label.c_str());
    }
    std::fclose(f);
}

// ---------- cycle detection (DFS back-edges in the subgraph) ----------
struct Cycle { std::vector<uint32_t> path; };

std::vector<Cycle> find_cycles(const Subgraph& sg) {
    const size_t N = sg.nodes.size();
    std::vector<std::vector<uint32_t>> adj(N);
    for (auto& e : sg.edges) { adj[e.first].push_back(e.second); adj[e.second].push_back(e.first); }

    std::vector<Cycle> cycles;
    std::vector<char> visited(N, 0);
    std::vector<int32_t> parent(N, -1);
    std::vector<uint32_t> depth(N, 0);

    for (uint32_t start = 0; start < N; start++) {
        if (visited[start]) continue;
        // Iterative DFS
        std::vector<std::pair<uint32_t, size_t>> stack;
        stack.push_back({start, 0});
        visited[start] = 1;
        while (!stack.empty()) {
            auto& [u, k] = stack.back();
            if (k < adj[u].size()) {
                uint32_t v = adj[u][k++];
                if ((int32_t)v == parent[u]) continue;  // skip the edge we came in on
                if (visited[v]) {
                    // Back-edge: found a cycle. Reconstruct path v → ... → u → v.
                    Cycle c;
                    c.path.push_back(v);
                    for (uint32_t x = u; x != v && (int32_t)x != -1;) {
                        c.path.push_back(x);
                        if (parent[x] == -1) break;
                        x = (uint32_t)parent[x];
                    }
                    // Dedupe: only keep small unique cycles (limit per call)
                    if (c.path.size() >= 3 && c.path.size() <= 30 && cycles.size() < 16) {
                        cycles.push_back(std::move(c));
                    }
                    continue;
                }
                visited[v] = 1;
                parent[v] = (int32_t)u;
                depth[v] = depth[u] + 1;
                stack.push_back({v, 0});
            } else {
                stack.pop_back();
            }
        }
    }
    return cycles;
}

// ---------- meta.json ----------
void emit_meta(const Subgraph& raw, const Subgraph& mA, const Subgraph& mB,
               const std::vector<Bubble>& bubs,
               const std::vector<uint32_t>& target_bubble_ids,
               const Posmap& pm,
               const std::vector<std::string>& warnings,
               const std::string& key,
               const std::string& path,
               const std::vector<std::pair<std::string, std::string>>& extra_kv) {
    FILE* f = xopen(path, "w");
    std::fprintf(f, "{\n");
    std::fprintf(f, "  \"key\": \"%s\",\n", key.c_str());
    std::fprintf(f, "  \"subgraph_size\": {\"raw_nodes\": %zu, \"modeA_nodes\": %zu, \"modeB_segments\": %zu},\n",
                 raw.nodes.size(), mA.nodes.size(), mB.nodes.size());
    std::fprintf(f, "  \"target_bubble_ids\": [");
    for (size_t i = 0; i < target_bubble_ids.size(); i++)
        std::fprintf(f, "%s%u", i ? ", " : "", target_bubble_ids[i]);
    std::fprintf(f, "],\n");

    // Per-super-bubble Mode B classification
    // (super-bubble = group of bubbles sharing the same (L, R) anchor pair, classified together)
    std::fprintf(f, "  \"modeB_calls\": [\n");
    for (size_t i = 0; i < mB.super_bubbles.size(); i++) {
        const auto& sb = mB.super_bubbles[i];
        // DUP evidence: tcov/ncov of L, R anchors + per-path shared-node DR stats.
        // L_anchor and R_anchor reside at sub-ids L_sub_id and R_sub_id in mB.
        const auto& L = mB.nodes[sb.L_sub_id];
        const auto& R = mB.nodes[sb.R_sub_id];
        double L_dr = L.ncov > 0 ? (double)L.tcov / (double)L.ncov : -1.0;
        double R_dr = R.ncov > 0 ? (double)R.tcov / (double)R.ncov : -1.0;
        // Mean DR across the shared interior portion of the n-path
        // Also: DUP size estimate = sum(len * (DR - 1)) for shared nodes with DR > 1.
        // For a clean tandem-DUP of length L bp, every base in the dup region has
        // tcov ≈ 2 * ncov (DR = 2), so the sum equals L.
        double n_path_dr_sum = 0.0; int n_path_dr_count = 0;
        double dup_size_estimate = 0.0;
        for (uint32_t sid : sb.n_path) {
            const auto& sn = mB.nodes[sid];
            if (sn.color == 3 && sn.ncov > 0) {
                double dr = (double)sn.tcov / (double)sn.ncov;
                n_path_dr_sum += dr;
                n_path_dr_count++;
                if (dr > 1.0) dup_size_estimate += (double)sn.seq.size() * (dr - 1.0);
            }
        }
        double n_path_mean_dr = n_path_dr_count ? n_path_dr_sum / n_path_dr_count : -1.0;
        // DUP hint: if no allele-level diff (sv_tag=="EQ") but the shared interior
        // shows elevated tumor depth, the bubble is most likely a tandem dup
        // whose junction sequence is captured elsewhere in the subgraph.
        bool dup_hint = (sb.sv_tag == "EQ" && dup_size_estimate >= 50.0);

        std::fprintf(f, "    {\"L_sub\": %u, \"R_sub\": %u, \"sv_tag\": \"%s\","
                     " \"diff_n\": %u, \"diff_t\": %u,"
                     " \"n_path_len\": %zu, \"t_path_len\": %zu,"
                     " \"shared_prefix\": %zu, \"shared_suffix\": %zu,"
                     " \"rc_applied\": %s,"
                     " \"L_anchor_dr\": %.3f, \"R_anchor_dr\": %.3f,"
                     " \"n_path_shared_mean_dr\": %.3f,"
                     " \"dup_size_estimate_bp\": %.0f,"
                     " \"dup_hint\": %s,"
                     " \"bubble_ids\": [",
                     sb.L_sub_id, sb.R_sub_id, sb.sv_tag.c_str(),
                     sb.diff_n, sb.diff_t,
                     sb.n_path.size(), sb.t_path.size(),
                     sb.shared_prefix, sb.shared_suffix,
                     sb.rc_applied ? "true" : "false",
                     L_dr, R_dr, n_path_mean_dr,
                     dup_size_estimate,
                     dup_hint ? "true" : "false");
        for (size_t j = 0; j < sb.bubble_ids_orig.size(); j++)
            std::fprintf(f, "%s%u", j ? "," : "", sb.bubble_ids_orig[j]);
        std::fprintf(f, "]}%s\n", i + 1 < mB.super_bubbles.size() ? "," : "");
    }
    std::fprintf(f, "  ],\n");

    // Cycles in the raw subgraph (key signal for DUP/INV: tandem-DUPs create a
    // junction edge from end-of-copy back to start-of-copy, which closes a loop
    // through the shared interior).
    auto raw_cycles = find_cycles(raw);
    std::fprintf(f, "  \"cycles\": [");
    for (size_t i = 0; i < raw_cycles.size(); i++) {
        std::fprintf(f, "%s[", i ? ", " : "");
        for (size_t j = 0; j < raw_cycles[i].path.size(); j++)
            std::fprintf(f, "%s\"%s\"", j ? "," : "", raw.nodes[raw_cycles[i].path[j]].label.c_str());
        std::fprintf(f, "]");
    }
    std::fprintf(f, "],\n");

    // Per-input-bubble grouping info (which super-bubble each contributing bubble belongs to)
    std::fprintf(f, "  \"input_bubbles\": [\n");
    for (size_t i = 0; i < mB.bubbles.size(); i++) {
        const auto& bi = mB.bubbles[i];
        const Bubble& orig = bubs[bi.bubble_id_orig];
        std::fprintf(f, "    {\"bubble_id\": %u, \"type\": \"%c\","
                     " \"n_nodes_raw\": %zu, \"n_anchors_raw\": %zu}%s\n",
                     bi.bubble_id_orig, bi.type,
                     orig.node_ids.size(), orig.anchor_ids.size(),
                     i + 1 < mB.bubbles.size() ? "," : "");
    }
    std::fprintf(f, "  ],\n");

    // Extra key-value pairs
    for (size_t i = 0; i < extra_kv.size(); i++) {
        std::fprintf(f, "  \"%s\": %s%s\n",
                     extra_kv[i].first.c_str(), extra_kv[i].second.c_str(),
                     i + 1 < extra_kv.size() ? "," : "");
    }

    std::fprintf(f, "  ,\n  \"warnings\": [");
    for (size_t i = 0; i < warnings.size(); i++) {
        std::fprintf(f, "%s\"%s\"", i ? ", " : "", warnings[i].c_str());
    }
    std::fprintf(f, "]\n");
    std::fprintf(f, "}\n");
    std::fclose(f);
    (void)pm;
}

// ---------- bubble-pickers ----------
struct TruthRec { std::string chrom; uint32_t pos, end; int32_t svlen; std::string svtype, sv_id; };
void load_truth(const std::string& path, std::vector<TruthRec>& out) {
    std::ifstream fs(path);
    if (!fs) die("truth VCF open");
    std::string line;
    while (std::getline(fs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> t;
        size_t b = 0;
        for (size_t i = 0; i <= line.size(); i++) if (i == line.size() || line[i] == '\t') { t.push_back(line.substr(b, i - b)); b = i + 1; }
        if (t.size() < 8) continue;
        TruthRec r{};
        r.chrom = t[0]; r.pos = (uint32_t)std::stoul(t[1]) - 1; r.sv_id = t[2];
        std::map<std::string, std::string> kv;
        size_t kvb = 0;
        for (size_t i = 0; i <= t[7].size(); i++) if (i == t[7].size() || t[7][i] == ';') { auto p = t[7].substr(kvb, i - kvb); auto eq = p.find('='); if (eq != std::string::npos) kv[p.substr(0, eq)] = p.substr(eq+1); kvb = i+1; }
        r.svtype = kv.count("SVTYPE") ? kv["SVTYPE"] : "";
        if (r.svtype.empty()) continue;
        r.svlen = kv.count("SVLEN") ? std::stoi(kv["SVLEN"]) : 0;
        r.end = kv.count("END") ? (uint32_t)std::stoul(kv["END"]) - 1 : r.pos;
        out.push_back(r);
    }
}

// Pick bubble(s) whose anchors fall within [lo - 500, hi + 500] on the given chrom.
void pick_bubbles_by_region(const std::vector<Bubble>& bubs,
                            const std::vector<uint64_t>& node_to_sig,
                            const Posmap& pm,
                            const std::string& chrom, uint32_t lo, uint32_t hi,
                            std::vector<uint32_t>& out) {
    // find chrom_idx
    int chrom_idx = -1;
    for (size_t i = 0; i < pm.chrom_names.size(); i++) if (pm.chrom_names[i] == chrom) { chrom_idx = (int)i; break; }
    if (chrom_idx < 0) { std::fprintf(stderr, "[inspect] chrom %s not in posmap\n", chrom.c_str()); return; }

    for (uint32_t bid = 0; bid < bubs.size(); bid++) {
        AnchorAnalysis aa;
        analyze_anchors(bubs[bid], node_to_sig, pm, aa);
        if (!aa.anchorable) continue;
        if (aa.chrom_idx != chrom_idx) continue;
        // overlap test on [aa.L_pos, aa.R_pos] vs [lo, hi]
        if (aa.R_pos < lo || aa.L_pos > hi) continue;
        out.push_back(bid);
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr,
                     "usage:\n"
                     "  %s <dir> bubble <id>\n"
                     "  %s <dir> region <chr>:<start>-<end>\n"
                     "  %s <dir> sv     <truth_vcf> <sv_id>\n",
                     argv[0], argv[0], argv[0]);
        return 1;
    }
    const std::string dir = argv[1];
    const std::string subcmd = argv[2];
    const std::string arg3 = argv[3];

    std::fprintf(stderr, "[inspect] loading graph.bin + posmap.bin + counts.bin + seqs.bin (slow) ...\n");
    Graph g;       load_graph(dir + "/graph.bin", g);
    std::fprintf(stderr, "[inspect]   graph: K=%lu E=%lu\n", (unsigned long)g.K, (unsigned long)g.edges.size());
    Posmap pm;     load_posmap(dir + "/posmap.bin", pm);
    std::fprintf(stderr, "[inspect]   posmap: chroms=%zu sigs=%zu\n", pm.chrom_names.size(), pm.sig_pos.size());
    std::vector<uint64_t> node_to_sig;
    build_node_to_sig(dir + "/counts.bin", g.K, node_to_sig);
    std::fprintf(stderr, "[inspect]   counts -> node_to_sig built\n");
    Seqs seqs;     load_seqs(dir + "/seqs.bin", seqs);
    std::fprintf(stderr, "[inspect]   seqs: K=%lu total=%lu\n", (unsigned long)seqs.K, (unsigned long)seqs.total_payload);

    std::vector<std::vector<int32_t>> adj(g.K);
    for (auto& e : g.edges) { adj[e.first].push_back((int32_t)e.second); adj[e.second].push_back((int32_t)e.first); }
    std::fprintf(stderr, "[inspect]   adjacency built\n");

    PosIndex pos_index;
    build_pos_index(g.color, node_to_sig, pm, pos_index);
    size_t total_pos_entries = 0;
    for (auto& kv : pos_index) total_pos_entries += kv.second.size();
    std::fprintf(stderr, "[inspect]   pos_index: %zu chroms, %zu (pos,node) entries (color=3 only)\n",
                 pos_index.size(), total_pos_entries);

    std::vector<Bubble> bubs;
    decompose_with_ids(adj, g.color, bubs);
    std::fprintf(stderr, "[inspect]   decompose: %zu bubbles\n", bubs.size());

    // Resolve target bubble(s)
    std::vector<uint32_t> targets;
    std::string key;
    std::vector<std::pair<std::string, std::string>> extra_kv;

    if (subcmd == "bubble") {
        uint32_t bid = (uint32_t)std::stoul(arg3);
        if (bid >= bubs.size()) die("bubble id out of range");
        targets.push_back(bid);
        key = "bubble_" + arg3;
    } else if (subcmd == "region") {
        auto colon = arg3.find(':');
        auto dash = arg3.find('-');
        if (colon == std::string::npos || dash == std::string::npos) die("region format chr:start-end");
        std::string chrom = arg3.substr(0, colon);
        uint32_t lo = (uint32_t)std::stoul(arg3.substr(colon + 1, dash - colon - 1));
        uint32_t hi = (uint32_t)std::stoul(arg3.substr(dash + 1));
        pick_bubbles_by_region(bubs, node_to_sig, pm, chrom, lo, hi, targets);
        key = "region_" + sanitize(arg3);
    } else if (subcmd == "sv") {
        if (argc < 5) die("sv requires <truth_vcf> <sv_id>");
        std::string vcf = arg3, sv_id = argv[4];
        std::vector<TruthRec> truth;
        load_truth(vcf, truth);
        const TruthRec* found = nullptr;
        for (const auto& r : truth) if (r.sv_id == sv_id) { found = &r; break; }
        if (!found) die("sv_id not in truth VCF");
        // Pick bubbles in [found.pos - 500, found.end + 500] on found.chrom
        uint32_t lo = found->pos > 500 ? found->pos - 500 : 0;
        uint32_t hi = found->end + 500;
        pick_bubbles_by_region(bubs, node_to_sig, pm, found->chrom, lo, hi, targets);
        key = "sv_" + sv_id;
        char buf[128];
        std::snprintf(buf, sizeof(buf), "{\"chrom\": \"%s\", \"pos\": %u, \"end\": %u, \"svtype\": \"%s\", \"svlen\": %d}",
                      found->chrom.c_str(), found->pos, found->end, found->svtype.c_str(), found->svlen);
        extra_kv.emplace_back("truth_sv", buf);
    } else die("unknown subcmd");

    if (targets.empty()) { std::fprintf(stderr, "[inspect] no bubbles matched\n"); return 2; }
    std::fprintf(stderr, "[inspect] targets: %zu bubble(s)\n", targets.size());

    mkdir_p(dir + "/inspect");
    std::string out_dir = dir + "/inspect/" + key;
    mkdir_p(out_dir);

    Subgraph raw;
    build_subgraph(bubs, targets, adj, g.color, g.ncov, g.tcov, node_to_sig, pm, seqs, pos_index, raw);
    std::fprintf(stderr, "[inspect] subgraph: %zu nodes, %zu edges\n", raw.nodes.size(), raw.edges.size());

    emit_gfa(raw, out_dir + "/raw.gfa", false);

    std::vector<std::string> warnings;
    Subgraph mA;
    run_mode_a(raw, mA, warnings);
    std::fprintf(stderr, "[inspect] mode A: %zu super-nodes, %zu edges\n", mA.nodes.size(), mA.edges.size());
    emit_gfa(mA, out_dir + "/modeA.gfa", false);

    Subgraph mB;
    run_mode_b(mA, mB, warnings);
    std::fprintf(stderr, "[inspect] mode B: %zu segments after refinement\n", mB.nodes.size());
    emit_gfa(mB, out_dir + "/modeB.gfa", true);

    emit_meta(raw, mA, mB, bubs, targets, pm, warnings, key, out_dir + "/meta.json", extra_kv);
    std::fprintf(stderr, "[inspect] wrote %s/{raw,modeA,modeB}.gfa + meta.json\n", out_dir.c_str());
    return 0;
}
