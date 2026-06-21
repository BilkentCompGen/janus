// Repeat resolution (X-cut) over the compacted super-graph. The fill pass recorded, per coverage-
// flagged repeat run, the in->run->out spans reads threaded through it (dir/threads.bin). Here we
// split a super-node into one copy per (in,out) pair ONLY when the reads pair its in-flanks to its
// out-flanks 1:1 AND that pairing covers the node's ACTUAL super-graph topology. The reads carry
// the pairing, so an observed A1->C->B2 stays A1->C1->B2; we never invent A1->C->B1.
//
// A node is left collapsed (and logged to resolve_ambiguous.tsv) whenever the topology shows in>=2
// AND out>=2 but the threading does not pin it down: a missing pairing (repeat longer than a read,
// so some in/out was never linked), a one-to-many fan, or a cross-product tangle.
#include "graph.h"
#include "svp.h"
#include "unordered_dense.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>

static const uint32_t RNONE = UINT32_MAX;

static bool read_threads(const std::string& path, std::vector<RepTriple>& out) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    uint64_t n;
    if (std::fread(&n, 8, 1, f) != 1) { std::fprintf(stderr, "[FATAL] threads.bin header\n"); std::abort(); }
    out.resize(n);
    if (n && std::fread(out.data(), sizeof(RepTriple), n, f) != n) { std::fprintf(stderr, "[FATAL] threads.bin body\n"); std::abort(); }
    std::fclose(f);
    return true;
}

// Does the read-threading pin C down? The full (in,out) pairs must be a 1:1 matching whose ins are
// exactly C's topological in-set and whose outs are exactly its out-set. Returns the matched pairs.
static const char* match_topology(const ankerl::unordered_dense::map<uint64_t, uint32_t>& pairs,
                                  const std::vector<uint32_t>& in_set, const std::vector<uint32_t>& out_set,
                                  std::vector<std::pair<uint32_t, uint32_t>>& matched) {
    matched.clear();
    ankerl::unordered_dense::map<uint32_t, uint32_t> indeg, outdeg;
    for (const auto& kv : pairs) {
        const uint32_t a = (uint32_t)(kv.first >> 32), b = (uint32_t)kv.first;
        if (a == RNONE || b == RNONE) continue;                 // unspanned read -> incomplete
        indeg[a]++; outdeg[b]++;
        matched.push_back({a, b});
    }
    if (matched.size() != in_set.size() || in_set.size() != out_set.size()) return "partial/asymmetric";
    for (uint32_t a : in_set)  if (indeg[a]  != 1) return "uncovered_in";    // an in-flank with !=1 read pairing
    for (uint32_t b : out_set) if (outdeg[b] != 1) return "uncovered_out";   // an out-flank with !=1 read pairing
    return nullptr;                                              // perfect, topology-covering matching
}

void resolve_repeats(const std::string& dir, std::vector<SuperNode>& snodes, std::vector<uint8_t>& blocks,
                     std::vector<CovRun>& runs, std::vector<IEdge>& super_edges,
                     const std::vector<uint32_t>& node_to_super) {
    std::vector<RepTriple> raw;
    if (!read_threads(dir + "/threads.bin", raw)) return;       // capture not enabled -> no-op

    // observed super (in,out) pairs per super mid-node
    ankerl::unordered_dense::map<uint32_t, ankerl::unordered_dense::map<uint64_t, uint32_t>> bymid;
    for (const RepTriple& t : raw) {
        const uint32_t smid = node_to_super[t.mid];
        const uint32_t sin  = t.in  == RNONE ? RNONE : node_to_super[t.in];
        const uint32_t sout = t.out == RNONE ? RNONE : node_to_super[t.out];
        if (sin == smid || sout == smid) continue;              // edge internal to the unitig
        bymid[smid][((uint64_t)sin << 32) | sout] += t.count;
    }

    // super-graph in/out adjacency (super_edges is already unique, so each entry is a distinct flank)
    const uint32_t S = (uint32_t)snodes.size();
    ankerl::unordered_dense::map<uint32_t, std::vector<uint32_t>> in_set, out_set;
    for (const IEdge& e : super_edges) { out_set[e.first].push_back(e.second); in_set[e.second].push_back(e.first); }

    std::vector<uint64_t> boff(S + 1, 0), roff(S + 1, 0);        // block/runs offsets of the ORIGINAL nodes
    for (uint32_t i = 0; i < S; i++) { boff[i + 1] = boff[i] + snodes[i].block_len; roff[i + 1] = roff[i] + snodes[i].runs_len; }

    FILE* fa = std::fopen((dir + "/resolve_ambiguous.tsv").c_str(), "w");
    std::fprintf(fa, "super_mid\tn_in\tn_out\tin_supers\tout_supers\tthreaded_pairs\treason\n");
    ankerl::unordered_dense::set<uint64_t> rm;                  // edges to drop: (from<<32)|to
    std::vector<IEdge> add;
    uint64_t resolved = 0, copies = 0, ambiguous = 0, amb_dumped = 0;
    std::vector<std::pair<uint32_t, uint32_t>> matched;

    for (const auto& kv : bymid) {
        const uint32_t C = kv.first;
        static const std::vector<uint32_t> EMPTY;
        const auto iit = in_set.find(C), oit = out_set.find(C);
        const std::vector<uint32_t>& is = iit == in_set.end() ? EMPTY : iit->second;
        const std::vector<uint32_t>& os = oit == out_set.end() ? EMPTY : oit->second;
        if (is.size() < 2 || os.size() < 2) continue;           // not a branching repeat in the graph

        const char* fail = match_topology(kv.second, is, os, matched);
        if (!fail) {                                            // reads pin it -> SPLIT
            resolved++;
            for (const auto& p : matched) {
                const uint32_t ci = (uint32_t)snodes.size();
                const SuperNode sc = snodes[C];                 // copy carries C's block/runs/colour
                std::vector<uint8_t> bcopy(blocks.begin() + boff[C], blocks.begin() + boff[C] + sc.block_len);
                std::vector<CovRun>  rcopy(runs.begin() + roff[C], runs.begin() + roff[C] + sc.runs_len);
                blocks.insert(blocks.end(), bcopy.begin(), bcopy.end());
                runs.insert(runs.end(), rcopy.begin(), rcopy.end());
                snodes.push_back(sc);
                add.push_back({p.first, ci}); add.push_back({ci, p.second});
                copies++;
            }
            for (uint32_t a : is) rm.insert(((uint64_t)a << 32) | C);
            for (uint32_t b : os) rm.insert(((uint64_t)C << 32) | b);
        } else {                                                // topology branches but reads don't pin it
            ambiguous++;
            if (amb_dumped++ < 200000) {
                std::string si, so, sp;
                for (uint32_t a : is) si += std::to_string(a) + ",";
                for (uint32_t b : os) so += std::to_string(b) + ",";
                for (const auto& pr : kv.second) { const uint32_t a = (uint32_t)(pr.first >> 32), b = (uint32_t)pr.first;
                    sp += (a == RNONE ? "NONE" : std::to_string(a)) + ">" + (b == RNONE ? "NONE" : std::to_string(b)) + ":" + std::to_string(pr.second) + " "; }
                std::fprintf(fa, "%u\t%zu\t%zu\t%s\t%s\t%s\t%s\n", C, is.size(), os.size(), si.c_str(), so.c_str(), sp.c_str(), fail);
            }
        }
    }
    std::fclose(fa);

    // apply the edge surgery: drop C's old flank edges, add the per-copy edges, re-sort/unique
    std::vector<IEdge> ne; ne.reserve(super_edges.size() + add.size());
    for (const IEdge& e : super_edges) if (!rm.count(((uint64_t)e.first << 32) | e.second)) ne.push_back(e);
    ne.insert(ne.end(), add.begin(), add.end());
    std::sort(ne.begin(), ne.end()); ne.erase(std::unique(ne.begin(), ne.end()), ne.end());
    super_edges.swap(ne);

    std::fprintf(stderr,
        "[resolve] X-cut: split %lu repeat super-nodes -> %lu copies (S %u -> %zu) | left collapsed (2-in/2-out, threading incomplete) %lu -> resolve_ambiguous.tsv\n",
        (unsigned long)resolved, (unsigned long)copies, S, snodes.size(), (unsigned long)ambiguous);
}
