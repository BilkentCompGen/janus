// bubble_recall.cpp — per-bubble TP/FP precision + signature recall for the read-based
// LCP-Cortex caller against a paired assembly oracle (rule 1: assemblies are GT only).
//
// Inputs:  counts.bin  (svp lcp output: sig -> n,t)
//          graph.bin   (svp graph output: K, color, ncov, tcov, E, edges)
//          normal.fa   (germline reference assembly)
//          tumor.fa    (truth tumor assembly, e.g. simulation/truth/wgs/out.fa)
//
// Output:
//   Signature-level (same as tools/recall.cpp):
//     truth somatic set S = (tumor cores ∪ tumor-RC cores) \ (normal cores ∪ normal-RC cores)
//     called somatic set C = { sig : counts.n==0 && counts.t>=min_mult }
//     RECALL = |C ∩ S| / |S|,  PRECISION = |C ∩ S| / |C|
//   Bubble-level:
//     For each connected component of non-shared nodes in the graph, classify as
//     I (tumor-only), D (normal-only), or C (mixed) by the colour bits seen.
//     A bubble is TP if ANY of its node sigs is in S; FP otherwise.
//     Precision_bubble = TP_I / (TP_I + FP_I).
//
// usage: bubble_recall <counts.bin> <graph.bin> <normal.fa> <tumor.fa>
//                      [--level 5] [--flank 10] [--min-mult 3]

#include "svp.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

extern "C" {
#include "encoding.h"
}

static char comp(char c) { switch (c & 0xDF) { case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C'; } return 'N'; }

// Both-strand LCP-core signature set of a FASTA (reads are unstranded; a read core can
// match either the forward or the reverse-complement core of the assembly).
static void add_sigs(const char* fa, int level, int flank, std::unordered_set<uint64_t>& out) {
    std::string seq; std::vector<uint64_t> rs, re;
    read_fasta_flat(fa, flank, seq, rs, re);
    std::fprintf(stderr, "[bubble_recall] %s: %zu bp, %zu records; forward cores...\n", fa, seq.size(), rs.size());
    std::vector<Core> cs; std::vector<uint64_t> sig;
    extract_cores(seq, rs, re, level, cs); compute_sigs(seq, cs, flank, 0, rs, re, sig);
    for (uint64_t x : sig) out.insert(x);
    std::fprintf(stderr, "[bubble_recall] %s: reverse-complement cores...\n", fa);
    std::string rseq; std::vector<uint64_t> rrs, rre; const std::string bar(2 * flank + 1, 'N');
    for (size_t r = 0; r < rs.size(); r++) { if (!rrs.empty()) rseq += bar; rrs.push_back(rseq.size()); for (uint64_t p = re[r]; p > rs[r]; p--) rseq += comp(seq[p - 1]); rre.push_back(rseq.size()); }
    extract_cores(rseq, rrs, rre, level, cs); compute_sigs(rseq, cs, flank, 0, rrs, rre, sig);
    for (uint64_t x : sig) out.insert(x);
    std::fprintf(stderr, "[bubble_recall] %s: %zu distinct core sigs (both strands)\n", fa, out.size());
}

using IEdge = std::pair<uint32_t, uint32_t>;

static FILE* xopen(const char* p, const char* m) { FILE* f = std::fopen(p, m); if (!f) { std::fprintf(stderr, "[FATAL] open %s\n", p); std::abort(); } return f; }
static void xfread(void* p, size_t sz, size_t n, FILE* f) { if (std::fread(p, sz, n, f) != n) { std::fprintf(stderr, "[FATAL] short read\n"); std::abort(); } }

int main(int argc, char** argv) {
    if (argc < 5) { std::fprintf(stderr, "usage: bubble_recall <counts.bin> <graph.bin> <normal.fa> <tumor.fa> [--level 5] [--flank 10] [--min-mult 3]\n"); return 2; }
    const char* cp = argv[1]; const char* gp = argv[2]; const char* nfa = argv[3]; const char* tfa = argv[4];
    int level = 5, flank = 10, min_mult = 3;
    for (int i = 5; i + 1 < argc; i += 2) {
        std::string k = argv[i]; const char* v = argv[i + 1];
        if (k == "--level") level = std::atoi(v);
        else if (k == "--flank") flank = std::atoi(v);
        else if (k == "--min-mult") min_mult = std::atoi(v);
        else { std::fprintf(stderr, "[FATAL] arg %s\n", argv[i]); std::abort(); }
    }
    LCP_INIT();

    // 1. Truth somatic sig set S = T \ R (both strands)
    std::unordered_set<uint64_t> R, T;
    add_sigs(nfa, level, flank, R);
    add_sigs(tfa, level, flank, T);
    std::unordered_set<uint64_t> S; S.reserve(T.size());
    for (uint64_t x : T) if (!R.count(x)) S.insert(x);
    std::fprintf(stderr, "[bubble_recall] truth somatic |S| = %zu\n", S.size());

    // 2. Walk counts.bin in original order, apply the same prune as run_graph, recover
    //    node-id -> sig in graph.bin order (i'th passing sig = node id i).
    FILE* fc = xopen(cp, "rb");
    uint64_t N; xfread(&N, 8, 1, fc);
    std::vector<uint64_t> node_sigs; node_sigs.reserve(N / 2);
    uint64_t total_called = 0, called_hit_S = 0;
    SigCount c;
    for (uint64_t i = 0; i < N; i++) {
        xfread(&c, sizeof(c), 1, fc);
        if ((c.n > c.t ? c.n : c.t) < (uint32_t)min_mult) continue;
        node_sigs.push_back(c.sig);
        if (c.n == 0 && c.t >= (uint32_t)min_mult) { total_called++; if (S.count(c.sig)) called_hit_S++; }
    }
    std::fclose(fc);

    // 3. Load graph.bin and confirm node count matches.
    FILE* fg = xopen(gp, "rb");
    uint64_t K; xfread(&K, 8, 1, fg);
    if (K != node_sigs.size()) { std::fprintf(stderr, "[FATAL] graph.bin K=%lu but counts.bin filtered=%zu (min_mult mismatch with the graph stage that wrote graph.bin?)\n", (unsigned long)K, node_sigs.size()); std::abort(); }
    std::vector<uint8_t> color(K); xfread(color.data(), 1, K, fg);
    std::vector<uint32_t> ncov(K), tcov(K); xfread(ncov.data(), 4, K, fg); xfread(tcov.data(), 4, K, fg);
    uint64_t E; xfread(&E, 8, 1, fg);
    std::vector<IEdge> ie(E); xfread(ie.data(), sizeof(IEdge), E, fg);
    std::fclose(fg);

    // 4. Adjacency.
    std::vector<std::vector<int32_t>> adj(K);
    for (const IEdge& e : ie) { adj[e.first].push_back((int32_t)e.second); adj[e.second].push_back((int32_t)e.first); }

    // 5. BFS connected components of colour != 3; per-bubble TP/FP.
    std::vector<char> vis(K, 0);
    uint64_t b_tot = 0, b_I = 0, b_D = 0, b_C = 0, b_I_tp = 0, b_D_tp = 0, b_C_tp = 0;
    uint64_t n_I_nodes = 0, n_I_nodes_in_S = 0;
    std::queue<uint32_t> q;
    for (uint32_t s = 0; s < K; s++) {
        if (color[s] == 3 || vis[s]) continue;
        q.push(s); vis[s] = 1;
        uint8_t seen = 0; uint32_t nn = 0; bool has_S = false;
        while (!q.empty()) {
            const uint32_t u = q.front(); q.pop();
            seen |= color[u]; nn++;
            if (S.count(node_sigs[u])) has_S = true;
            for (int32_t v : adj[u]) if (!vis[v] && color[v] != 3) { vis[v] = 1; q.push((uint32_t)v); }
        }
        b_tot++;
        if (seen == 2) { b_I++; if (has_S) b_I_tp++; n_I_nodes += nn; }
        else if (seen == 1) { b_D++; if (has_S) b_D_tp++; }
        else { b_C++; if (has_S) b_C_tp++; }
    }
    // I-type per-node share of S — informational (sigs in I-bubbles ∩ S).
    for (uint32_t i = 0; i < K; i++) if (color[i] == 2 && S.count(node_sigs[i])) n_I_nodes_in_S++;

    // 6. Report.
    std::printf("== Signature-level (same metric as tools/recall) ==\n");
    std::printf("  min_mult=%d  |S| truth-somatic=%zu\n", min_mult, S.size());
    std::printf("  |C| called-somatic (n=0,t>=min_mult) = %lu\n", (unsigned long)total_called);
    std::printf("  |C ∩ S| = %lu\n", (unsigned long)called_hit_S);
    std::printf("  RECALL    = %.4f%%\n", S.empty() ? 0.0 : 100.0 * (double)called_hit_S / (double)S.size());
    std::printf("  PRECISION = %.4f%%\n", total_called ? 100.0 * (double)called_hit_S / (double)total_called : 0.0);
    std::printf("\n== Bubble-level ==\n");
    std::printf("  bubbles total = %lu  (I=%lu  D=%lu  C=%lu)\n", (unsigned long)b_tot, (unsigned long)b_I, (unsigned long)b_D, (unsigned long)b_C);
    std::printf("  I (tumor-only / somatic candidate):\n");
    std::printf("    TP-bubbles (>=1 sig in S) = %lu / %lu  ->  PRECISION_bubble = %.4f%%\n",
                (unsigned long)b_I_tp, (unsigned long)b_I, b_I ? 100.0 * (double)b_I_tp / (double)b_I : 0.0);
    std::printf("    nodes in I-bubbles = %lu  (in S: %lu  -> %.4f%% of I-bubble nodes are real somatic)\n",
                (unsigned long)n_I_nodes, (unsigned long)n_I_nodes_in_S,
                n_I_nodes ? 100.0 * (double)n_I_nodes_in_S / (double)n_I_nodes : 0.0);
    std::printf("  D (normal-only):   TP=%lu / %lu (informational; not somatic)\n", (unsigned long)b_D_tp, (unsigned long)b_D);
    std::printf("  C (mixed):         TP=%lu / %lu\n", (unsigned long)b_C_tp, (unsigned long)b_C);
    return 0;
}
