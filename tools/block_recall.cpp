// block_recall.cpp — per-super-node (Mode A "block") TP/FP precision + signature recall
// against the paired assembly oracle. Same truth set as tools/bubble_recall; differs only
// in the called-side aggregation: super-nodes from modea.bin instead of BFS bubbles.
//
// Inputs:  counts.bin  (post-prune sigs + bytes — NEW schema)
//          modea.bin   (super-nodes + chain_ids + blocks + runs + super_edges)
//          normal.fa   (germline ref)
//          tumor.fa    (truth tumor assembly)
//
// usage: block_recall <counts.bin> <modea.bin> <normal.fa> <tumor.fa>
//                     [--level 5] [--flank 10]

#include "svp.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

extern "C" {
#include "encoding.h"
}

static char comp(char c) { switch (c & 0xDF) { case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C'; } return 'N'; }

static void add_sigs(const char* fa, int level, int flank, std::unordered_set<uint64_t>& out) {
    std::string seq; std::vector<uint64_t> rs, re;
    read_fasta_flat(fa, flank, seq, rs, re);
    std::fprintf(stderr, "[block_recall] %s: %zu bp, %zu records; forward cores...\n", fa, seq.size(), rs.size());
    std::vector<Core> cs; std::vector<uint64_t> sig;
    extract_cores(seq, rs, re, level, cs); compute_sigs(seq, cs, flank, 0, rs, re, sig);
    for (uint64_t x : sig) out.insert(x);
    std::fprintf(stderr, "[block_recall] %s: reverse-complement cores...\n", fa);
    std::string rseq; std::vector<uint64_t> rrs, rre; const std::string bar(2 * flank + 1, 'N');
    for (size_t r = 0; r < rs.size(); r++) { if (!rrs.empty()) rseq += bar; rrs.push_back(rseq.size()); for (uint64_t p = re[r]; p > rs[r]; p--) rseq += comp(seq[p - 1]); rre.push_back(rseq.size()); }
    extract_cores(rseq, rrs, rre, level, cs); compute_sigs(rseq, cs, flank, 0, rrs, rre, sig);
    for (uint64_t x : sig) out.insert(x);
    std::fprintf(stderr, "[block_recall] %s: %zu distinct core sigs (both strands)\n", fa, out.size());
}

struct SuperNode { uint32_t chain_len, block_len, runs_len; uint8_t color; };

static FILE* xopenf(const char* p, const char* m) { FILE* f = std::fopen(p, m); if (!f) { std::fprintf(stderr, "[FATAL] open %s\n", p); std::abort(); } return f; }
static void xfread(void* p, size_t sz, size_t n, FILE* f) { if (std::fread(p, sz, n, f) != n) { std::fprintf(stderr, "[FATAL] short read\n"); std::abort(); } }

int main(int argc, char** argv) {
    if (argc < 5) { std::fprintf(stderr, "usage: block_recall <counts.bin> <modea.bin> <normal.fa> <tumor.fa> [--level 5] [--flank 10]\n"); return 2; }
    const char* cp = argv[1]; const char* mp = argv[2]; const char* nfa = argv[3]; const char* tfa = argv[4];
    int level = 5, flank = 10;
    for (int i = 5; i + 1 < argc; i += 2) {
        std::string k = argv[i]; const char* v = argv[i + 1];
        if (k == "--level") level = std::atoi(v);
        else if (k == "--flank") flank = std::atoi(v);
        else { std::fprintf(stderr, "[FATAL] arg %s\n", argv[i]); std::abort(); }
    }
    LCP_INIT();

    // 1. Truth somatic set S = T \ R (both strands).
    std::unordered_set<uint64_t> R, T;
    add_sigs(nfa, level, flank, R);
    add_sigs(tfa, level, flank, T);
    std::unordered_set<uint64_t> S; S.reserve(T.size());
    for (uint64_t x : T) if (!R.count(x)) S.insert(x);
    std::fprintf(stderr, "[block_recall] truth somatic |S| = %zu\n", S.size());

    // 2. counts.bin (new schema): u64 K | K * SigCount | u64 total | u16 seqlen[K] | u8 payload[total].
    //    Only the SigCount section is needed; skip the rest.
    FILE* fc = xopenf(cp, "rb");
    uint64_t K; xfread(&K, 8, 1, fc);
    std::vector<SigCount> sigs(K); xfread(sigs.data(), sizeof(SigCount), K, fc);
    std::fclose(fc);

    // Sig-level metric (matches tools/bubble_recall): in the new pipeline counts.bin
    // is already post-prune so every record is a "called" candidate — color=2 means n==0.
    uint64_t total_called = 0, called_hit_S = 0;
    for (uint64_t i = 0; i < K; i++) {
        if (sigs[i].n != 0) continue;
        total_called++;
        if (S.count(sigs[i].sig)) called_hit_S++;
    }

    // 3. modea.bin: u64 S | S*SuperNode | u64 tc | u32 chain_ids[tc] | (rest skipped)
    FILE* fm = xopenf(mp, "rb");
    uint64_t Smn; xfread(&Smn, 8, 1, fm);
    std::vector<SuperNode> snodes(Smn); xfread(snodes.data(), sizeof(SuperNode), Smn, fm);
    uint64_t tc; xfread(&tc, 8, 1, fm);
    std::vector<uint32_t> chain_ids(tc); xfread(chain_ids.data(), 4, tc, fm);
    std::fclose(fm);

    // 4. Per-color-2 super-node: TP iff any chain member's sig is in S.
    uint64_t sc1 = 0, sc2 = 0, sc3 = 0, sc2_tp = 0, sc1_tp = 0, sc3_tp = 0;
    uint64_t nodes_in_c2 = 0, nodes_in_c2_S = 0;
    size_t off = 0;
    for (uint64_t i = 0; i < Smn; i++) {
        const uint32_t m = snodes[i].chain_len; const uint8_t c = snodes[i].color;
        bool has_S = false;
        for (uint32_t j = 0; j < m; j++) {
            const uint64_t sig = sigs[chain_ids[off + j]].sig;
            if (S.count(sig)) has_S = true;
            if (c == 2) { nodes_in_c2++; if (S.count(sig)) nodes_in_c2_S++; }
        }
        off += m;
        if (c == 1) { sc1++; if (has_S) sc1_tp++; }
        else if (c == 2) { sc2++; if (has_S) sc2_tp++; }
        else { sc3++; if (has_S) sc3_tp++; }
    }

    // 5. Report.
    std::printf("== Signature-level (post-prune; called = n==0 in counts.bin) ==\n");
    std::printf("  |S| truth-somatic         = %zu\n", S.size());
    std::printf("  |C| called-somatic (n=0)  = %lu\n", (unsigned long)total_called);
    std::printf("  |C ∩ S|                   = %lu\n", (unsigned long)called_hit_S);
    std::printf("  RECALL    = %.4f%%\n", S.empty() ? 0.0 : 100.0 * (double)called_hit_S / (double)S.size());
    std::printf("  PRECISION = %.4f%%\n", total_called ? 100.0 * (double)called_hit_S / (double)total_called : 0.0);
    std::printf("\n== Block-level (Mode A super-nodes) ==\n");
    std::printf("  super-nodes total = %lu  (color=1 N-only=%lu  color=2 T-only=%lu  color=3 shared=%lu)\n",
                (unsigned long)Smn, (unsigned long)sc1, (unsigned long)sc2, (unsigned long)sc3);
    std::printf("  color=2 (tumor-only / somatic candidate):\n");
    std::printf("    TP-blocks (>=1 sig in S) = %lu / %lu  ->  PRECISION_block = %.4f%%\n",
                (unsigned long)sc2_tp, (unsigned long)sc2, sc2 ? 100.0 * (double)sc2_tp / (double)sc2 : 0.0);
    std::printf("    color=2 chain nodes      = %lu  (in S: %lu  -> RECALL_block = %.4f%%)\n",
                (unsigned long)nodes_in_c2, (unsigned long)nodes_in_c2_S,
                S.empty() ? 0.0 : 100.0 * (double)nodes_in_c2_S / (double)S.size());
    std::printf("  color=1 informational:     TP=%lu / %lu\n", (unsigned long)sc1_tp, (unsigned long)sc1);
    std::printf("  color=3 informational:     TP=%lu / %lu\n", (unsigned long)sc3_tp, (unsigned long)sc3);
    return 0;
}
