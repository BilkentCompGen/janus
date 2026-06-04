// truth_sig_localizer.cpp — for each truth somatic sig (∈ T \ R), find its super-node
// in modea.bin and classify the local bubble shape. Answers: do truth sigs land in
// std_both bubbles (canonical DEL/INS) or in T_only_no_N (allele isolated, no parallel
// N path)?
//
// usage: truth_sig_localizer <counts.bin> <modea.bin> <normal.fa> <tumor.fa> [--level 5] [--flank 10]

#include "svp.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
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
    std::fprintf(stderr, "[loc] %s: %zu bp, %zu records; fwd cores...\n", fa, seq.size(), rs.size());
    std::vector<Core> cs; std::vector<uint64_t> sig;
    extract_cores(seq, rs, re, level, cs); compute_sigs(seq, cs, flank, 0, rs, re, sig);
    for (uint64_t x : sig) out.insert(x);
    std::fprintf(stderr, "[loc] %s: RC cores...\n", fa);
    std::string rseq; std::vector<uint64_t> rrs, rre; const std::string bar(2 * flank + 1, 'N');
    for (size_t r = 0; r < rs.size(); r++) { if (!rrs.empty()) rseq += bar; rrs.push_back(rseq.size()); for (uint64_t p = re[r]; p > rs[r]; p--) rseq += comp(seq[p - 1]); rre.push_back(rseq.size()); }
    extract_cores(rseq, rrs, rre, level, cs); compute_sigs(rseq, cs, flank, 0, rrs, rre, sig);
    for (uint64_t x : sig) out.insert(x);
    std::fprintf(stderr, "[loc] %s: %zu distinct core sigs\n", fa, out.size());
}

struct SuperNode { uint32_t chain_len, block_len, runs_len; uint8_t color; };

int main(int argc, char** argv) {
    if (argc < 5) { std::fprintf(stderr, "usage: truth_sig_localizer <counts.bin> <modea.bin> <normal.fa> <tumor.fa> [--level L] [--flank F]\n"); return 2; }
    const char* cp = argv[1]; const char* mp = argv[2]; const char* nfa = argv[3]; const char* tfa = argv[4];
    int level = 5, flank = 10;
    for (int i = 5; i + 1 < argc; i += 2) {
        std::string k = argv[i]; const char* v = argv[i + 1];
        if (k == "--level") level = std::atoi(v);
        else if (k == "--flank") flank = std::atoi(v);
        else { std::fprintf(stderr, "[FATAL] arg %s\n", argv[i]); std::abort(); }
    }
    LCP_INIT();

    // Truth somatic set S = T \ R
    std::unordered_set<uint64_t> R, T;
    add_sigs(nfa, level, flank, R);
    add_sigs(tfa, level, flank, T);
    std::unordered_set<uint64_t> S; S.reserve(T.size());
    for (uint64_t x : T) if (!R.count(x)) S.insert(x);
    std::fprintf(stderr, "[loc] |S| truth-somatic = %zu\n", S.size());

    // counts.bin → sig→node_id
    FILE* fc = std::fopen(cp, "rb"); if (!fc) { std::fprintf(stderr,"[FATAL] %s\n",cp); std::abort(); }
    uint64_t K; std::fread(&K, 8, 1, fc);
    std::vector<SigCount> sigs(K); std::fread(sigs.data(), sizeof(SigCount), K, fc);
    std::fclose(fc);
    std::unordered_map<uint64_t, uint32_t> sig2node; sig2node.reserve(K);
    for (uint32_t i = 0; i < K; i++) sig2node.emplace(sigs[i].sig, i);

    // modea.bin → node→super, super_color
    FILE* fm = std::fopen(mp, "rb"); if (!fm) { std::fprintf(stderr,"[FATAL] %s\n",mp); std::abort(); }
    uint64_t Smn; std::fread(&Smn, 8, 1, fm);
    std::vector<SuperNode> snodes(Smn); std::fread(snodes.data(), sizeof(SuperNode), Smn, fm);
    uint64_t tc; std::fread(&tc, 8, 1, fm);
    std::vector<uint32_t> chain_ids(tc); std::fread(chain_ids.data(), 4, tc, fm);
    std::fclose(fm);
    std::vector<uint32_t> node2super(K, UINT32_MAX);
    size_t off = 0;
    for (uint32_t s = 0; s < Smn; s++) {
        for (uint32_t j = 0; j < snodes[s].chain_len; j++) node2super[chain_ids[off + j]] = s;
        off += snodes[s].chain_len;
    }

    // For each truth somatic sig, classify
    uint64_t in_graph = 0, not_in_graph = 0;
    uint64_t c1 = 0, c2 = 0, c3 = 0;
    uint64_t super_c1 = 0, super_c2 = 0, super_c3 = 0;
    for (uint64_t s : S) {
        auto it = sig2node.find(s);
        if (it == sig2node.end()) { not_in_graph++; continue; }
        in_graph++;
        const uint32_t nid = it->second;
        const uint32_t n = sigs[nid].n, t = sigs[nid].t;
        const uint8_t  c = (n && t) ? 3 : n ? 1 : 2;
        if (c == 1) c1++; else if (c == 2) c2++; else c3++;
        const uint32_t sid = node2super[nid];
        const uint8_t  sc = snodes[sid].color;
        if (sc == 1) super_c1++; else if (sc == 2) super_c2++; else super_c3++;
    }
    std::printf("== truth somatic sig localization ==\n");
    std::printf("  |S|                          = %zu\n", S.size());
    std::printf("  truth sigs FOUND in graph    = %lu  (%.2f%%)\n", (unsigned long)in_graph, 100.0 * in_graph / S.size());
    std::printf("  truth sigs NOT in graph      = %lu  (%.2f%% — recall miss)\n", (unsigned long)not_in_graph, 100.0 * not_in_graph / S.size());
    std::printf("\n  node-color distribution of truth-somatic sigs found in graph:\n");
    std::printf("    color=1 (normal-only)  = %lu  (%.2f%%)  ← SHOULD BE ZERO if signature pipeline is symmetric\n", (unsigned long)c1, 100.0 * c1 / in_graph);
    std::printf("    color=2 (tumor-only)   = %lu  (%.2f%%)  ← expected ~100%%\n", (unsigned long)c2, 100.0 * c2 / in_graph);
    std::printf("    color=3 (shared)       = %lu  (%.2f%%)  ← false-color-3 (would mask the SV)\n", (unsigned long)c3, 100.0 * c3 / in_graph);
    std::printf("\n  super-node-color (where the truth sig landed after Mode A):\n");
    std::printf("    super color=1 = %lu  (%.2f%%)\n", (unsigned long)super_c1, 100.0 * super_c1 / in_graph);
    std::printf("    super color=2 = %lu  (%.2f%%)\n", (unsigned long)super_c2, 100.0 * super_c2 / in_graph);
    std::printf("    super color=3 = %lu  (%.2f%%)\n", (unsigned long)super_c3, 100.0 * super_c3 / in_graph);
    return 0;
}
