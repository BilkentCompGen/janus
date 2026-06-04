// truth_sig_dump.cpp — dump |S| = (T \ R) as raw u64 records to a file.
//
// usage: truth_sig_dump <normal.fa> <tumor.fa> <out.bin> [--level 5] [--flank 10]

#include "svp.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_set>
#include <vector>

extern "C" {
#include "encoding.h"
}

static char comp(char c) { switch (c & 0xDF) { case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C'; } return 'N'; }

static void add_sigs(const char* fa, int level, int flank, std::unordered_set<uint64_t>& out) {
    std::string seq; std::vector<uint64_t> rs, re;
    read_fasta_flat(fa, flank, seq, rs, re);
    std::fprintf(stderr, "[dump] %s: %zu bp, %zu records; fwd cores...\n", fa, seq.size(), rs.size());
    std::vector<Core> cs; std::vector<uint64_t> sig;
    extract_cores(seq, rs, re, level, cs); compute_sigs(seq, cs, flank, 0, rs, re, sig);
    for (uint64_t x : sig) out.insert(x);
    std::fprintf(stderr, "[dump] %s: RC cores...\n", fa);
    std::string rseq; std::vector<uint64_t> rrs, rre; const std::string bar(2 * flank + 1, 'N');
    for (size_t r = 0; r < rs.size(); r++) { if (!rrs.empty()) rseq += bar; rrs.push_back(rseq.size()); for (uint64_t p = re[r]; p > rs[r]; p--) rseq += comp(seq[p - 1]); rre.push_back(rseq.size()); }
    extract_cores(rseq, rrs, rre, level, cs); compute_sigs(rseq, cs, flank, 0, rrs, rre, sig);
    for (uint64_t x : sig) out.insert(x);
    std::fprintf(stderr, "[dump] %s: %zu distinct sigs\n", fa, out.size());
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: truth_sig_dump <normal.fa> <tumor.fa> <out.bin> [--level L] [--flank F]\n"); return 2; }
    const char* nfa = argv[1]; const char* tfa = argv[2]; const char* op = argv[3];
    int level = 5, flank = 10;
    for (int i = 4; i + 1 < argc; i += 2) {
        std::string k = argv[i]; const char* v = argv[i + 1];
        if (k == "--level") level = std::atoi(v);
        else if (k == "--flank") flank = std::atoi(v);
        else { std::fprintf(stderr, "[FATAL] arg %s\n", argv[i]); std::abort(); }
    }
    LCP_INIT();
    std::unordered_set<uint64_t> R, T;
    add_sigs(nfa, level, flank, R);
    add_sigs(tfa, level, flank, T);
    std::vector<uint64_t> S; S.reserve(T.size());
    for (uint64_t x : T) if (!R.count(x)) S.push_back(x);
    std::fprintf(stderr, "[dump] |S| = %zu  -> %s\n", S.size(), op);
    FILE* f = std::fopen(op, "wb"); if (!f) { std::fprintf(stderr, "[FATAL] open %s\n", op); std::abort(); }
    if (std::fwrite(S.data(), 8, S.size(), f) != S.size()) { std::fprintf(stderr, "[FATAL] write\n"); std::abort(); }
    std::fclose(f);
    return 0;
}
