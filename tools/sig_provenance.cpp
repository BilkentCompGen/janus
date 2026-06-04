// sig_provenance.cpp — for every sig in counts.bin, record how many times it appears
// as an LCP core sig in (normal.fa fwd, normal.fa RC, tumor.fa fwd, tumor.fa RC).
// Output schema: u64 K | K × Prov{u32 n_fwd, u32 n_rc, u32 t_fwd, u32 t_rc}.
//
// usage: sig_provenance <counts.bin> <normal.fa> <tumor.fa> <out.bin> [--level 5] [--flank 10]

#include "svp.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <string>
#include <unordered_map>
#include <vector>

extern "C" {
#include "encoding.h"
}

static char comp(char c) {
    switch (c & 0xDF) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
    }
    return 'N';
}

struct Prov { uint32_t n_fwd, n_rc, t_fwd, t_rc; };

static void build_rc(const std::string& seq,
                     const std::vector<uint64_t>& rs,
                     const std::vector<uint64_t>& re,
                     int flank,
                     std::string& rseq,
                     std::vector<uint64_t>& rrs,
                     std::vector<uint64_t>& rre) {
    const std::string bar(2 * flank + 1, 'N');
    for (size_t r = 0; r < rs.size(); r++) {
        if (!rrs.empty()) rseq += bar;
        rrs.push_back(rseq.size());
        for (uint64_t p = re[r]; p > rs[r]; p--) rseq += comp(seq[p - 1]);
        rre.push_back(rseq.size());
    }
}

static void count_pass(const std::string& seq,
                       const std::vector<uint64_t>& rs,
                       const std::vector<uint64_t>& re,
                       int level, int flank,
                       std::unordered_map<uint64_t, Prov>& prov,
                       uint32_t Prov::*field) {
    std::vector<Core> cs; std::vector<uint64_t> sig;
    extract_cores(seq, rs, re, level, cs);
    compute_sigs(seq, cs, flank, 0, rs, re, sig);
    uint64_t hits = 0;
    for (uint64_t x : sig) {
        auto it = prov.find(x);
        if (it == prov.end()) continue;
        (it->second.*field)++;
        hits++;
    }
    std::fprintf(stderr, "[prov]   %zu cores extracted, %lu hit called set\n",
                 sig.size(), (unsigned long)hits);
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: sig_provenance <counts.bin> <normal.fa> <tumor.fa> <out.bin> [--level L] [--flank F]\n");
        return 2;
    }
    const char* cp = argv[1]; const char* nfa = argv[2]; const char* tfa = argv[3]; const char* op = argv[4];
    int level = 5, flank = 10;
    for (int i = 5; i + 1 < argc; i += 2) {
        std::string k = argv[i]; const char* v = argv[i + 1];
        if (k == "--level") level = std::atoi(v);
        else if (k == "--flank") flank = std::atoi(v);
        else { std::fprintf(stderr, "[FATAL] arg %s\n", argv[i]); std::abort(); }
    }
    LCP_INIT();

    FILE* fc = std::fopen(cp, "rb");
    if (!fc) { std::fprintf(stderr, "[FATAL] open %s\n", cp); std::abort(); }
    uint64_t K; std::fread(&K, 8, 1, fc);
    std::vector<SigCount> sigs(K);
    std::fread(sigs.data(), sizeof(SigCount), K, fc);
    std::fclose(fc);
    std::fprintf(stderr, "[prov] loaded counts.bin: K=%lu sigs\n", (unsigned long)K);

    std::unordered_map<uint64_t, Prov> prov;
    prov.reserve(K * 2);
    for (uint64_t i = 0; i < K; i++) prov.emplace(sigs[i].sig, Prov{0, 0, 0, 0});
    std::fprintf(stderr, "[prov] provenance map initialized (%zu entries)\n", prov.size());

    {
        std::string seq; std::vector<uint64_t> rs, re;
        std::fprintf(stderr, "[prov] reading normal.fa (%s)...\n", nfa);
        read_fasta_flat(nfa, flank, seq, rs, re);
        std::fprintf(stderr, "[prov]   %zu bp, %zu records; FWD pass...\n", seq.size(), rs.size());
        count_pass(seq, rs, re, level, flank, prov, &Prov::n_fwd);
        std::fprintf(stderr, "[prov]   building RC...\n");
        std::string rseq; std::vector<uint64_t> rrs, rre;
        build_rc(seq, rs, re, flank, rseq, rrs, rre);
        std::fprintf(stderr, "[prov]   RC pass...\n");
        count_pass(rseq, rrs, rre, level, flank, prov, &Prov::n_rc);
    }
    {
        std::string seq; std::vector<uint64_t> rs, re;
        std::fprintf(stderr, "[prov] reading tumor.fa (%s)...\n", tfa);
        read_fasta_flat(tfa, flank, seq, rs, re);
        std::fprintf(stderr, "[prov]   %zu bp, %zu records; FWD pass...\n", seq.size(), rs.size());
        count_pass(seq, rs, re, level, flank, prov, &Prov::t_fwd);
        std::fprintf(stderr, "[prov]   building RC...\n");
        std::string rseq; std::vector<uint64_t> rrs, rre;
        build_rc(seq, rs, re, flank, rseq, rrs, rre);
        std::fprintf(stderr, "[prov]   RC pass...\n");
        count_pass(rseq, rrs, rre, level, flank, prov, &Prov::t_rc);
    }

    // Dump in counts.bin order.
    std::vector<Prov> arr(K);
    for (uint64_t i = 0; i < K; i++) {
        auto it = prov.find(sigs[i].sig);
        arr[i] = (it != prov.end()) ? it->second : Prov{0, 0, 0, 0};
    }
    FILE* fo = std::fopen(op, "wb");
    if (!fo) { std::fprintf(stderr, "[FATAL] open %s\n", op); std::abort(); }
    std::fwrite(&K, 8, 1, fo);
    std::fwrite(arr.data(), sizeof(Prov), K, fo);
    std::fclose(fo);
    std::fprintf(stderr, "[prov] dumped %lu records -> %s\n", (unsigned long)K, op);

    // Summary classification of called sigs.
    uint64_t both = 0, n_only = 0, t_only = 0, phantom = 0;
    uint64_t any_fwd = 0, any_rc = 0, both_strand = 0;
    uint64_t rep_n_fwd = 0, rep_n_rc = 0, rep_t_fwd = 0, rep_t_rc = 0;
    uint64_t max_rep = 0;
    for (uint64_t i = 0; i < K; i++) {
        const auto& p = arr[i];
        const bool n = (p.n_fwd + p.n_rc) > 0;
        const bool t = (p.t_fwd + p.t_rc) > 0;
        const bool f = (p.n_fwd + p.t_fwd) > 0;
        const bool r = (p.n_rc + p.t_rc) > 0;
        if (n && t) both++;
        else if (n) n_only++;
        else if (t) t_only++;
        else phantom++;
        if (f) any_fwd++;
        if (r) any_rc++;
        if (f && r) both_strand++;
        if (p.n_fwd > 1) rep_n_fwd++;
        if (p.n_rc  > 1) rep_n_rc++;
        if (p.t_fwd > 1) rep_t_fwd++;
        if (p.t_rc  > 1) rep_t_rc++;
        const uint32_t mx = std::max({p.n_fwd, p.n_rc, p.t_fwd, p.t_rc});
        if (mx > max_rep) max_rep = mx;
    }
    std::fprintf(stderr, "\n[prov] called sigs (K=%lu) classification by sample presence:\n", (unsigned long)K);
    std::fprintf(stderr, "  BOTH (n & t)        = %lu (%.2f%%)\n", (unsigned long)both,    100.0 * both    / K);
    std::fprintf(stderr, "  N_only              = %lu (%.2f%%)\n", (unsigned long)n_only,  100.0 * n_only  / K);
    std::fprintf(stderr, "  T_only              = %lu (%.2f%%)\n", (unsigned long)t_only,  100.0 * t_only  / K);
    std::fprintf(stderr, "  TRUE PHANTOM        = %lu (%.2f%%)\n", (unsigned long)phantom, 100.0 * phantom / K);
    std::fprintf(stderr, "\n[prov] strand presence (across all samples):\n");
    std::fprintf(stderr, "  any FWD             = %lu (%.2f%%)\n", (unsigned long)any_fwd,    100.0 * any_fwd    / K);
    std::fprintf(stderr, "  any RC              = %lu (%.2f%%)\n", (unsigned long)any_rc,     100.0 * any_rc     / K);
    std::fprintf(stderr, "  both FWD and RC     = %lu (%.2f%%)\n", (unsigned long)both_strand, 100.0 * both_strand / K);
    std::fprintf(stderr, "\n[prov] repetition (sig appears >1 time in one (sample, strand)):\n");
    std::fprintf(stderr, "  n_fwd > 1           = %lu (%.2f%%)\n", (unsigned long)rep_n_fwd, 100.0 * rep_n_fwd / K);
    std::fprintf(stderr, "  n_rc  > 1           = %lu (%.2f%%)\n", (unsigned long)rep_n_rc,  100.0 * rep_n_rc  / K);
    std::fprintf(stderr, "  t_fwd > 1           = %lu (%.2f%%)\n", (unsigned long)rep_t_fwd, 100.0 * rep_t_fwd / K);
    std::fprintf(stderr, "  t_rc  > 1           = %lu (%.2f%%)\n", (unsigned long)rep_t_rc,  100.0 * rep_t_rc  / K);
    std::fprintf(stderr, "  max repetition (any single bucket) = %u\n", (unsigned)max_rep);
    return 0;
}
