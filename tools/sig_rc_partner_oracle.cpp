// sig_rc_partner_oracle.cpp — POSITIONAL pairing oracle.
//
// For each fwd LCP core in normal.fa (or tumor.fa), compute the expected RC-core
// position in the RC sequence: a core at seq[st, en) → rseq[rrs[r] + re[r] - en,
//                                                          rrs[r] + re[r] - st).
// If an RC core exists at exactly that position, the pair (fwd_sig, rc_sig) is
// recorded as a TRUE PARTNER.
//
// Then for each sig in counts.bin:
//   - Look up via fwd_sig → rc_sig table (and vice versa)
//   - Check if the partner sig is also in counts.bin
//   - Report pairing rate stratified by provenance class
//
// usage: sig_rc_partner_oracle <counts.bin> <provenance.bin> <normal.fa> <tumor.fa> [--level 5] [--flank 10]

#include "svp.h"
#include <algorithm>
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

static char comp(char c) {
    switch (c & 0xDF) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
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

static void pair_sample(const char* fa, int level, int flank,
                        std::vector<std::pair<uint64_t,uint64_t>>& pairs_out,
                        uint64_t& exact_pairs, uint64_t& no_match) {
    std::string seq; std::vector<uint64_t> rs, re;
    read_fasta_flat(fa, flank, seq, rs, re);
    std::fprintf(stderr, "[oracle] %s: %zu bp, %zu records\n", fa, seq.size(), rs.size());

    std::vector<Core> fwd_cs; std::vector<uint64_t> fwd_sig;
    extract_cores(seq, rs, re, level, fwd_cs);
    compute_sigs(seq, fwd_cs, flank, 0, rs, re, fwd_sig);
    std::fprintf(stderr, "[oracle]   fwd cores: %zu\n", fwd_cs.size());

    std::string rseq; std::vector<uint64_t> rrs, rre;
    build_rc(seq, rs, re, flank, rseq, rrs, rre);
    std::vector<Core> rc_cs; std::vector<uint64_t> rc_sig;
    extract_cores(rseq, rrs, rre, level, rc_cs);
    compute_sigs(rseq, rc_cs, flank, 0, rrs, rre, rc_sig);
    std::fprintf(stderr, "[oracle]   rc cores: %zu\n", rc_cs.size());

    size_t fi = 0, ri = 0;
    for (size_t r = 0; r < rs.size(); r++) {
        while (fi < fwd_cs.size() && fwd_cs[fi].start < rs[r]) fi++;
        size_t f_lo = fi;
        while (fi < fwd_cs.size() && fwd_cs[fi].end <= re[r]) fi++;
        size_t f_hi = fi;
        while (ri < rc_cs.size() && rc_cs[ri].start < rrs[r]) ri++;
        size_t r_lo = ri;
        while (ri < rc_cs.size() && rc_cs[ri].end <= rre[r]) ri++;
        size_t r_hi = ri;
        // Build position-keyed map: (start, end) → rc_idx for RC cores in this chromosome.
        std::unordered_map<uint64_t, size_t> rc_key_to_idx;
        rc_key_to_idx.reserve(r_hi - r_lo);
        for (size_t k = r_lo; k < r_hi; k++) {
            const uint64_t key = (rc_cs[k].start << 32) ^ rc_cs[k].end;
            rc_key_to_idx[key] = k;
        }
        for (size_t k = f_lo; k < f_hi; k++) {
            const uint64_t st = fwd_cs[k].start, en = fwd_cs[k].end;
            const uint64_t expected_rc_start = rrs[r] + (re[r] - en);
            const uint64_t expected_rc_end   = rrs[r] + (re[r] - st);
            const uint64_t key = (expected_rc_start << 32) ^ expected_rc_end;
            auto it = rc_key_to_idx.find(key);
            if (it != rc_key_to_idx.end()) {
                pairs_out.emplace_back(fwd_sig[k], rc_sig[it->second]);
                exact_pairs++;
            } else {
                no_match++;
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: sig_rc_partner_oracle <counts.bin> <provenance.bin> <normal.fa> <tumor.fa> [--level L] [--flank F]\n");
        return 2;
    }
    const char* cp = argv[1]; const char* pp = argv[2]; const char* nfa = argv[3]; const char* tfa = argv[4];
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
    std::vector<SigCount> sigs(K); std::fread(sigs.data(), sizeof(SigCount), K, fc);
    std::fclose(fc);
    std::unordered_set<uint64_t> sigset; sigset.reserve(K * 2);
    std::unordered_map<uint64_t, uint32_t> sig_to_idx; sig_to_idx.reserve(K * 2);
    for (uint32_t i = 0; i < K; i++) {
        sigset.insert(sigs[i].sig);
        sig_to_idx.emplace(sigs[i].sig, i);
    }
    std::fprintf(stderr, "[oracle] counts.bin K = %lu\n", (unsigned long)K);

    FILE* fp = std::fopen(pp, "rb");
    if (!fp) { std::fprintf(stderr, "[FATAL] open %s\n", pp); std::abort(); }
    uint64_t Kp; std::fread(&Kp, 8, 1, fp);
    std::vector<Prov> prov(Kp); std::fread(prov.data(), sizeof(Prov), Kp, fp);
    std::fclose(fp);
    std::fprintf(stderr, "[oracle] provenance Kp = %lu\n", (unsigned long)Kp);

    // Collect (fwd_sig, rc_sig) pairs from both assemblies.
    std::vector<std::pair<uint64_t, uint64_t>> all_pairs;
    uint64_t n_exact_n = 0, n_no_n = 0, n_exact_t = 0, n_no_t = 0;
    pair_sample(nfa, level, flank, all_pairs, n_exact_n, n_no_n);
    std::fprintf(stderr, "[oracle] normal: exact=%lu no=%lu\n", (unsigned long)n_exact_n, (unsigned long)n_no_n);
    pair_sample(tfa, level, flank, all_pairs, n_exact_t, n_no_t);
    std::fprintf(stderr, "[oracle] tumor:  exact=%lu no=%lu\n", (unsigned long)n_exact_t, (unsigned long)n_no_t);
    std::fprintf(stderr, "[oracle] total pairs collected: %lu\n", (unsigned long)all_pairs.size());

    // Build partner lookup: fwd_sig → rc_sig, rc_sig → fwd_sig.
    std::unordered_map<uint64_t, uint64_t> fwd2rc, rc2fwd;
    fwd2rc.reserve(all_pairs.size() * 2);
    rc2fwd.reserve(all_pairs.size() * 2);
    for (const auto& p : all_pairs) {
        fwd2rc[p.first]  = p.second;
        rc2fwd[p.second] = p.first;
    }
    std::fprintf(stderr, "[oracle] partner table sizes: fwd2rc=%zu rc2fwd=%zu\n", fwd2rc.size(), rc2fwd.size());

    // For each sig in counts.bin, look up its RC partner.
    auto cls = [&](uint64_t i) -> int {
        const Prov& p = prov[i];
        const bool n = (p.n_fwd + p.n_rc) > 0;
        const bool t = (p.t_fwd + p.t_rc) > 0;
        if (n && t) return 0;
        if (n)      return 1;
        if (t)      return 2;
        return 3;
    };
    const char* CN[] = {"BOTH", "N_only", "T_only", "TRUE_PHANTOM"};

    uint64_t tot[4] = {0,0,0,0};
    uint64_t partner_in_table[4] = {0,0,0,0};
    uint64_t partner_in_counts[4] = {0,0,0,0};
    for (uint64_t i = 0; i < K; i++) {
        const int c = cls(i);
        tot[c]++;
        const uint64_t sig = sigs[i].sig;
        uint64_t partner;
        auto itf = fwd2rc.find(sig);
        if (itf != fwd2rc.end()) { partner = itf->second; partner_in_table[c]++; }
        else {
            auto itr = rc2fwd.find(sig);
            if (itr != rc2fwd.end()) { partner = itr->second; partner_in_table[c]++; }
            else continue;
        }
        if (sigset.count(partner)) partner_in_counts[c]++;
    }
    std::printf("\n== Positional RC pairing oracle ==\n");
    std::printf("  pairs in oracle table: %lu  (from %lu fwd cores normal+tumor)\n",
                (unsigned long)all_pairs.size(), (unsigned long)(n_exact_n + n_exact_t + n_no_n + n_no_t));
    std::printf("\n  class           n_sigs       has_partner_oracle    partner_in_counts.bin (of those with oracle partner)\n");
    for (int c = 0; c < 4; c++) {
        if (tot[c] == 0) continue;
        const double pct_tbl = 100.0 * partner_in_table[c] / tot[c];
        const double pct_cts = partner_in_table[c] ? 100.0 * partner_in_counts[c] / partner_in_table[c] : 0.0;
        std::printf("  %-14s  %-10lu   %lu (%.2f%%)        %lu (%.2f%%)\n",
                    CN[c], (unsigned long)tot[c],
                    (unsigned long)partner_in_table[c], pct_tbl,
                    (unsigned long)partner_in_counts[c], pct_cts);
    }
    return 0;
}
