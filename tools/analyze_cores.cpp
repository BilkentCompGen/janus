#include "../src/svp.h"
#include "../src/MurmurHash3.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <string>

extern "C" {
#include "encoding.h"
}

static void compute_three_sigs(const std::string& seq, const std::vector<Core>& cs, int m, std::vector<uint64_t>& kL, std::vector<uint64_t>& kF, std::vector<uint64_t>& kC, std::vector<uint64_t>& kD) {
    const uint64_t N = cs.size();
    kL.resize(N); kF.resize(N); kC.resize(N); kD.resize(N);
    const uint64_t mu = (uint64_t)m, B = seq.size();
    for (uint64_t i = 0; i < N; i++) {
        const uint64_t L = cs[i].start > mu ? cs[i].start - mu : 0;
        const uint64_t R = cs[i].end + mu < B ? cs[i].end + mu : B;
        kL[i] = (uint64_t)cs[i].label;
        kF[i] = murmur64(seq.data() + L, (int)(R - L), 0);
        kC[i] = murmur64(seq.data() + cs[i].start, (int)(cs[i].end - cs[i].start), 0);
        kD[i] = murmur64(seq.data() + L, (int)(R - L), cs[i].label);
    }
}

struct Counts { uint64_t total, distinct, singletons, dup_keys, dup_cores; uint32_t max_rep; };

static Counts build_counts(const std::vector<uint64_t>& keys, std::unordered_map<uint64_t, uint32_t>& ct) {
    ct.clear();
    ct.reserve(keys.size() * 2);
    for (uint64_t k : keys) ct[k]++;
    Counts c{keys.size(), ct.size(), 0, 0, 0, 0};
    for (const auto& kv : ct) {
        if (kv.second == 1) c.singletons++;
        else { c.dup_keys++; c.dup_cores += kv.second; }
        if (kv.second > c.max_rep) c.max_rep = kv.second;
    }
    return c;
}

static void print_summary(const char* name, const Counts& c) {
    const double t = (double)c.total;
    std::fprintf(stdout, "  %-12s  cores=%-12lu  distinct=%-12lu (%.3f%%)  singletons=%-12lu (%.3f%%)  dup-keys=%-10lu  cores-in-dups=%-12lu  max-rep=%u\n", name, (unsigned long)c.total, (unsigned long)c.distinct, 100.0*(double)c.distinct/t, (unsigned long)c.singletons, 100.0*(double)c.singletons/t, (unsigned long)c.dup_keys, (unsigned long)c.dup_cores, c.max_rep);
}

static void print_hist(const char* name, const std::unordered_map<uint64_t, uint32_t>& ct) {
    uint64_t b[8] = {0};
    for (const auto& kv : ct) {
        const uint32_t c = kv.second;
        const int idx = (c == 1) ? 0 : (c == 2) ? 1 : (c <= 4) ? 2 : (c <= 8) ? 3 : (c <= 16) ? 4 : (c <= 64) ? 5 : (c <= 256) ? 6 : 7;
        b[idx]++;
    }
    std::fprintf(stdout, "  hist(%-12s)  rep=1: %-12lu  =2: %-10lu  3-4: %-10lu  5-8: %-9lu  9-16: %-9lu  17-64: %-9lu  65-256: %-9lu  257+: %lu\n", name, b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
}

static uint64_t cross_singletons(const std::unordered_map<uint64_t, uint32_t>& a, const std::unordered_map<uint64_t, uint32_t>& b) {
    const auto& small = a.size() < b.size() ? a : b;
    const auto& big   = a.size() < b.size() ? b : a;
    uint64_t n = 0;
    for (const auto& kv : small) {
        if (kv.second != 1) continue;
        const auto it = big.find(kv.first);
        if (it != big.end() && it->second == 1) n++;
    }
    return n;
}

static void parse_args(int argc, char** argv, int& level, int& flank) {
    level = 5; flank = 10;
    for (int i = 3; i < argc; i++) {
        const bool is_l = std::strcmp(argv[i], "--level") == 0;
        const bool is_f = std::strcmp(argv[i], "--flank") == 0;
        if (!(is_l || is_f) || i + 1 >= argc) { std::fprintf(stderr, "[FATAL] arg: %s\n", argv[i]); std::abort(); }
        const int v = std::atoi(argv[i + 1]);
        if (is_l) level = v; else flank = v;
        i++;
    }
}

static void load_and_extract(const char* path, int flank, int level, std::string& seq, std::vector<Core>& cs, const char* tag) {
    std::vector<uint64_t> rs, re;
    read_fasta_flat(path, flank, seq, rs, re);
    std::fprintf(stdout, "[%s] %s  buf=%lu bp  records=%lu\n", tag, path, seq.size(), rs.size());
    std::fflush(stdout);
    extract_cores(seq, rs, re, level, cs);
    std::fprintf(stdout, "[%s] cores=%lu  (level=%d, flank=%d)\n", tag, cs.size(), level, flank);
    std::fflush(stdout);
}

static void report_one_side(const char* tag, const char* path, const std::vector<uint64_t>& kL, const std::vector<uint64_t>& kF, const std::vector<uint64_t>& kC, const std::vector<uint64_t>& kD, std::unordered_map<uint64_t, uint32_t>& cL, std::unordered_map<uint64_t, uint32_t>& cF, std::unordered_map<uint64_t, uint32_t>& cC, std::unordered_map<uint64_t, uint32_t>& cD) {
    std::fprintf(stdout, "\n=== %s: %s ===\n", tag, path);
    print_summary("label-only",  build_counts(kL, cL));
    print_summary("flank-only",  build_counts(kF, cF));
    print_summary("core-only",   build_counts(kC, cC));
    print_summary("label+flank", build_counts(kD, cD));
    print_hist("label-only",  cL);
    print_hist("flank-only",  cF);
    print_hist("core-only",   cC);
    print_hist("label+flank", cD);
    std::fflush(stdout);
}

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: analyze_cores <normal.fa[.gz]> <tumor.fa[.gz]> [--level N=5] [--flank M=10]\n"); return 2; }
    int level, flank;
    parse_args(argc, argv, level, flank);
    LCP_INIT();

    std::string N, T;
    std::vector<Core> ncs, tcs;
    load_and_extract(argv[1], flank, level, N, ncs, "ref");
    load_and_extract(argv[2], flank, level, T, tcs, "tum");

    std::vector<uint64_t> nL, nF, nC, nD, tL, tF, tC, tD;
    compute_three_sigs(N, ncs, flank, nL, nF, nC, nD);
    compute_three_sigs(T, tcs, flank, tL, tF, tC, tD);

    std::unordered_map<uint64_t, uint32_t> cnL, cnF, cnC, cnD, ctL, ctF, ctC, ctD;
    report_one_side("ref", argv[1], nL, nF, nC, nD, cnL, cnF, cnC, cnD);
    report_one_side("tum", argv[2], tL, tF, tC, tD, ctL, ctF, ctC, ctD);

    std::fprintf(stdout, "\n=== anchor candidates (singletons in BOTH ref and tum) ===\n");
    std::fprintf(stdout, "  label-only:   %lu\n", (unsigned long)cross_singletons(cnL, ctL));
    std::fprintf(stdout, "  flank-only:   %lu\n", (unsigned long)cross_singletons(cnF, ctF));
    std::fprintf(stdout, "  core-only:    %lu\n", (unsigned long)cross_singletons(cnC, ctC));
    std::fprintf(stdout, "  label+flank:  %lu  <-- current svp\n", (unsigned long)cross_singletons(cnD, ctD));
    return 0;
}
