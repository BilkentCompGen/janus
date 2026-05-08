#include "svp.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

extern "C" {
#include "encoding.h"
}

static int parse_int(const char* s, const char* name) {
    char* end;
    const long v = std::strtol(s, &end, 10);
    if (*end || v <= 0 || v > 1000000) { std::fprintf(stderr, "[FATAL] bad %s: %s\n", name, s); std::abort(); }
    return (int)v;
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: svp <normal.fa[.gz]> <tumor.fa[.gz]> <out.tsv> [--level N=5] [--flank M=10]\n"); return 2; }
    int level = 5, flank = 10;
    for (int i = 4; i < argc; i++) {
        const bool is_level = std::strcmp(argv[i], "--level") == 0;
        const bool is_flank = std::strcmp(argv[i], "--flank") == 0;
        if (!(is_level || is_flank) || i + 1 >= argc) { std::fprintf(stderr, "[FATAL] unknown arg: %s\n", argv[i]); std::abort(); }
        const int v = parse_int(argv[i + 1], is_level ? "level" : "flank");
        if (is_level) level = v; else flank = v;
        i++;
    }

    LCP_INIT();
    std::string N, T;
    std::vector<uint64_t> n_rs, n_re, t_rs, t_re;
    read_fasta_flat(argv[1], flank, N, n_rs, n_re);
    read_fasta_flat(argv[2], flank, T, t_rs, t_re);
    std::fprintf(stderr, "[INFO] N=%lu (%lu rec)  T=%lu (%lu rec)\n", N.size(), n_rs.size(), T.size(), t_rs.size());

    std::vector<Core> ncs, tcs;
    extract_cores(N, n_rs, n_re, level, ncs);
    extract_cores(T, t_rs, t_re, level, tcs);
    std::fprintf(stderr, "[INFO] cores N=%lu T=%lu\n", ncs.size(), tcs.size());

    std::vector<uint64_t> nsig, tsig;
    compute_sigs(N, ncs, flank, nsig);
    compute_sigs(T, tcs, flank, tsig);

    std::vector<Anchor> anchors;
    unique_anchors(nsig, tsig, anchors);
    std::fprintf(stderr, "[INFO] anchors=%lu\n", anchors.size());

    std::vector<int> kept;
    chain(anchors, kept);
    std::fprintf(stderr, "[INFO] chain=%lu\n", kept.size());

    std::vector<Bubble> bubs;
    emit_bubbles(N, ncs, T, tcs, anchors, kept, bubs);
    std::fprintf(stderr, "[INFO] bubbles=%lu\n", bubs.size());

    write_bubbles(argv[3], bubs);
    return 0;
}
