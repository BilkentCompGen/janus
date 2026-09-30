#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include "cmd.h"

// lift: the pieces' local alignments, made against reference windows named <chromosome>:<start>, back in chromosome
// coordinates; a piece's unaligned sides are the clip stage's job
int cmd_lift(const char* paf, const char* ref, const char* out) {
    std::map<std::string, uint64_t> len;
    char nm[256];
    unsigned long l, z;
    FILE* fi = std::fopen((std::string(ref) + ".fai").c_str(), "r");
    while (std::fscanf(fi, "%255s %lu %lu %lu %lu", nm, &l, &z, &z, &z) == 5) len[nm] = l;
    std::fclose(fi);
    fi = std::fopen(paf, "r");
    FILE* fo = std::fopen(out, "w");
    char* ln = nullptr;
    size_t cap = 0;
    for (ssize_t n; (n = getline(&ln, &cap, fi)) > 0;) {
        char q[256], s[4], tn[512];
        unsigned qn, qs, qe, tl, ts, te;
        int k;
        std::sscanf(ln, "%255s %u %u %u %3s %511s %u %u %u%n", q, &qn, &qs, &qe, s, tn, &tl, &ts, &te, &k);
        char* c = std::strrchr(tn, ':');
        *c = 0;
        uint32_t lo = std::strtoul(c + 1, nullptr, 10);
        std::fprintf(fo, "%s\t%u\t%u\t%u\t%s\t%s\t%lu\t%lu\t%lu%s", q, qn, qs, qe, s, tn, (unsigned long)len[tn], (unsigned long)(ts + lo), (unsigned long)(te + lo), ln + k);
    }
    std::free(ln);
    std::fclose(fi);
    std::fclose(fo);
    return 0;
}
