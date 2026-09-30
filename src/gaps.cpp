#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include "cmd.h"
#include "files.h"

static const uint32_t MING = 50;    // shortest unaligned stretch of a piece worth aligning on its own (an SV-sized insert)
static const uint32_t MAXG = 20000; // longest one; beyond that the stretch is not a junction insert

// clip: the unaligned stretches of the pieces (between two alignments, or a clipped end) of MING..MAXG bp, as a FASTA named
// <piece>|<start>|<end>. Such a stretch is too short to survive as a supplementary alignment of a 15 kb piece, but on its own
// (short-read mode) it places the partner breakpoint of a complex junction.
int cmd_clip(const char* in, const char* paf, const char* out) {
    seqs f;
    read_fasta((std::string(in) + ".fa").c_str(), &f);
    std::map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> cov;
    char q[256];
    unsigned qs, qe;
    FILE* fi = std::fopen(paf, "r");
    char* ln = nullptr;
    size_t cap = 0;
    for (ssize_t n; (n = getline(&ln, &cap, fi)) > 0;) if (std::sscanf(ln, "%255s %*u %u %u", q, &qs, &qe) == 3) cov[q].push_back({qs, qe});
    std::free(ln);
    std::fclose(fi);
    FILE* fo = std::fopen(out, "w");
    for (size_t i = 0; i < f.name.size(); i++) {
        auto it = cov.find(f.name[i]);
        if (it == cov.end()) continue;
        std::vector<std::pair<uint32_t, uint32_t>>& c = it->second;
        std::sort(c.begin(), c.end());
        std::vector<std::pair<uint32_t, uint32_t>> g;
        uint32_t at = 0, L = f.seq[i].size();
        for (auto& p : c) { if (p.first >= at + MING) g.push_back({at, p.first}); at = std::max(at, p.second); }
        if (L >= at + MING) g.push_back({at, L});
        for (auto& p : g) if (p.second - p.first <= MAXG) std::fprintf(fo, ">%s|%u|%u\n%s\n", f.name[i].c_str(), p.first, p.second, f.seq[i].substr(p.first, p.second - p.first).c_str());
    }
    std::fclose(fo);
    return 0;
}

// patch: the pieces' own alignments (own), then the stretches' alignments back in piece coordinates (query name, length,
// start, end)
int cmd_patch(const char* in, const char* paf, const char* own, const char* out) {
    seqs f;
    read_fasta((std::string(in) + ".fa").c_str(), &f);
    std::map<std::string, uint32_t> len;
    for (size_t i = 0; i < f.name.size(); i++) len[f.name[i]] = f.seq[i].size();
    FILE* fo = std::fopen(out, "w");
    for (const char* src : {own, paf}) {
        FILE* fi = std::fopen(src, "r");
        char* ln = nullptr;
        size_t cap = 0;
        for (ssize_t n; (n = getline(&ln, &cap, fi)) > 0;) {
            if (src == own) { std::fputs(ln, fo); continue; }
            char q[300];
            unsigned qs, qe;
            int k;
            std::sscanf(ln, "%299s %*u %u %u%n", q, &qs, &qe, &k);
            char* b = std::strrchr(q, '|');
            *b = 0;
            char* a = std::strrchr(q, '|');
            *a = 0;
            uint32_t off = std::strtoul(a + 1, nullptr, 10);
            std::fprintf(fo, "%s\t%u\t%u\t%u%s", q, len[q], off + qs, off + qe, ln + k);
        }
        std::free(ln);
        std::fclose(fi);
    }
    std::fclose(fo);
    return 0;
}
