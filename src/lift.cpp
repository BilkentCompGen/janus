#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include "cmd.h"
#include "files.h"

static const uint32_t EXPL = 100; // piece bases the local alignments may leave uncovered before the piece is realigned genome-wide

int cmd_lift(const char* paf, const char* on, const char* ref, const char* out) {
    std::map<std::string, uint64_t> len;
    std::map<std::string, std::vector<std::string>> rows;
    std::map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> cov;
    char nm[256];
    unsigned long l, z;
    FILE* fi = std::fopen((std::string(ref) + ".fai").c_str(), "r");
    while (std::fscanf(fi, "%255s %lu %lu %lu %lu", nm, &l, &z, &z, &z) == 5) len[nm] = l;
    std::fclose(fi);
    fi = std::fopen(paf, "r");
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
        rows[q].push_back(std::string(q) + "\t" + std::to_string(qn) + "\t" + std::to_string(qs) + "\t" + std::to_string(qe) + "\t" + s + "\t" + tn + "\t" + std::to_string(len[tn]) + "\t" + std::to_string(ts + lo) + "\t" + std::to_string(te + lo) + (ln + k));
        cov[q].push_back({qs, qe});
    }
    std::free(ln);
    std::fclose(fi);
    seqs f;
    read_fasta(on, &f);
    std::string o(out);
    FILE *fp = std::fopen((o + ".paf").c_str(), "w"), *fm = std::fopen((o + ".more.fa").c_str(), "w");
    for (size_t j = 0; j < f.name.size(); j++) {
        std::vector<std::pair<uint32_t, uint32_t>>& v = cov[f.name[j]];
        std::sort(v.begin(), v.end());
        uint32_t got = 0, end = 0;
        for (auto& x : v) { if (x.second > end) got += x.second - std::max(x.first, end); end = std::max(end, x.second); }
        if (f.seq[j].size() - got >= EXPL) std::fprintf(fm, ">%s\n%s\n", f.name[j].c_str(), f.seq[j].c_str());
        else for (std::string& r : rows[f.name[j]]) std::fputs(r.c_str(), fp);
    }
    std::fclose(fp);
    std::fclose(fm);
    return 0;
}
