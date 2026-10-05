#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <map>
#include "cmd.h"
#include "files.h"
#include "io.h"

static const uint32_t MING = 50;
static const uint32_t MAXG = 20000;

int cmd_lift(const Args& a) {
    Fai fai(a.data);
    Lines in(a.in);
    Paf p;
    FILE* fo = std::fopen(a.out, "w");
    while (in.next()) {
        p.parse(in.s);
        char* c = std::strrchr(p.t, ':');
        *c = 0;
        uint32_t lo = std::strtoul(c + 1, nullptr, 10);
        std::unordered_map<std::string, FaiEntry>::const_iterator it = fai.m.find(p.t);
        unsigned long len = it == fai.m.end() ? 0 : it->second.len;
        std::fprintf(fo, "%s\t%u\t%u\t%u\t%s\t%s\t%lu\t%lu\t%lu%s", p.q, p.ql, p.qs, p.qe, p.st, p.t, len, (unsigned long)(p.ts + lo), (unsigned long)(p.te + lo), in.s + p.k9);
    }
    std::fclose(fo);
    return 0;
}

int cmd_clip(const Args& a) {
    Fasta f((std::string(a.in) + ".fa").c_str());
    std::map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> cov;
    {
        Lines in(a.reads);
        Paf p;
        while (in.next())
            if (p.parse(in.s) >= 4) cov[p.q].push_back({p.qs, p.qe});
    }
    FILE* fo = std::fopen(a.out, "w");
    for (size_t i = 0; i < f.name.size(); i++) {
        std::map<std::string, std::vector<std::pair<uint32_t, uint32_t>>>::iterator it = cov.find(f.name[i]);
        if (it == cov.end()) continue;
        std::vector<std::pair<uint32_t, uint32_t>>& c = it->second;
        std::sort(c.begin(), c.end());
        std::vector<std::pair<uint32_t, uint32_t>> g;
        uint32_t at = 0, L = f.seq[i].size();
        for (std::pair<uint32_t, uint32_t>& p : c) {
            if (p.first >= at + MING) g.push_back({at, p.first});
            at = std::max(at, p.second);
        }
        if (L >= at + MING) g.push_back({at, L});
        for (std::pair<uint32_t, uint32_t>& p : g)
            if (p.second - p.first <= MAXG) std::fprintf(fo, ">%s|%u|%u\n%s\n", f.name[i].c_str(), p.first, p.second, f.seq[i].substr(p.first, p.second - p.first).c_str());
    }
    std::fclose(fo);
    return 0;
}

int cmd_patch(const Args& a) {
    Fasta f((std::string(a.in) + ".fa").c_str());
    std::map<std::string, uint32_t> len;
    for (size_t i = 0; i < f.name.size(); i++) len[f.name[i]] = f.seq[i].size();
    FILE* fo = std::fopen(a.out, "w");
    {
        Lines own(a.data);
        while (own.next()) std::fputs(own.s, fo);
    }
    Lines in(a.reads);
    Paf p;
    while (in.next()) {
        p.parse(in.s);
        char* b = std::strrchr(p.q, '|');
        *b = 0;
        char* c = std::strrchr(p.q, '|');
        *c = 0;
        uint32_t off = std::strtoul(c + 1, nullptr, 10);
        std::fprintf(fo, "%s\t%u\t%u\t%u%s", p.q, len[p.q], off + p.qs, off + p.qe, in.s + p.k4);
    }
    std::fclose(fo);
    return 0;
}
