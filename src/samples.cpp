#include "samples.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

namespace {

void die(const std::string& m) { std::fprintf(stderr, "[FATAL] samples: %s\n", m.c_str()); std::abort(); }

size_t fields(char* line, char* out[3]) {
    size_t n = 0;
    for (char* p = line; *p && n < 3;) {
        while (*p == ' ' || *p == '\t') p++;
        if (!*p || *p == '\n' || *p == '\r') break;
        out[n++] = p;
        while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
        if (*p) *p++ = 0;
    }
    return n;
}

void check(std::vector<Sample>& s) {
    if (s.size() < 2) die("need one normal and at least one tumor row");
    if (s.size() > 64) die("at most 64 samples: the per-chunk carriage mask is 64-bit");
    size_t nn = 0;
    for (const Sample& x : s) nn += x.normal;
    if (nn != 1) die("need exactly one normal row, found " + std::to_string(nn));
    for (size_t i = 0; i < s.size(); i++) {
        if (s[i].name.empty()) die("empty sample name");

        if (s[i].name.find_first_of(";|/ \t") != std::string::npos)
            die("sample name may not contain ; | / space or tab: " + s[i].name);
        for (size_t j = 0; j < i; j++)
            if (s[i].name == s[j].name) die("duplicate sample name: " + s[i].name);
        if (access(s[i].path.c_str(), R_OK)) die("cannot read " + s[i].path);
    }
}

}

std::vector<Sample> samples_read(const std::string& tsv) {
    FILE* f = std::fopen(tsv.c_str(), "r");
    if (!f) die("open " + tsv);
    std::vector<Sample> raw;
    char* buf = nullptr;
    size_t cap = 0;
    bool hdr = false;
    while (getline(&buf, &cap, f) > 0) {
        char* c[3];
        if (buf[0] == '#') continue;
        const size_t nf = fields(buf, c);
        if (nf == 0) continue;
        if (nf != 3) die("three columns required (name role path), got " + std::to_string(nf));
        if (!hdr) {
            hdr = true;
            if (std::strcmp(c[0], "name") || std::strcmp(c[1], "role") || std::strcmp(c[2], "path"))
                die("first row must be the header: name role path");
            continue;
        }
        const bool nrm = std::strcmp(c[1], "normal") == 0;
        if (!nrm && std::strcmp(c[1], "tumor")) die("role must be normal or tumor: " + std::string(c[1]));
        raw.push_back({c[0], c[2], (uint8_t)nrm});
    }
    std::free(buf);
    std::fclose(f);
    check(raw);

    std::vector<Sample> out;
    for (const Sample& x : raw) if (x.normal) out.push_back(x);
    for (const Sample& x : raw) if (!x.normal) out.push_back(x);
    std::fprintf(stderr, "[samples] %s: 1 normal + %zu tumor\n", tsv.c_str(), out.size() - 1);
    return out;
}

void samples_lock(const std::vector<Sample>& s, const std::string& tsv) {
    if (!access(tsv.c_str(), R_OK)) {
        const std::vector<Sample> old = samples_read(tsv);
        if (old.size() != s.size()) die("sheet has " + std::to_string(s.size()) + " rows, " + tsv +
                                        " was built from " + std::to_string(old.size()));
        for (size_t i = 0; i < s.size(); i++)
            if (old[i].name != s[i].name || old[i].path != s[i].path || old[i].normal != s[i].normal)
                die("sheet row " + std::to_string(i) + " differs from " + tsv +
                    ": the artifacts in this dir were built from a different sample order");
        return;
    }
    FILE* f = std::fopen((tsv + ".tmp").c_str(), "w");
    if (!f) die("open " + tsv + ".tmp");
    std::fprintf(f, "name\trole\tpath\n");
    for (const Sample& x : s) std::fprintf(f, "%s\t%s\t%s\n", x.name.c_str(), x.normal ? "normal" : "tumor", x.path.c_str());
    std::fclose(f);
    if (std::rename((tsv + ".tmp").c_str(), tsv.c_str())) die("rename " + tsv);
}
