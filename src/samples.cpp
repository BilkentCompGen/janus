#include "samples.h"
#include "log.h"

#include <stdlib.h>
#include <cstring>
#include <unistd.h>

static size_t fields(char* line, char* out[3]) {
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

static void check(std::vector<Sample>& s) {
    if (s.size() < 2) die("samples", "need one normal and at least one tumor row");
    if (s.size() > 64) die("samples", "at most 64 samples: the per-chunk carriage mask is 64-bit");
    size_t nn = 0;
    for (const Sample& x : s) nn += x.normal;
    if (nn != 1) die("samples", "need exactly one normal row, found %zu", nn);
    for (size_t i = 0; i < s.size(); i++) {
        if (s[i].name.empty()) die("samples", "empty sample name");
        if (s[i].name.find_first_of(";|/ \t") != std::string::npos)
            die("samples", "sample name may not contain ; | / space or tab: %s", s[i].name.c_str());
        for (size_t j = 0; j < i; j++)
            if (s[i].name == s[j].name) die("samples", "duplicate sample name: %s", s[i].name.c_str());
        if (access(s[i].path.c_str(), R_OK)) die("samples", "cannot read %s", s[i].path.c_str());
    }
}

std::vector<Sample> samples_read(const std::string& tsv) {
    FILE* f = xopen(tsv.c_str(), "r");
    std::vector<Sample> raw;
    char* buf = nullptr;
    size_t cap = 0;
    bool hdr = false;
    while (getline(&buf, &cap, f) > 0) {
        char* c[3];
        if (buf[0] == '#') continue;
        const size_t nf = fields(buf, c);
        if (nf == 0) continue;
        if (nf != 3) die("samples", "three columns required (name role path), got %zu", nf);
        if (!hdr) {
            hdr = true;
            if (std::strcmp(c[0], "name") || std::strcmp(c[1], "role") || std::strcmp(c[2], "path"))
                die("samples", "first row must be the header: name role path");
            continue;
        }
        const bool nrm = std::strcmp(c[1], "normal") == 0;
        if (!nrm && std::strcmp(c[1], "tumor")) die("samples", "role must be normal or tumor: %s", c[1]);
        raw.push_back({c[0], c[2], (uint8_t)nrm});
    }
    std::free(buf);
    std::fclose(f);
    check(raw);

    std::vector<Sample> out;
    for (const Sample& x : raw) if (x.normal) out.push_back(x);
    for (const Sample& x : raw) if (!x.normal) out.push_back(x);
    say("samples", "%s: 1 normal + %zu tumor", tsv.c_str(), out.size() - 1);
    return out;
}

void samples_lock(const std::vector<Sample>& s, const std::string& tsv) {
    if (!access(tsv.c_str(), R_OK)) {
        const std::vector<Sample> old = samples_read(tsv);
        if (old.size() != s.size()) die("samples", "sheet has %zu rows, %s was built from %zu",
                                        s.size(), tsv.c_str(), old.size());
        for (size_t i = 0; i < s.size(); i++)
            if (old[i].name != s[i].name || old[i].path != s[i].path || old[i].normal != s[i].normal)
                die("samples", "sheet row %zu differs from %s: the artifacts in this dir were built from a different sample order",
                    i, tsv.c_str());
        return;
    }
    FILE* f = xopen((tsv + ".tmp").c_str(), "w");
    std::fprintf(f, "name\trole\tpath\n");
    for (const Sample& x : s)
        std::fprintf(f, "%s\t%s\t%s\n", x.name.c_str(), x.normal ? "normal" : "tumor", x.path.c_str());
    atomic_close(f, (tsv + ".tmp").c_str(), tsv.c_str());
}
