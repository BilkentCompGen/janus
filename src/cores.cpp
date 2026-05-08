#include "svp.h"
#include "MurmurHash3.h"

#include <climits>
#include <cstdio>
#include <cstdlib>

extern "C" {
#include "lps.h"
}

void extract_cores(const std::string& seq, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, int level, std::vector<Core>& out) {
    out.clear();
    for (size_t r = 0; r < rec_start.size(); r++) {
        const uint64_t s = rec_start[r], e = rec_end[r];
        if (e - s > (uint64_t)INT_MAX) { std::fprintf(stderr, "[FATAL] record %zu length %lu exceeds INT_MAX\n", r, (unsigned long)(e - s)); std::abort(); }
        struct lps p;
        init_lps_offset(&p, seq.data() + s, (int)(e - s), s);
        lps_deepen(&p, level);
        const size_t base = out.size();
        out.resize(base + p.size);
        for (int i = 0; i < p.size; i++) out[base + i] = {p.cores[i].start, p.cores[i].end, p.cores[i].label};
        free_lps(&p);
    }
}

void compute_sigs(const std::string& seq, const std::vector<Core>& cs, int m, std::vector<uint64_t>& out) {
    out.resize(cs.size());
    const uint64_t n = seq.size();
    const uint64_t mu = (uint64_t)m;
    for (size_t i = 0; i < cs.size(); i++) {
        const uint64_t L = cs[i].start > mu ? cs[i].start - mu : 0;
        const uint64_t R = cs[i].end + mu < n ? cs[i].end + mu : n;
        out[i] = murmur64(seq.data() + L, (int)(R - L), cs[i].label);
    }
}
