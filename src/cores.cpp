#include "svp.h"
#include "MurmurHash3.h"

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>

extern "C" {
#include "lps.h"
}

static inline uint64_t fmix64(uint64_t k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

void extract_cores(const std::string& seq, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, int level, std::vector<Core>& out) {
    out.clear();
    for (size_t r = 0; r < rec_start.size(); r++) {
        const uint64_t s = rec_start[r], e = rec_end[r];
        if (e - s > (uint64_t)INT_MAX) { std::fprintf(stderr, "[FATAL] record %zu length %lu exceeds INT_MAX\n", r, (unsigned long)(e - s)); std::abort(); }
        int n; const uint64_t *st, *en; const uint32_t *lb;
        extract_fast(seq.data() + s, (int)(e - s), s, level, &n, &st, &en, &lb);
        const size_t base = out.size();
        out.resize(base + (size_t)n);
        for (int i = 0; i < n; i++) out[base + i] = {st[i], en[i], lb[i]};
    }
}

static void flank_sigs(const std::string& seq, const std::vector<Core>& cs, int m, std::vector<uint64_t>& out) {
    out.resize(cs.size());
    const uint64_t n = seq.size(), mu = (uint64_t)m;
    for (size_t i = 0; i < cs.size(); i++) {
        const uint64_t L = cs[i].start > mu ? cs[i].start - mu : 0;
        const uint64_t R = cs[i].end + mu < n ? cs[i].end + mu : n;
        out[i] = murmur64(seq.data() + L, (int)(R - L), cs[i].label);
    }
}

static void mark_collisions(const std::vector<uint64_t>& sig, std::vector<char>& collide) {
    collide.assign(sig.size(), 0);
    std::unordered_map<uint64_t, int32_t> seen;
    seen.reserve(sig.size() * 2);
    for (size_t i = 0; i < sig.size(); i++) {
        int32_t& v = seen[sig[i]];
        if (v == 0) { v = (int32_t)(i + 1); continue; }
        collide[i] = 1;
        if (v > 0) { collide[v - 1] = 1; v = -1; }
    }
}

static void rescue_colliders(const std::vector<Core>& cs, const std::vector<uint64_t>& s, const std::vector<char>& collide, int R, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, std::vector<uint64_t>& out) {
    const uint64_t GOLD = 0x9E3779B97F4A7C15ULL;
    size_t ci = 0;
    for (size_t k = 0; k < rec_start.size(); k++) {
        while (ci < cs.size() && cs[ci].start < rec_start[k]) ci++;
        const size_t lo = ci;
        while (ci < cs.size() && cs[ci].end <= rec_end[k]) ci++;
        const size_t hi = ci;
        for (size_t i = lo; i < hi; i++) {
            if (!collide[i]) continue;
            uint64_t h = (uint64_t)R * GOLD;
            for (int j = (int)i - R; j <= (int)i + R; j++) {
                const uint64_t v = (j < (int)lo || j >= (int)hi) ? 0ULL : s[j];
                h = fmix64(h ^ v);
            }
            out[i] = h;
        }
    }
}

void compute_sigs(const std::string& seq, const std::vector<Core>& cs, int m, int context_radius, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, std::vector<uint64_t>& out) {
    flank_sigs(seq, cs, m, out);
    if (context_radius == 0 || cs.empty()) return;
    std::vector<char> collide;
    mark_collisions(out, collide);
    const std::vector<uint64_t> s = out;
    rescue_colliders(cs, s, collide, context_radius, rec_start, rec_end, out);
}
