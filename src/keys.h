#pragma once
#include <algorithm>
#include <cstdint>
#include <vector>
#include "files.h"
#include "reads.h"

class Combos {
public:
    uint64_t S, Z, ZT;
    uint8_t full;
    std::vector<uint64_t> seen, key, best;
    std::vector<uint32_t> n;
    std::vector<uint8_t> branch, keeps, on;
    std::vector<uint32_t> seg, fork;
    Combos(const Graph& g, const Adjacency& a);
    uint64_t keep(const std::vector<uint8_t>& sd);
    void table(uint64_t sites);
    void grow();
    bool kept(uint64_t k) const;
    void find(Chain& c, uint32_t m) const;
    void tally(Chain& c, uint32_t m, uint64_t s, uint64_t rid);
    static uint32_t in(uint64_t k) { return (k >> 32) - 1; }
    static uint32_t out(uint64_t k) { return k & 0xffffffff; }
    static uint64_t centrality(const Chain& c, const Hit& h, uint64_t rid) {
        uint32_t st = c.start(h.a), en = c.end(h.b);
        return (uint64_t)std::min(st, c.r.len - en) << 32 | (0xffffffffu - (uint32_t)rid);
    }
    static uint64_t central_read(uint64_t s, uint64_t score) { return s << 40 | (0xffffffffu - (uint32_t)score); }

private:
    const Graph& g;
    const Adjacency& a;
    uint32_t edge(uint32_t u, uint32_t v) const;
    static uint64_t home(uint64_t k, uint64_t z) { return (k * 0x9e3779b97f4a7c15ull) >> 20 & (z - 1); }
    void remember(uint64_t k);
    bool known(uint64_t k) const;
    void add(uint64_t k, uint64_t s, uint64_t sc);
};
