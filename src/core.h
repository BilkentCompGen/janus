#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "lps.h"

class Bases {
public:
    Bases();
    bool plain(const uint8_t* s, uint32_t n) const;
    void pack(const uint8_t* s, uint32_t n, uint8_t* d) const;

private:
    uint8_t enc[256], bad[256];
};

extern const Bases BASES;

struct Block;

uint64_t core_hash(const uint64_t* w, uint32_t k);

class Read {
public:
    uint32_t len, k, n, nf;
    std::vector<uint32_t> a, b;
    void squeeze(const uint8_t* rec);
    void load(const uint8_t* exc, int w);
    void store(std::vector<char>& out, int w) const;
    const uint8_t* next(const Block& bl, size_t* o, size_t* e, int w);
    uint64_t hash(uint32_t i);
    uint64_t head() const { return hs[0]; }
    uint32_t forward(uint32_t i) const { return i < nf ? a[i] : len - b[i]; }
    std::string bases(const uint8_t* rec, uint32_t lo, uint32_t hi, uint32_t rc) const;
    void letters(bool reverse, char* out) const;
    void starts(uint32_t* map) const;

private:
    std::vector<uint64_t> hs, rs, mk, word;
    std::vector<uint32_t> cum;
    std::vector<uint8_t> raw;
    uint32_t pad;
    uint32_t rank(uint32_t x) const;
    uint32_t cut(uint32_t i, uint64_t* out) const;
};

static const int LEV = 5;

class Parser {
public:
    Parser() : mem(nullptr) {}
    ~Parser();
    void cores(Read& r, const uint8_t* rec);

private:
    lcp_soa o[2];
    uint64_t* mem;
    uint32_t *lab, *lsz;
    std::vector<uint32_t> map;
    std::vector<char> hs;
    void grow(Read& r);
    uint32_t parse(Read& r, uint32_t at);
};
