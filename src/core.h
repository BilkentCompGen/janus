#pragma once
#include <cstddef>
#include <cstdint>

struct rdcore {
    uint64_t *hs, *rs, *mk; uint8_t *raw, *u; uint32_t *cum, *a, *b;
    uint32_t len, k, pad, n, nf, cap, ccap;
};

size_t exc_size(const char* p, size_t n);
void collapse(rdcore* r, const uint8_t* rec);
void load_cores(rdcore* r, const uint8_t* unit, int w);
uint32_t core_cut(const rdcore* r, uint32_t i, uint64_t* out);
uint64_t core_hash(const uint64_t* w, uint32_t k);
void free_core(rdcore* r);
