#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

struct block;

struct rdcore { std::vector<uint64_t> hs, rs, mk; std::vector<uint32_t> cum, a, b; std::vector<uint8_t> raw; uint32_t len, k, pad, n, nf; };

void collapse(rdcore* r, const uint8_t* rec);
void load_cores(rdcore* r, const uint8_t* u, int w);
const uint8_t* next_read(rdcore* r, const block* b, size_t* o, size_t* e, int w);
uint32_t core_cut(const rdcore* r, uint32_t i, uint64_t* out);
uint64_t core_hash(const uint64_t* w, uint32_t k);
