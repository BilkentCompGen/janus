#pragma once
#include <cstdint>
#include <vector>
#include "files.h"
#include "reads.h"

struct hit { uint64_t key; uint32_t a, b; };

struct keys { uint64_t S, Z; const gph* g; const adj* a; std::vector<uint32_t> seg; std::vector<uint8_t> br, ks, on; std::vector<uint64_t> key, best; std::vector<uint32_t> n; };

void make_keys(keys* x, const gph* g, const adj* a);
uint64_t keep_keys(keys* x, const std::vector<uint8_t>& sd);
void make_table(keys* x, uint64_t np);
bool kept(const keys* x, uint64_t key);
void chain_keys(const keys* x, rth* c, uint32_t m, std::vector<hit>& out);
void tally(keys* x, rth* c, uint32_t m, uint64_t s, uint64_t rid, std::vector<hit>& h);
