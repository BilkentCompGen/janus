#pragma once
#include <cstdint>
#include <vector>
#include "files.h"
#include "reads.h"

struct hit { uint64_t key; uint32_t a, b; };

struct keys { uint64_t S, Z; const gph* g; const adj* a; std::vector<uint32_t> unit, off, len; std::vector<uint8_t> br, ku, kb; std::vector<uint64_t> k, best; std::vector<uint32_t> n; };

void make_keys(keys* x, const gph* g, const adj* a, const std::vector<uint8_t>& drop);
void chain_keys(const keys* x, rth* c, uint32_t m, std::vector<hit>& out);
void count_keys(keys* x, const char* reads, int t);
