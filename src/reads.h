#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "core.h"

struct rth { rdcore r; std::vector<uint64_t> h; std::vector<uint32_t> nd, id, ix, ed; };

struct slot { uint64_t h; uint32_t v; };

struct rpass;

typedef void (*read_fn)(rpass* x, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid);

struct rpass { std::vector<rth> t; const std::vector<uint64_t>* g; std::vector<slot> tab; uint64_t mask; read_fn fn; void* ctx; int f, wb; };

void each_read(rpass* x, const std::vector<uint64_t>& g, const char* reads, int t, read_fn fn, void* ctx);
void make_map(rpass* x);
void lookup(const rpass* x, rth* c);
void chains(rpass* x, rth* c, const uint8_t* rec, uint64_t rid);
std::string bases(const rdcore* r, const uint8_t* rec, uint32_t lo, uint32_t hi, uint32_t rc);
bool same(const std::string& t, const std::string& n);
