#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "core.h"

struct want { std::vector<uint32_t> p; uint64_t mask, score; uint32_t pad, mul, lp, rp; uint8_t lk; std::string s; };

struct rth { rdcore r; std::vector<uint64_t> h; std::vector<uint32_t> nd, id, ix, ed; };

struct rpass;

typedef void (*read_fn)(rpass* x, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid);

struct rpass { std::vector<rth> t; const std::vector<uint64_t>* g; std::vector<uint32_t> off; read_fn fn; void* ctx; int f, wb, ids; };

void each_read(rpass* x, const std::vector<uint64_t>& g, const char* reads, const char* ids, int t, read_fn fn, void* ctx);
void fetch(std::vector<want>& w, const std::vector<uint64_t>& g, const char* reads, int t);
std::string bases(const rdcore* r, const uint8_t* rec, uint32_t lo, uint32_t hi, uint32_t rc);
bool same(const std::string& t, const std::string& n);
