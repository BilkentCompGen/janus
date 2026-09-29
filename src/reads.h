#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "core.h"

struct want { std::vector<uint32_t> p; uint64_t mask, score; uint32_t pad, mul, lp, rp; uint8_t lk; std::string s; };

struct rth { rdcore r; std::vector<uint64_t> h; std::vector<uint32_t> nd, id, ix, ed; };

struct slot { uint64_t h; uint32_t v; };

struct rpass;

typedef void (*read_fn)(rpass* x, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid);

struct rpass { std::vector<rth> t; const std::vector<uint64_t>* g; std::vector<slot> tab; uint64_t mask; read_fn fn; void* ctx; int f, wb; };

struct fctx { std::vector<want>* w; std::vector<uint32_t> fo, fi; };

void each_read(rpass* x, const std::vector<uint64_t>& g, const char* reads, int t, read_fn fn, void* ctx);
void make_map(rpass* x);
void lookup(const rpass* x, rth* c);
void chains(rpass* x, rth* c, const uint8_t* rec, uint64_t rid);
void wanted(fctx* f, std::vector<want>& w, uint64_t K);
void match(const fctx* f, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, int s);
std::string bases(const rdcore* r, const uint8_t* rec, uint32_t lo, uint32_t hi, uint32_t rc);
bool same(const std::string& t, const std::string& n);
