#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

static const uint32_t NONE = ~0u; // no node, no edge, no slot
static const uint32_t PF = 16; // prefetch distance of loops with one random access per step

struct bin_header { char magic[4]; uint32_t ver; uint64_t nreads, total; uint8_t resv[8]; };

struct block { int id; const char* in[2]; size_t inlen[2]; std::vector<char> out; void* ctx; uint64_t rec; };

struct bits { std::vector<char>* out; uint64_t acc; int pos; };

typedef void (*work_fn)(block* b);
typedef void (*par_fn)(void* ctx, int i);

bin_header header(char kind, uint32_t ver, uint64_t n, uint64_t total);
FILE* open_in(const char* path, bin_header* h);
FILE* open_tmp(const char* path);
void finish(FILE* fo, const char* path, const bin_header* h);
void save(const char* path, const bin_header* h, const void* p, size_t n);
void par(int n, par_fn f, void* ctx);
uint8_t code(int c);
uint64_t be64(const uint8_t* p);
uint32_t rbits(const uint8_t* p, uint64_t* pos, int w);
void wbits(bits* s, uint64_t v, int w);
void wflush(bits* s);
size_t svb_bytes(uint32_t tag);
size_t svb_size(const char* p, size_t n);
size_t exc_size(const char* p, size_t n);
std::vector<std::string> split(const char* s);
uint64_t stream_run(FILE* fi, FILE* fe, FILE* fo, int t, work_fn wf, void* ctx);
