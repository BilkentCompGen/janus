#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

struct bin_header { char magic[4]; uint32_t ver; uint64_t nreads, total; uint8_t resv[8]; };

struct code_table { uint8_t t[256]; code_table(); };

extern const code_table CODE;

struct block { int id; const char* in[2]; size_t inlen[2]; std::vector<char> out; void* ctx; uint64_t rec; };

struct bits { std::vector<char>* out; uint64_t acc; int pos; };

typedef size_t (*size_fn)(const char* p, size_t n);
typedef void (*work_fn)(block* b);
typedef void (*par_fn)(void* ctx, int i);

FILE* open_in(const char* path, bin_header* h);
FILE* open_tmp(const char* path);
void finish(FILE* fo, const char* path, const bin_header* h);
void par(int n, par_fn f, void* ctx);
uint64_t be64(const uint8_t* p);
uint32_t rbits(const uint8_t* p, uint64_t* pos, int w);
void wbits(bits* s, uint64_t v, int w);
void wflush(bits* s);
size_t nlines(const char* p, size_t n, int k);
size_t svb_size(const char* p, size_t n);
std::vector<std::string> split(const char* s);
uint64_t stream_run(FILE* fi, FILE* f2, FILE* fo, int t, size_fn sf, size_fn s2, work_fn wf, void* ctx);
