#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdio>

struct bin_header { char magic[4]; uint32_t ver; uint64_t nreads, total; uint8_t resv[8]; };

struct code_table { uint8_t t[256]; code_table(); };

extern const code_table CODE;

struct block { int id; const char* in; size_t inlen; const char* in2; size_t inlen2; char* out; size_t outlen, outcap; char* buf; size_t bufcap; void* ctx; };

struct bits { uint64_t acc; int pos; };

typedef size_t (*size_fn)(const char* p, size_t n);
typedef void (*work_fn)(block* b);

const uint8_t* map_file(const char* path);
FILE* open_tmp(const char* path);
void commit_file(FILE* fo, const char* path);
void bput(block* b, const void* d, size_t n);
void bgrow(block* b, size_t n);
void wbits(block* b, bits* s, uint64_t v, int w);
void wflush(block* b, bits* s);
size_t nlines(const char* p, size_t n, int k);
size_t svb_size(const char* p, size_t n);
uint64_t stream_run(FILE* fi, FILE* f2, FILE* fo, size_t hdrsz, int t, size_fn sf, size_fn s2, work_fn wf, void* ctx);
