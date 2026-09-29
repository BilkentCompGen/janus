#pragma once
#include <cstdint>
#include <vector>

struct cnt { std::vector<uint64_t> hs, e; std::vector<uint16_t> n, ec; };

struct gph { uint64_t K, E, S; std::vector<uint64_t> g, e; std::vector<uint16_t> c, ec; std::vector<uint8_t> sh; };

struct adj { std::vector<uint32_t> oo, on, io, in; };

void read_cnt(const char* path, cnt* x);
void read_gph(const char* path, gph* x);
void write_gph(const char* path, const gph* x);
void make_adj(const gph* x, adj* a);
void masks(const gph* x, std::vector<uint64_t>& m);
std::vector<uint64_t> read_set(const char* path);
void write_set(const char* path, const char* magic, const std::vector<uint64_t>& v, uint64_t aux);
