#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct cnt { std::vector<uint64_t> hs, e; std::vector<uint16_t> n, ec; };

struct gph { uint64_t K, E, S; std::vector<uint64_t> g, e; std::vector<uint16_t> c; std::vector<uint8_t> nv; std::vector<uint32_t> st; };

struct adj { std::vector<uint32_t> oo, on, io, in; };

struct seqs { std::vector<std::string> name, seq; };

void read_cnt(const char* path, cnt* x);
void read_gph(const char* path, gph* x);
void write_gph(const char* path, const gph* x);
void make_adj(const gph* x, adj* a);
void read_fasta(const char* path, seqs* f);
