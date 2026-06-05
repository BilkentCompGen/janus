#pragma once
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

// Shared LCP primitives. The reads-only caller (graph.cpp) uses extract_cores + compute_sigs;
// read_fasta_flat serves the assembly-based GT/analysis tools only (never a caller input).
struct Core { uint64_t start, end; uint32_t label; };

void read_fasta_flat(const std::string& path, int flank, std::string& out, std::vector<uint64_t>& rec_start, std::vector<uint64_t>& rec_end);
void read_fasta_flat_named(const std::string& path, int flank, std::string& out, std::vector<uint64_t>& rec_start, std::vector<uint64_t>& rec_end, std::vector<std::string>& names);
void extract_cores(const std::string& seq, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, int level, std::vector<Core>& out);
void compute_sigs(const std::string& seq, const std::vector<Core>& cs, int m, int context_radius, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, std::vector<uint64_t>& out);

// Shared graph-stage structs (counts.bin record + directed node-id edge), hoisted from
// graph.cpp's anon namespace so modea.cpp can read graph.bin / counts.bin without duplication.
struct SigCount { uint64_t sig; uint32_t n, t; };
using IEdge = std::pair<uint32_t, uint32_t>;

// Fail-fast checkpoint IO helpers. Defined in graph.cpp; declared here so modea.cpp can use them.
FILE* xopen(const std::string& p, const char* m);
void  wr(FILE* f, const void* p, size_t n);
void  rd(FILE* f, void* p, size_t n);
void  atomic_close(FILE* f, const std::string& tmp, const std::string& dst);
