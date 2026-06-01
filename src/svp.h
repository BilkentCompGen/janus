#pragma once
#include <cstdint>
#include <string>
#include <vector>

// Shared LCP primitives. The reads-only caller (graph.cpp) uses extract_cores + compute_sigs;
// read_fasta_flat serves the assembly-based GT/analysis tools only (never a caller input).
struct Core { uint64_t start, end; uint32_t label; };

void read_fasta_flat(const std::string& path, int flank, std::string& out, std::vector<uint64_t>& rec_start, std::vector<uint64_t>& rec_end);
void extract_cores(const std::string& seq, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, int level, std::vector<Core>& out);
void compute_sigs(const std::string& seq, const std::vector<Core>& cs, int m, int context_radius, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, std::vector<uint64_t>& out);
