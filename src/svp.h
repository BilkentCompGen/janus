#pragma once
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

struct Core   { uint64_t start, end; uint32_t label; };
struct Anchor { uint32_t n_idx, t_idx; };
struct Bubble { uint64_t n_start, n_end, t_start, t_end; std::string_view n_seq, t_seq; const char* type; };

void read_fasta_flat(const std::string& path, int flank, std::string& out, std::vector<uint64_t>& rec_start, std::vector<uint64_t>& rec_end);
void extract_cores(const std::string& seq, const std::vector<uint64_t>& rec_start, const std::vector<uint64_t>& rec_end, int level, std::vector<Core>& out);
void compute_sigs(const std::string& seq, const std::vector<Core>& cs, int m, std::vector<uint64_t>& out);
void unique_anchors(const std::vector<uint64_t>& nsig, const std::vector<uint64_t>& tsig, std::vector<Anchor>& out);
void chain(const std::vector<Anchor>& anchors, std::vector<int>& kept);
void emit_bubbles(const std::string& N, const std::vector<Core>& ncs, const std::string& T, const std::vector<Core>& tcs, const std::vector<Anchor>& anchors, const std::vector<int>& kept, std::vector<Bubble>& out);
void write_bubbles(const std::string& path, const std::vector<Bubble>& bubs);
