#pragma once
#include "log.h"
#include <cstdint>
#include <vector>

using cov_t = uint8_t;
constexpr cov_t COV_MAX = 255;

struct SigCount { uint64_t sig; uint32_t slot; };
struct CEdge { uint32_t src, dst; };

extern "C" void extract(const char* str, int len, int target_level,
                         int* n_out, const uint64_t** starts_out,
                         const uint64_t** ends_out, const uint64_t** sigs_out);

extern "C" void extract2(const char* str, int len, int target_level,
                          int* nf, const uint64_t** stf, const uint64_t** enf, const uint64_t** sgf,
                          int* nr, const uint64_t** str_, const uint64_t** enr, const uint64_t** sgr);
