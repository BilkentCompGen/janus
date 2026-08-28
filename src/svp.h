#pragma once
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

// SoA LCP extraction over biolcp (src/extract.cpp). Each core's signature is XXH3-64 over the
// core's own bytes (no flank), computed once here and returned in sigs_out. end is exclusive.
extern "C" void extract(const char* str, int len, uint64_t offset, int target_level,
                             int* n_out, const uint64_t** starts_out,
                             const uint64_t** ends_out, const uint64_t** sigs_out);
uint64_t core_sig(const char* key, int len);

// Shared graph-stage structs (survivor record + directed node-id edge), hoisted from
// graph.cpp's anon namespace so compact.cpp can read graph.bin without duplication.
struct SigCount { uint64_t sig; uint32_t n, t; };
using IEdge = std::pair<uint32_t, uint32_t>;

// dir/threads.bin record: a coverage-flagged repeat run was entered from raw node `in`, spans
// the run anchored at `mid`, and left to `out` (UINT32_MAX = read boundary), seen `count` times.
struct RepTriple { uint32_t mid, in, out, count; };

// Super-graph node schema: compact.cpp's build_super_graph fills these in RAM, simplify.cpp's
// run_simplify consumes them and writes the same node layout to simplify.bin — shared here so
// the build and the writer stay byte-identical. pad zero-inits the 3 alignment bytes after
// color (sizeof is 12); without it those bytes are uninitialized heap and leak to disk, making
// simplify.bin non-reproducible.
struct SuperNode { uint32_t block_len, runs_len; uint8_t color; uint8_t pad[3] = {}; };
struct CovRun    { uint32_t end_pos, n_cov, t_cov; };

// Fail-fast checkpoint IO helpers. Defined in graph.cpp; declared here so compact.cpp can use them.
FILE* xopen(const std::string& p, const char* m);
void  wr(FILE* f, const void* p, size_t n);
void  rd(FILE* f, void* p, size_t n);
void  atomic_close(FILE* f, const std::string& tmp, const std::string& dst);

// KMP longest k such that A[la-k:la] == B[0:k]. Defined in compact.cpp; shared with simplify.cpp.
uint32_t suffix_prefix(const uint8_t* A, uint32_t la, const uint8_t* B, uint32_t lb);
