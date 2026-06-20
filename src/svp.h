#pragma once
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

// SoA LCP extraction over biolcp (src/extract_fast.cpp). The reads-only caller (graph.cpp)
// calls this directly and hashes each core's ±flank window with its own emit_sig. end is
// exclusive; label is biolcp's native content fingerprint (unused by the caller — the
// signature is the murmur of the window, not the label).
extern "C" void extract_fast(const char* str, int len, uint64_t offset, int target_level,
                             int* n_out, const uint64_t** starts_out,
                             const uint64_t** ends_out, const uint32_t** labels_out);

// Shared graph-stage structs (survivor record + directed node-id edge), hoisted from
// graph.cpp's anon namespace so compact.cpp can read graph.bin without duplication.
struct SigCount { uint64_t sig; uint32_t n, t; };
using IEdge = std::pair<uint32_t, uint32_t>;

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
