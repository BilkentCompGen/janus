#pragma once
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

// SoA LCP extraction over biolcp (src/extract.cpp). The sequence is homopolymer-compressed first,
// so a core's signature is MurmurHash64A over its COMPRESSED bytes -- run-length slips, which are
// what long-read indel errors overwhelmingly are, no longer change a locus's identity. starts/ends
// are still RAW coordinates (end exclusive), so callers' positions and margins are unaffected.
extern "C" void extract(const char* str, int len, uint64_t offset, int target_level,
                             int* n_out, const uint64_t** starts_out,
                             const uint64_t** ends_out, const uint64_t** sigs_out);
uint64_t core_sig(const char* key, int len);

// One kept core of one parse: its span in the read plus its signature. (Was refidx.h's, kept here now
// that the reference index is gone -- a read's cores go in as extract() produced them, both strands.)
struct Span { uint32_t st, en; uint64_t sig; };

// Shared graph-stage structs (survivor record + directed node-id edge), hoisted from
// graph.cpp's anon namespace so other stages can read graph.bin without duplication.
struct SigCount { uint64_t sig; uint32_t n, t; };
// Counted directed edge, the graph.bin v3 on-disk edge record: src->dst node ids with
// per-colour traversal counts (ncov = #normal reads crossing the adjacency, tcov = #tumor reads).
struct CEdge { uint32_t src, dst, ncov, tcov; };

// Fail-fast checkpoint IO helpers. Defined in graph.cpp; refinement reads graph.bin back through rd().
FILE* xopen(const std::string& p, const char* m);
void  wr(FILE* f, const void* p, size_t n);
void  rd(FILE* f, void* p, size_t n);
void  atomic_close(FILE* f, const std::string& tmp, const std::string& dst);
