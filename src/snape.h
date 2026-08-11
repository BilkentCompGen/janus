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

// One kept core of one parse: its span in the read plus its signature. A read's cores go in as
// extract() produced them, both strands.
struct Span { uint32_t st, en; uint64_t sig; };

// Shared graph-stage structs (survivor record + directed node-id edge), hoisted from
// graph.cpp's anon namespace so other stages can read graph.bin without duplication.
struct SigCount { uint64_t sig; uint32_t n, t; };
// Counted directed edge, the graph.bin on-disk edge record: src->dst node ids with
// per-colour traversal counts (ncov = #normal reads crossing the adjacency, tcov = #tumor reads).
struct CEdge { uint32_t src, dst, ncov, tcov; };

// One reads input, FASTQ or .svb, walked sequentially. Defined in graph.cpp, which already had the
// .svb reader, so that every stage that walks the reads uses the SAME one -- bubbles and junctions
// used to take FASTQ only, which is the entire reason the pipeline needed a shell to hand them
// `<(snape dump ...)`. A .svb record has no name, so one is synthesised from its index; that name only
// ever has to be unique per read, which is what `snape call` counts support with.
struct Reads {
    FILE* f = nullptr;
    bool svb = false;
    uint64_t idx = 0;
};
Reads reads_open(const std::string& path);
bool reads_next(Reads& r, std::string& seq, std::string* name);
void reads_close(Reads& r);

// Fail-fast checkpoint IO helpers, defined in graph.cpp: a half-written checkpoint must never be
// left where a resume would trust it.
FILE* xopen(const std::string& p, const char* m);
void  wr(FILE* f, const void* p, size_t n);
void  rd(FILE* f, void* p, size_t n);
void  atomic_close(FILE* f, const std::string& tmp, const std::string& dst);
