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

// Both strands of one read from ONE compression: run for run, hpc(rc(s)) is rc(hpc(s)), so the
// reverse strand is a reversal of the compressed forward string and its run map. The reverse
// coordinates are positions in the read's reverse complement, exactly as if the caller had
// reverse-complemented the read itself and called extract() on it. The graph build reads every
// base twice, so this halves the compression work; single-strand callers still use extract().
extern "C" void extract2(const char* str, int len, int target_level,
                              int* nf, const uint64_t** stf, const uint64_t** enf, const uint64_t** sgf,
                              int* nr, const uint64_t** str_, const uint64_t** enr, const uint64_t** sgr);
uint64_t core_sig(const char* key, int len);

// Per-sample coverage count. WIDTH IS A MEASURED CHOICE, not a guess.
//
// A count's MAGNITUDE is read in exactly one place that reaches an output: the chunk VAF
// (junctions.cpp). Everywhere else it is a zero test (`cov[0]==0` is somatic) or a floor test
// (min_mult), and saturation at the top can change neither. The VAF in turn feeds
// `keep[c] = vaf(c) >= min_vaf` with min_vaf defaulting to 0, so it drops nothing.
//
// MEASURED on the real hcc1954 graph (70.4 M nodes, 115,774 chunks), counts clipped at 255:
//   somatic cores 207,948 -> 207,948   chunks 115,774 -> 115,774   no-anchor 968 -> 968
// i.e. the chunk set, which is the ONLY thing that flows downstream, is unchanged. What does move
// is the VAF histogram (VAF>=0.02: 58,601 -> 66,312), because the VAF's denominator is a germline
// ANCHOR count -- local depth, not copy number -- and anchors in repeat-adjacent sequence pass 255.
// That histogram is printed and discarded.
//
// It is kept at 8 bits because the study this is built for sequences clonally expanded SINGLE
// CELLS: there is no within-line cell fraction for a VAF to estimate, and the between-line question
// -- which lines carry the junction -- is answered exactly by the support pattern, with no ratio and
// no denominator. 21 samples cost 4.2 GB of graph.bin at this width against 7.3 GB at 16 bits.
using cov_t = uint8_t;
constexpr cov_t COV_MAX = 255;

// Survivor record: the sig IS the node identity; its S counts live in a separate flat arena at
// stride S, which `slot` indexes. Putting the counts in the record (or in a hash map's value) would
// either malloc per entry or pay a fixed MAX_SAMPLES per entry; this pays 4 bytes.
struct SigCount { uint64_t sig; uint32_t slot; };
// Directed adjacency, node ids, with its own slot into the edge count arena.
struct CEdge { uint32_t src, dst; };      // no per-sample counts: nothing reads them

// One reads input, FASTQ or .svb, walked sequentially. Defined in graph.cpp, which already had the
// .svb reader, so that every stage that walks the reads uses the SAME one -- bubbles and junctions
// used to take FASTQ only, which is the entire reason the pipeline needed a shell to hand them
// `<(snape dump ...)`. A .svb record has no name, so one is synthesised from its index; that name only
// ever has to be unique per read, which is what `snape call` counts support with.
struct Reads {
    FILE* f = nullptr;
    bool svb = false;
    uint64_t idx = 0;
    std::string tag;                 // "" or "<sample>|", prepended to every name this reader hands out
};
// `sample` is the line's name from the sheet, or "" for a reader whose names nothing downstream
// distinguishes (pick's normal, bubbles' collect). With N tumour lines the tag is what stops two
// .svb files' synthesised r0, r1, ... from colliding: a contig is named after the read it was cut
// from, and `snape call` counts support by distinct contig name.
Reads reads_open(const std::string& path, const std::string& sample);
bool reads_next(Reads& r, std::string& seq, std::string* name);
void reads_close(Reads& r);

// Fail-fast checkpoint IO helpers, defined in graph.cpp: a half-written checkpoint must never be
// left where a resume would trust it.
FILE* xopen(const std::string& p, const char* m);
void  wr(FILE* f, const void* p, size_t n);
void  rd(FILE* f, void* p, size_t n);
void  atomic_close(FILE* f, const std::string& tmp, const std::string& dst);
