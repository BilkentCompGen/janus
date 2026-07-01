#pragma once
#include "svp.h"
#include <cstdint>
#include <string>
#include <vector>

// 2-colour LCP-core graph built FROM READS (never an assembly), as 3 resumable stages
// checkpointed under a work <dir>:
//   graph  : reads -> bloom-count -> stream-fill -> self-contained <dir>/graph.bin
//            (colour + node-id edges + inline core-byte payload)
//   refine : graph.bin -> <dir>/refined.bin (topology + counts error-correction)
//   call   : refined.bin -> <dir>/calls.tsv + calls.altsigs.bin (somatic variant calls)
// Colour 1 = normal-only, 2 = tumor-only (candidate somatic), 3 = shared backbone.

// reads -> bloom-count -> stream-fill -> self-contained graph.bin.
void run_graph_reads(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int margin, int min_mult, int threads);
// Topology + counts error-correction: dir/graph.bin -> dir/refined.bin (same v2 format) (refinement.cpp).
void run_refine(const std::string& dir);
// Somatic variant calling over the refined core graph -> dir/calls.* (caller.cpp).
void run_call(const std::string& dir);
// Coordinate-reporting: map the reference-free calls onto ref.fasta (a coordinate ORACLE, never
// the tumor/normal assemblies) -> dir/calls.vcf (anchor.cpp). Explicit/opt-in; not part of `run`.
void run_anchor(const std::string& dir, const std::string& ref);
