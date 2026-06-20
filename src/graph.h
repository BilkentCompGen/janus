#pragma once
#include "svp.h"
#include <cstdint>
#include <string>
#include <vector>

// 2-colour LCP-core graph built FROM READS (never an assembly), as 2 resumable stages
// checkpointed under a work <dir>:
//   graph    : reads -> bloom-count -> stream-fill -> self-contained <dir>/graph.bin
//              (colour + node-id edges + inline ±flank window payload)
//   simplify : graph -> in-RAM unitig compaction (compact.cpp) then bubble-zip to
//              base-pair resolution -> <dir>/simplify.bin  (the compacted graph is never
//              written to disk)
// Colour 1 = normal-only, 2 = tumor-only (candidate somatic), 3 = shared backbone.

// merged: reads -> bloom-count -> stream-fill -> self-contained graph.bin.
void run_graph_reads(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int flank, int margin, int min_mult, int threads);
// Unitig compaction, built straight into RAM for run_simplify (defined in compact.cpp).
void build_super_graph(const std::string& dir, std::vector<SuperNode>& snodes,
                       std::vector<uint8_t>& blocks, std::vector<CovRun>& runs, std::vector<IEdge>& super_edges);
void run_simplify(const std::string& dir);
