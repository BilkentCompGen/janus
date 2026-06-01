#pragma once
#include "svp.h"
#include <cstdint>
#include <string>
#include <vector>

// 2-colour LCP-core graph built FROM READS (never an assembly), as 3 resumable stages
// checkpointed under a work <dir>:
//   lcp     : reads -> per-signature read multiplicity + within-read edges
//             -> <dir>/counts.bin, <dir>/edges.bin       (the expensive, gz-bound LCP pass)
//   graph   : counts+edges -> prune sequencing-error nodes (max-mult < min_mult), colour,
//             resolve to node-id edges -> <dir>/graph.bin
//   bubbles : graph -> permissive Cortex bubble decomposition -> <out_tsv>
// Colour 1 = normal-only, 2 = tumor-only (candidate somatic), 3 = shared backbone.
// type 'D' = normal-only, 'I' = tumor-only (somatic), 'C' = mixed.
struct CortexBubble { char type; uint32_t n_nodes; uint64_t normal_cov, tumor_cov; };

void run_lcp(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int flank, int margin, int threads);
void run_graph(const std::string& dir, int min_mult);
void run_bubbles(const std::string& dir, const std::string& out_tsv);
