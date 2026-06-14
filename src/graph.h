#pragma once
#include "svp.h"
#include <cstdint>
#include <string>
#include <vector>

// 2-colour LCP-core graph built FROM READS (never an assembly), as 5 resumable stages
// checkpointed under a work <dir>:
//   lcp     : reads -> per-sig multiplicity + within-read directed edges, prune by min_mult,
//             second LCP pass fills window bytes per kept sig
//             -> <dir>/counts.bin (post-prune sigs + inline bytes), <dir>/edges.bin
//   graph   : counts+edges -> colour, resolve to directed node-id edges -> <dir>/graph.bin
//   modeb   : graph+counts -> in-memory unitig compaction (Mode A) then bubble-zip to
//             base-pair resolution -> <dir>/modeb.bin  (Mode A is never written to disk)
//   bubbles : graph -> permissive Cortex bubble decomposition -> <out_tsv>
// Colour 1 = normal-only, 2 = tumor-only (candidate somatic), 3 = shared backbone.
// type 'D' = normal-only, 'I' = tumor-only (somatic), 'C' = mixed.
struct CortexBubble { char type; uint32_t n_nodes; uint64_t normal_cov, tumor_cov; };

void run_lcp(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int flank, int margin, int min_mult, int threads);
void run_graph(const std::string& dir);
// merged: reads -> bloom-count -> stream-fill -> self-contained graph.bin (no counts.bin/edges.bin).
void run_graph_reads(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int flank, int margin, int min_mult, int threads);
// Mode A unitig compaction, built straight into RAM for run_modeb (defined in modea.cpp).
void build_super_graph(const std::string& dir, std::vector<SuperNode>& snodes,
                       std::vector<uint8_t>& blocks, std::vector<CovRun>& runs, std::vector<IEdge>& super_edges);
void run_modeb(const std::string& dir);
void run_bubbles(const std::string& dir, const std::string& out_tsv);
