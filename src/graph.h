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
// Topology-only type (from core counts): 'D'=DEL 'I'=INS 'P'=1+1-shape 'S'=N+N 'C'=other.
// Real variant call (base-pair diff) is done by tools/call from the sig lists below.
struct CortexBubble {
    char type;
    uint32_t n_nodes_n, n_nodes_t;
    uint64_t normal_cov, tumor_cov;
    std::vector<uint64_t> anchors;     // color=3 anchor sigs around the bubble
    std::vector<uint64_t> sigs_n;      // color=1 chain sigs (normal-allele cores)
    std::vector<uint64_t> sigs_t;      // color=2 chain sigs (tumor-allele cores)
};

void run_lcp(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int flank, int margin, int threads);
void run_graph(const std::string& dir, int min_mult);
void run_bubbles(const std::string& dir, const std::string& out_tsv);
