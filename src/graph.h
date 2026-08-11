#pragma once
#include "snape.h"
#include <cstdint>
#include <string>
#include <vector>

// 2-colour LCP-core graph built FROM READS (never an assembly), checkpointed under a work <dir>:
//   reads -> counted cores -> self-contained <dir>/graph.bin
//   (per-node murmur sig + ncov/tcov + directed edges with per-colour counts; no sequence)
// Colour is derived from the counts on demand: ncov>0 && tcov>0 = shared, ncov>0 = normal-only,
// else tumor-only (candidate somatic).

// One pass over both samples' reads -> exact node and edge counts -> self-contained graph.bin.
void run_graph_reads(const std::string& normal_fq, const std::string& tumor_fq, const std::string& dir, int level, int margin, int min_mult, int threads);
