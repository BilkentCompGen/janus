#pragma once
#include "samples.h"
#include "snape.h"
#include <cstdint>
#include <string>
#include <vector>

// (N+1)-colour LCP-core graph built FROM READS (never an assembly), checkpointed under a work <dir>:
//   reads -> counted cores -> self-contained <dir>/graph.bin
//   (per-node murmur sig + one count plane per sample + directed edges with the same planes)
//
// Colour is derived from the counts on demand and never stored: plane 0 is the normal, so cov[0]==0
// is somatic, and the SUPPORT PATTERN of a somatic core -- which tumour lines carry it -- is which
// of planes 1..S-1 are non-zero. With one tumour that reduces to the old shared / normal-only /
// tumor-only trichotomy exactly.
void run_graph_reads(const std::vector<Sample>& samples, const std::string& dir,
                     int level, int margin, int min_mult, int threads);
