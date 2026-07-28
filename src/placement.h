#pragma once
#include <cstdint>
#include <vector>
#include "svp.h"        // Span (via refidx.h)
#include "posmap.h"

// Lookup-based read placement (NO dynamic-programming alignment). eff = kept cores in read order;
// ref_merge is a passthrough under posaware, so eff carries raw cores with read-bp offsets in .st.
// Fills locus[i] (packed via posmap::pack_locus, or posmap::NOVEL):
//   anchor (1 reference locus)  -> that locus, PINNED, never moved. An SV read's flanking anchors stay
//                                  at their true distant loci, so the read-order edge across the
//                                  breakpoint is the junction -- the SV signal survives.
//   repeat (>=2 loci)           -> the candidate colinear with the nearest anchor (same contig+strand,
//                                  |dref-dread|<=band); else the smallest packed locus (deterministic).
//   novel  (0 loci)             -> posmap::NOVEL (kept as a node keyed by sig alone: somatic signal).
// Deterministic (no hash-order dependence). Runs ONCE, in pass 2 only.
void place(const std::vector<Span>& eff, const posmap::Map& pm, std::vector<uint64_t>& locus);
