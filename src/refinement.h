#pragma once
#include <string>

// Topology + counts error-correction over the mm-floored core graph (src/refinement.cpp):
//   <dir>/graph.bin -> <dir>/refined.bin, same v3 format.
//
// Four passes, each judging a structure by its read count RELATIVE TO LOCAL COVERAGE rather than
// against a global floor -- a global cut is what min-mult already does, and pushing it harder
// distorts real bubbles:
//   1. tips      short dead-end chains whose count is far below the junction they hang off, plus
//                low-count isolated islands
//   2. low-cov   short side-sprouts carrying a small fraction of a well-covered neighbour
//   3. bubbles   the minority arm of a short bubble -- but a balanced arm is a germline het and a
//                tumour-only arm above a VAF floor is a candidate somatic, so only a clearly
//                imbalanced non-somatic arm is popped
//   4. edges     node counts are mm-floored at build time but EDGE counts are raw, so a "tumour-only"
//                edge can rest on a single read. Applies the same floor per colour, plus a relative
//                gate against endpoint coverage.
//
// Measured on the 30x simulations the somatic chunk set shrinks ~2% (hard) / ~7.6% (complex) with no
// change in recall or precision -- the dropped chunks produced no calls either way. It is wired in
// for what pass 4 addresses at other depths: at high coverage a single-read edge is a far larger
// share of the somatic false positives than it is at 30x.
//
// Tunable through the environment; see run_refine for the full list (REF_TIP_MAX, REF_EDGE_MM, ...).
void run_refine(const std::string& dir);
