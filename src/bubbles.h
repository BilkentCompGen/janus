#pragma once

// Sub-SV bubbles (src/bubbles.cpp): the somatic cores that are an SNV or a small edit rather than a
// junction, written out as a signature list so `svp junctions` never builds a contig for them.
//
// Takes its arguments positionally, as the `svp bubbles` subcommand:
//   <graph.bin> <normal> <tumour> [threads=16] [k=21] [w=22] [maxlen=40] [max_dlen=12] [max_unmatched=20]
// Writes <graph dir>/drop.sigs. Returns 0 on success, non-zero on a usage error. LCP_INIT() must have
// run first -- it re-parses the reads and the parse must match the one that built the graph.
int cmd_bubbles(int argc, char** argv);
