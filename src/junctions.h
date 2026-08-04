#pragma once

// Junction contigs: the somatic material of a graph, carved out of the reads that carry it
// (src/junctions.cpp).
//
// A CHUNK is a connected component of the graph's somatic subgraph (ncov == 0, tcov > 0) -- one
// junction candidate, found without a reference. For each chunk the reads are scanned once and the
// first `per_chunk` reads carrying it are cut down to the novel span plus `pad_cores` of anchor
// either side, which is what an aligner needs to place the breakpoint and nothing more.
//
// Takes its arguments positionally, as the `svp junctions` subcommand:
//   <graph.bin> <reads.fastq> <out.fa> [threads=16] [pad_cores=10] [per_chunk=5] [min_flank=300]
//   [min_vaf=0]
// Returns 0 on success, non-zero on a usage error. LCP_INIT() must have run first -- it re-parses the
// reads and the parse must match the one that built the graph.
int cmd_junctions(int argc, char** argv);
