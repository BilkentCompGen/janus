#pragma once

// Breakpoints from an aligned contig BAM (src/call.cpp) -- the stage after minimap2, in snape itself
// instead of Severus.
//
// Takes its arguments positionally, as the `snape call` subcommand:
//   <ctg.bam> <normal.bam|-> <out.vcf> [min_mapq=10] [min_aln=0] [min_sv=50] [cluster=50]
//   [min_supp=5] [shred_n=6] [shred_win=300] [max_mm=0.01] [shred_rep=0.40] [anchor_rep=0.50]
//   [normal_dirty=0.25]
// Returns 0 on success, non-zero on a usage error.
int cmd_call(int argc, char** argv);
