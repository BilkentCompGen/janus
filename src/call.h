#pragma once

// Breakpoints from an aligned contig BAM (src/call.cpp) -- the stage after minimap2, in svp itself
// instead of Severus.
//
// Takes its arguments positionally, as the `svp call` subcommand:
//   <ctg.bam> <out.vcf> [min_mapq=10] [min_aln=1000] [min_sv=50] [cluster=50] [min_supp=3]
// Returns 0 on success, non-zero on a usage error.
int cmd_call(int argc, char** argv);
