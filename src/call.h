#pragma once

// Breakpoints from aligned contig BAMs (src/call.cpp) -- the stage after minimap2, in snape itself
// instead of Severus.
//
// Takes its arguments positionally, as the `snape call` subcommand:
//   <ctg.bam|name=bam,name=bam,...> <normal.bam|-> <out.vcf> [min_mapq=10] [min_aln=0] [min_sv=50]
//   [cluster=50] [min_supp=5] [shred_n=6] [shred_win=300] [max_mm=0.01] [shred_rep=0.40]
//   [anchor_rep=0.50] [normal_dirty=0.25] [supp_assign=3]
//
// The first argument is one contig set per tumour line. A bare path is one line named "t", which is
// the sample column this VCF has always carried, so one line runs the same code and prints what it
// always printed. `min_supp` is the DISCOVERY floor and is per line -- five contigs from one line,
// not five spread over five; `supp_assign` is the lower floor at which a line is then called
// present. With one line assignment cannot bite, because discovery already demanded more.
//
// Writes <out.vcf> with one sample column per line, plus <out>.<name>.vcf projections beside it.
// Returns 0 on success, non-zero on a usage error.
int cmd_call(int argc, char** argv);
