#pragma once

// svp pick <ctg.fa> <calls.vcf> <normal reads> <out.fastq> [threads] [min_shared]
//
// The normal reads worth aligning: the ones sharing LCP cores with the GERMLINE ANCHORS of the
// contigs the candidate calls name. `svp call` only ever asks the normal about those loci, so
// aligning the whole sample to answer it is work the graph can avoid.
int cmd_pick(int argc, char** argv);
