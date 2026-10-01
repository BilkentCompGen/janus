#!/usr/bin/env bash
# nanomonsv.sh <tumor.bam> <normal.bam> <outdir> <threads>: nanomonsv 0.9.0 (venv lib/callers/nanomonsv; racon and
# minimap2 from ~/.local/bin), parse both samples, get with the PacBio HiFi preset (--qv25) and the VNTR bed as the
# simple-repeat bed -> calls.vcf = PASS records of t.nanomonsv.result.vcf
source /home/yuio/snape/bench/callers/common.sh
M=$L/nanomonsv/bin/nanomonsv
sort -k1,1 -k2,2n $VNTR | bgzip > $O/repeat.bed.gz && tabix -p bed $O/repeat.bed.gz
mkdir -p $O/t $O/n
$M parse $T $O/t/t --reference_fasta $REF
$M parse $N $O/n/n --reference_fasta $REF
$M get $O/t/t $T $REF --control_prefix $O/n/n --control_bam $N --processes $P --qv25 --simple_repeat_bed $O/repeat.bed.gz
pass $O/t/t.nanomonsv.result.vcf $O/calls.vcf
done_msg
