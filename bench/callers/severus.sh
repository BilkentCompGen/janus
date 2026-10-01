#!/usr/bin/env bash
# severus.sh <tumor.bam> <normal.bam> <outdir> <threads>: Severus 1.7 (venv lib/callers/severus), unphased BAMs (as in
# our CASTLE reproduction), VNTR bed -> calls.vcf = PASS records of somatic_SVs/severus_somatic.vcf
source /home/yuio/snape/bench/callers/common.sh
$L/severus/bin/severus --target-bam $T --control-bam $N --out-dir $O/severus --vntr-bed $VNTR -t $P
pass $O/severus/somatic_SVs/severus_somatic.vcf $O/calls.vcf
done_msg
