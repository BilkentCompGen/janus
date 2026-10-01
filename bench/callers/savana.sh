#!/usr/bin/env bash
# savana.sh <tumor.bam> <normal.bam> <outdir> <threads>: SAVANA 1.3.8 (venv lib/callers/savana), default run + classify
# with the PacBio thresholds (--pb) -> calls.vcf = PASS records of <sample>.classified.somatic.vcf
source /home/yuio/snape/bench/callers/common.sh
$L/savana/bin/savana --tumour $T --normal $N --ref $REF --outdir $O/savana --threads $P --pb --sample s --overwrite
pass $O/savana/s.classified.somatic.vcf $O/calls.vcf
done_msg
