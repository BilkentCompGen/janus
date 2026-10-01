#!/usr/bin/env bash
# sniffles2.sh <tumor.bam> <normal.bam> <outdir> <threads>: Sniffles2 2.8.0 (venv lib/callers/sniffles) as the CASTLE
# panel ran it: one .snf per sample (tandem-repeat bed), joint calling normal + tumor, then the PASS calls supported by
# the tumor only (SUPP_VEC=01) -> calls.vcf
source /home/yuio/snape/bench/callers/common.sh
S=$L/sniffles/bin/sniffles
$S --input $N --snf $O/n.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
$S --input $T --snf $O/t.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
$S --input $O/n.snf $O/t.snf --vcf $O/joint.vcf --threads $P --reference $REF --allow-overwrite
awk -F'\t' '/^#/ || ($7=="PASS" && $8 ~ /SUPP_VEC=01(;|$)/)' $O/joint.vcf > $O/calls.vcf
done_msg
