#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
S=$L/sniffles/bin/sniffles
$S --input $N --snf $O/n.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
$S --input $T --snf $O/t.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
$S --input $O/n.snf $O/t.snf --vcf $O/joint.vcf --threads $P --reference $REF --allow-overwrite
awk -F'\t' '/^#/ || ($7=="PASS" && $8 ~ /SUPP_VEC=01(;|$)/)' $O/joint.vcf > $O/calls.vcf
done_msg
