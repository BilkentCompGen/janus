#!/bin/bash
# sniffles_somatic.sh <pair>: Sniffles2 as the CASTLE panel ran it: one .snf per sample, joint calling of normal + tumor,
# then the PASS calls supported by the tumor only (SUPP_VEC=01) -> sniffles/<pair>/sniffles.vcf
S=/home/yuio/snape; R=$S/ref; H=$S/hifi/$1; L=$S/lib/callers/sniffles/bin/sniffles; D=$S/bench/callers/sniffles/$1; T=16
set -e
for s in n t; do $L --input $H/$s.bam --snf $D/$s.snf --threads $T --reference $R/hg38_parmask.fa --tandem-repeats $R/vntr.trf.bed --allow-overwrite; done
$L --input $D/n.snf $D/t.snf --vcf $D/joint.vcf --threads $T --reference $R/hg38_parmask.fa --allow-overwrite
awk -F'\t' '/^#/ || ($7=="PASS" && $8 ~ /SUPP_VEC=01(;|$)/)' $D/joint.vcf > $D/sniffles.vcf
