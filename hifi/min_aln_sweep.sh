#!/usr/bin/env bash
# min_aln is a segment gate, and a BND pays it TWICE: a deletion needs one surviving alignment of a
# contig, a breakend needs two. Measured on hcc1954's 246 missed BNDs, 48.4% have >=3 alignments at
# BOTH breakpoints at MAPQ>=10 and lose them only to the 1000 bp floor. This asks what the floor is
# actually worth, with everything else held at the committed defaults and ONE frozen binary.
set -uo pipefail
cd /home/yuio/svp
SVP=/home/yuio/tmp/svp-d34b6f1
S=hcc1954
JX=/home/yuio/tmp/jx/$S
N=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam
OUT=/home/yuio/tmp/minaln
mkdir -p $OUT

for A in 0 50 100 150; do
    echo "########## min_aln=$A  $(date +%T) ##########"
    $SVP call "$JX/ctg.bam" "$N" "$OUT/a$A.vcf" 10 $A 2>&1 | /usr/bin/grep -E "aday|elenen"
    .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$OUT/a$A.vcf" \
        --min_size 50 --tolerance 500 --out_dir "$OUT/m$A" >/dev/null 2>&1
    awk -F'\t' -v a=$A 'NR>1 && $1!~/^Minda/{printf "  >> min_aln=%-5s TP %5s  FN %5s  FP %6s  prec %.3f  rec %.3f  F1 %.3f\n",a,$2,$3,$4,$5,$6,$7}' \
        "$OUT/m$A/results/overall.tsv"
done
echo "########## bitti $(date +%T) ##########"
