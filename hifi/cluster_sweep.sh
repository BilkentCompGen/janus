#!/usr/bin/env bash
# Two junctions are one event when both sides agree within `cluster` bp, and an event needs
# `min_supp` distinct contigs. `svp junctions` emits at most 5 contigs per chunk, so 3 of 5 must
# agree on BOTH breakpoints within 50 bp -- and a real cancer breakpoint carries microhomology, which
# is exactly what makes an aligner place the same junction a few bases apart in different reads.
# Same contig BAM and one binary throughout, so the only thing moving is the pair below.
set -uo pipefail
cd /home/yuio/svp
JX=/home/yuio/tmp/jx/hcc1954
N=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam
OUT=/home/yuio/tmp/clsweep
mkdir -p $OUT

for CL in 50 100 200 400; do
    for SUPP in 2 3; do
        ./svp call "$JX/ctg.bam" "$N" "$OUT/c${CL}_s${SUPP}.vcf" 10 0 50 $CL $SUPP >/dev/null 2>&1
        .venv/bin/minda truthset --base hifi/hcc1954/truth_clean.vcf --vcfs "$OUT/c${CL}_s${SUPP}.vcf" \
            --min_size 50 --tolerance 500 --out_dir "$OUT/m${CL}_${SUPP}" >/dev/null 2>&1
        awk -F'\t' -v c=$CL -v s=$SUPP 'NR>1 && $1!~/^Minda/{printf "  cluster=%-4s supp>=%s  TP %5s  FN %5s  FP %6s  prec %.3f  rec %.3f  F1 %.3f\n",c,s,$2,$3,$4,$5,$6,$7}' \
            "$OUT/m${CL}_${SUPP}/results/overall.tsv"
    done
done
echo "bitti $(date +%T)"
