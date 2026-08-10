#!/usr/bin/env bash
# `svp bubbles` drops a somatic arm when it is the same length as its germline sibling and colinear
# with it -- an SNV or a sub-SV edit. "Same length" is max_dlen MINIMIZERS, and at 12 that is roughly
# a 130-190 bp difference, far above the 50 bp floor the caller reports at. Calibrated on `hard`,
# where it cost 2 junctions; on real reads it removes 63.1% of hcc1937's somatic cores, and what
# hcc1937 misses is exactly small: its missed deletions have a median of 127 bp, its missed
# insertions 144 bp, and 60-73% of both are under 300 bp.
#
# A 50 bp indel is about 3 minimizers of length difference after homopolymer compression, so the
# threshold that separates "below the SV floor" from "an SV we are meant to report" is 2, not 12.
set -uo pipefail
cd /home/yuio/svp
S=${S:-hcc1937}
DLEN=${DLEN:-2}
G=/home/yuio/svp/outputs/svp/$S
JX=/home/yuio/tmp/jx/${S}_dlen$DLEN
REF=/home/datasets/references/genomes/hg38/hg38.fa
N=/home/yuio/tmp/jx/$S/normal2M.bam
[ "$S" = hcc1954 ] && N=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam
mkdir -p "$JX"

if [ ! -s "$JX/drop.sigs" ]; then
    echo "[$S dlen=$DLEN bubbles] $(date +%T)"
    taskset -c 0-15 ./svp bubbles "$G/graph.bin" <(./svp dump hifi/$S/n.svb) \
        <(./svp dump hifi/$S/t.svb) 16 21 22 40 $DLEN 20 2>&1 | tail -2
    cp "$G/drop.sigs" "$JX/drop.sigs"
else
    cp "$JX/drop.sigs" "$G/drop.sigs"
fi

if [ ! -s "$JX/ctg.fa" ]; then
    echo "[$S junctions] $(date +%T)"
    taskset -c 0-15 ./svp junctions "$G/graph.bin" <(./svp dump hifi/$S/t.svb) \
        "$JX/ctg.fa" 16 50 5 300 0 2>&1 | tail -2
fi
if [ ! -s "$JX/ctg.bam" ]; then
    echo "[$S align] $(date +%T)"
    taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no "$REF" "$JX/ctg.fa" 2> "$JX/mm2.log" \
        | samtools sort -@ 4 -m 3G -T "$JX/st" -o "$JX/ctg.bam" -
    samtools index -@ 8 "$JX/ctg.bam"
fi
for SP in 3 4 5; do
    ./svp call "$JX/ctg.bam" "$N" "$JX/c$SP.vcf" 10 0 50 50 $SP >/dev/null 2>&1
    .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/c$SP.vcf" \
        --min_size 50 --tolerance 500 --out_dir "$JX/m$SP" >/dev/null 2>&1
    awk -F'\t' -v s=$S -v d=$DLEN -v p=$SP 'NR>1 && $1!~/^Minda/{printf "  >> %-9s dlen=%s supp>=%s  TP %5s FN %5s FP %5s  prec %.3f  rec %.3f  F1 %.3f\n",s,d,p,$2,$3,$4,$5,$6,$7}' \
        "$JX/m$SP/results/overall.tsv"
done
echo "[$S dlen=$DLEN] bitti $(date +%T)"
