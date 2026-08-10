#!/usr/bin/env bash
# `svp junctions` emits at most per_chunk contigs per chunk, so at 5 a support floor of 5 means
# "every contig there agrees" and one ragged contig kills the call -- which is why raising min_supp
# buys precision only by giving up recall. This raises the CAP instead, so support counts evidence
# rather than a fraction of five, and re-sweeps the floor on top of it.
set -uo pipefail
cd /home/yuio/svp
S=hcc1954
G=/home/yuio/svp/outputs/svp/$S
JX=/home/yuio/tmp/jx/${S}_pc10
N=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam
REF=/home/datasets/references/genomes/hg38/hg38.fa
mkdir -p "$JX"

if [ ! -s "$JX/ctg.fa" ]; then
    echo "[junctions per_chunk=10] $(date +%T)"
    taskset -c 0-15 ./svp junctions "$G/graph.bin" <(./svp dump hifi/$S/t.svb) "$JX/ctg.fa" 16 25 10 300 0 2>&1 | tail -2
fi
if [ ! -s "$JX/ctg.bam" ]; then
    echo "[align] $(date +%T)"
    taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no "$REF" "$JX/ctg.fa" 2> "$JX/mm2.log" \
        | samtools sort -@ 4 -m 3G -T "$JX/st" -o "$JX/ctg.bam" -
    samtools index -@ 8 "$JX/ctg.bam"
fi

for SP in 3 4 5 6 8 10; do
    ./svp call "$JX/ctg.bam" "$N" "$JX/c$SP.vcf" 10 0 50 50 $SP 6 300 0.01 0.40 >/dev/null 2>&1
    .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/c$SP.vcf" \
        --min_size 50 --tolerance 500 --out_dir "$JX/m$SP" >/dev/null 2>&1
    awk -F'\t' -v p=$SP 'NR>1 && $1!~/^Minda/{printf "  per_chunk=10 supp>=%-3s TP %5s FP %5s  prec %.3f  rec %.3f  F1 %.3f\n",p,$2,$4,$5,$6,$7}' \
        "$JX/m$SP/results/overall.tsv"
done
echo "bitti $(date +%T)"
