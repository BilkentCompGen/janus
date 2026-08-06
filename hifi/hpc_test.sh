#!/usr/bin/env bash
# Does collapsing tandem repeats beyond homopolymers fix the real-data false positives?
#
# Measured on hcc1954: 41.7% of false-positive loci sit in a tandem repeat against 0.4% of true ones,
# and the periods are 2 (31%) and 5 (37%) -- period 1 is 0%, which is the homopolymer collapse already
# doing its job. So the question is whether widening that collapse removes the noise, or whether it
# starts merging distinct loci and makes things worse. Node count is the tell: if it drops far more
# than the sequence does, separate loci are sharing signatures.
#
# hcc1937, the cheapest pair with a truthset. SVP_HPC_PERIOD has to be set for graph, bubbles AND
# junctions -- all three parse through extract() and a mismatch between them is silent corruption.
set -uo pipefail
cd /home/yuio/svp
REF=/home/datasets/references/genomes/hg38/hg38.fa
S=hcc1937
N=hifi/$S/n.svb T=hifi/$S/t.svb

for P in 2 6 20; do
    export SVP_HPC_PERIOD=$P
    OUT=outputs/svp/${S}_p$P JX=/home/yuio/tmp/jx/${S}_p$P
    mkdir -p "$OUT" "$JX"
    echo "########## SVP_HPC_PERIOD=$P  $(date +%T) ##########"

    [ -s "$OUT/graph.bin" ] || /usr/bin/time -v -o "$JX/t_graph.txt" taskset -c 0-15 ./svp graph \
        "$N" "$T" "$OUT" --threads 16 --min-mult 3 2>&1 | /usr/bin/grep -E "nodes|edges|survivor" | tail -3
    [ -s "$OUT/drop.sigs" ] || taskset -c 0-15 ./svp bubbles "$OUT/graph.bin" \
        <(./svp dump "$N") <(./svp dump "$T") 16 21 22 40 12 20 2>&1 | tail -2
    [ -s "$JX/ctg.fa" ] || taskset -c 0-15 ./svp junctions "$OUT/graph.bin" <(./svp dump "$T") \
        "$JX/ctg.fa" 16 25 5 300 0 2>&1 | /usr/bin/grep -E "somatic ->|contigs \("
    [ -s "$JX/ctg.bam" ] || { taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no \
        -R "@RG\tID:t\tSM:t" "$REF" "$JX/ctg.fa" 2> "$JX/mm2.log" \
        | samtools sort -@ 4 -m 3G -T "$JX/st" -o "$JX/ctg.bam" - && samtools index -@ 8 "$JX/ctg.bam"; }
    ./svp call "$JX/ctg.bam" "$JX/calls.vcf" 2>&1 | tail -1
    .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/calls.vcf" \
        --min_size 50 --tolerance 500 --out_dir "$JX/minda" >/dev/null 2>&1
    awk -F'\t' 'NR>1 && $1!~/^Minda/{printf "  P='"$P"'  TP %s  FN %s  FP %s  prec %.3f  rec %.3f  F1 %.3f\n",$2,$3,$4,$5,$6,$7}' "$JX/minda/results/overall.tsv"
    ls -l --block-size=1 "$OUT/graph.bin" | awk '{printf "  graph.bin %.2f GB\n",$5/1e9}'
done
echo "########## bitti $(date +%T) ##########"
