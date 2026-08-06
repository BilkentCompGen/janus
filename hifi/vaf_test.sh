#!/usr/bin/env bash
# Does the chunk VAF gate `svp junctions` already implements explain the false positives?
#
# On hcc1954 only 18% of somatic chunks carry a tumour support above a tenth of their locus depth.
# A germline site leaking through the somatic definition would sit near 0.5, not below 0.1 -- so if
# these are errors rather than leaked germline, gating on VAF should cut the false positives and
# leave the true calls. If precision barely moves, the diagnosis is wrong and the leak is elsewhere.
set -uo pipefail
cd /home/yuio/svp
S=hcc1954
REF=/home/datasets/references/genomes/hg38/hg38.fa
G=outputs/svp/$S/graph.bin
T=hifi/$S/t.svb

for V in 0.05 0.10 0.20; do
    JX=/home/yuio/tmp/jx/${S}_vaf$V; mkdir -p "$JX"
    echo "=== min_vaf=$V  $(date +%T) ==="
    [ -s "$JX/ctg.fa" ] || taskset -c 0-15 ./svp junctions "$G" <(./svp dump "$T") "$JX/ctg.fa" \
        16 25 5 300 "$V" 2>&1 | /usr/bin/grep -E "contigs \(|chunks:"
    [ -s "$JX/ctg.bam" ] || { taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no \
        -R "@RG\tID:t\tSM:t" "$REF" "$JX/ctg.fa" 2> "$JX/mm2.log" \
        | samtools sort -@ 4 -m 3G -T "$JX/st" -o "$JX/ctg.bam" - && samtools index -@ 8 "$JX/ctg.bam"; }
    ./svp call "$JX/ctg.bam" "$JX/calls.vcf" 2>&1 | tail -1
    .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/calls.vcf" \
        --min_size 50 --tolerance 500 --out_dir "$JX/minda" >/dev/null 2>&1
    awk -F'\t' 'NR>1 && $1!~/^Minda/{printf "  min_vaf='"$V"'  TP %s  FN %s  FP %s  prec %.3f  rec %.3f  F1 %.3f\n",$2,$3,$4,$5,$6,$7}' "$JX/minda/results/overall.tsv"
done
echo "=== bitti $(date +%T) ==="
