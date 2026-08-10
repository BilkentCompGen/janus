#!/usr/bin/env bash
# Half of what we miss on hcc1954 is a PLACEMENT failure, not a detection failure: 34% of the misses
# have a contig at one end whose SA points to a third locus, and 16% have contigs at both ends that
# are never the same contig. Both are what more unique anchor sequence fixes. pad_cores is that
# anchor, and raising it cannot cost a contig -- the window is clamped to the read's own cores and
# min_flank tests the flank ACHIEVED, so a bigger pad only ever takes more of what the read has.
#
# Everything else is held: same graph.bin, same drop.sigs, per_chunk 5, and `svp call` at its
# committed defaults.
set -uo pipefail
cd /home/yuio/svp
S=${S:-hcc1954}
G=/home/yuio/svp/outputs/svp/$S
REF=/home/datasets/references/genomes/hg38/hg38.fa
N=/home/yuio/tmp/jx/$S/normal2M.bam
[ "$S" = hcc1954 ] && N=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam

for PAD in "$@"; do
    JX=/home/yuio/tmp/jx/${S}_pad$PAD
    mkdir -p "$JX"
    echo "########## $S pad=$PAD  $(date +%T) ##########"
    if [ ! -s "$JX/ctg.fa" ]; then
        taskset -c 0-15 ./svp junctions "$G/graph.bin" <(./svp dump hifi/$S/t.svb) \
            "$JX/ctg.fa" 16 $PAD 5 300 0 2>&1 | tail -2
    fi
    if [ ! -s "$JX/ctg.bam" ]; then
        echo "  [align] $(date +%T)"
        taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no "$REF" "$JX/ctg.fa" 2> "$JX/mm2.log" \
            | samtools sort -@ 4 -m 3G -T "$JX/st" -o "$JX/ctg.bam" -
        samtools index -@ 8 "$JX/ctg.bam"
    fi
    for SP in 3 4 5; do
        ./svp call "$JX/ctg.bam" "$N" "$JX/c$SP.vcf" 10 0 50 50 $SP >/dev/null 2>&1
        .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/c$SP.vcf" \
            --min_size 50 --tolerance 500 --out_dir "$JX/m$SP" >/dev/null 2>&1
        awk -F'\t' -v p=$PAD -v s=$SP 'NR>1 && $1!~/^Minda/{printf "  >> pad=%-3s supp>=%s  TP %5s FN %5s FP %5s  prec %.3f  rec %.3f  F1 %.3f\n",p,s,$2,$3,$4,$5,$6,$7}' \
            "$JX/m$SP/results/overall.tsv"
    done
done
echo "########## pad sweep bitti $(date +%T) ##########"
