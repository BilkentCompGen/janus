#!/usr/bin/env bash
# Does carrying the graph's own read counts into the caller fix what threshold tuning could not?
#
# `svp junctions` emits at most per_chunk contigs per chunk, so the caller's support statistic --
# distinct contig names -- measures that cap and not the evidence: a junction carried by 60 tumour
# reads and one carried by 3 both arrive as 5 records. The chunk's real read count and the locus depth
# now travel in the contig name, and the caller gates on those instead.
#
# hcc1937 and hcc1954: the cheapest pair and the one with the most truth. Contigs have to be rebuilt
# because the names changed, so this is junctions -> align -> call, not just call.
set -uo pipefail
cd /home/yuio/svp
REF=/home/datasets/references/genomes/hg38/hg38.fa

for S in hcc1937 hcc1954; do
    JX=/home/yuio/tmp/jx/${S}_depth; mkdir -p "$JX"
    echo "########## $S  $(date +%T) ##########"
    [ -s "$JX/ctg.fa" ] || taskset -c 0-15 ./svp junctions outputs/svp/$S/graph.bin \
        <(./svp dump hifi/$S/t.svb) "$JX/ctg.fa" 16 25 5 300 0 2>&1 | /usr/bin/grep -E "contigs \("
    head -1 "$JX/ctg.fa"
    [ -s "$JX/ctg.bam" ] || { taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no \
        -R "@RG\tID:t\tSM:t" "$REF" "$JX/ctg.fa" 2> "$JX/mm2.log" \
        | samtools sort -@ 4 -m 3G -T "$JX/st" -o "$JX/ctg.bam" - && samtools index -@ 8 "$JX/ctg.bam"; }

    # min_supp is now a READ count, so sweep it on the scale reads actually live on
    for MS in 3 5 10 20; do
        ./svp call "$JX/ctg.bam" "$JX/ms$MS.vcf" 10 1000 50 50 $MS 6 300 0.01 2>/dev/null
        .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/ms$MS.vcf" \
            --min_size 50 --tolerance 500 --out_dir "$JX/minda" >/dev/null 2>&1
        awk -F'\t' -v m=$MS 'NR>1 && $1!~/^Minda/{printf "  min_reads=%-3s TP %5s FN %5s FP %6s  prec %.3f  rec %.3f  F1 %.3f\n",m,$2,$3,$4,$5,$6,$7}' \
            "$JX/minda/results/overall.tsv"
    done
done
echo "########## bitti $(date +%T) ##########"
