#!/usr/bin/env bash
# Re-run everything downstream of the graph with the two changes:
#   - `svp bubbles` drops a somatic arm whose sibling arm is normal-only (one variant written twice)
#   - `svp call` has no segment length floor
# graph.bin is unchanged and reused; drop.sigs, ctg.fa, ctg.bam and the VCF are rebuilt into a
# separate <S>_gfix directory so the committed run stays on disk to compare against.
set -uo pipefail
cd /home/yuio/svp
SVP=/home/yuio/svp
REF=/home/datasets/references/genomes/hg38/hg38.fa

for S in "$@"; do
    G=$SVP/outputs/svp/$S
    JX=/home/yuio/tmp/jx/${S}_gfix
    OLD=/home/yuio/tmp/jx/$S
    N=$OLD/normal2M.bam
    [ "$S" = hcc1954 ] && N=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam
    mkdir -p "$JX"
    echo "########## $S  $(date +%T) ##########"

    if [ ! -s "$JX/drop.sigs" ]; then
        echo "  [bubbles] $(date +%T)"
        taskset -c 0-15 $SVP/svp bubbles "$G/graph.bin" <($SVP/svp dump hifi/$S/n.svb) \
            <($SVP/svp dump hifi/$S/t.svb) 16 21 22 40 12 20 2>&1 | /usr/bin/grep -E "bubble|kardes|cekirdek"
        cp "$G/drop.sigs" "$JX/drop.sigs"
    else
        echo "  [skip] bubbles"
        cp "$JX/drop.sigs" "$G/drop.sigs"
    fi

    if [ ! -s "$JX/ctg.fa" ]; then
        echo "  [junctions] $(date +%T)"
        taskset -c 0-15 $SVP/svp junctions "$G/graph.bin" <($SVP/svp dump hifi/$S/t.svb) \
            "$JX/ctg.fa" 16 25 5 300 0 2>&1 | tail -3
    fi

    if [ ! -s "$JX/ctg.bam" ]; then
        echo "  [align] $(date +%T)"
        taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no "$REF" "$JX/ctg.fa" 2> "$JX/mm2.log" \
            | samtools sort -@ 4 -m 3G -T "$JX/st" -o "$JX/ctg.bam" -
        samtools index -@ 8 "$JX/ctg.bam"
    fi

    $SVP/svp call "$JX/ctg.bam" "$N" "$JX/calls.vcf" 2>&1 | /usr/bin/grep -E "elenen|aday"
    .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/calls.vcf" \
        --min_size 50 --tolerance 500 --out_dir "$JX/minda" >/dev/null 2>&1
    awk -F'\t' -v s=$S 'NR>1 && $1!~/^Minda/{printf "  >> %-9s TP %5s  FN %5s  FP %6s  prec %.3f  rec %.3f  F1 %.3f\n",s,$2,$3,$4,$5,$6,$7}' \
        "$JX/minda/results/overall.tsv"
done
echo "########## bitti $(date +%T) ##########"
