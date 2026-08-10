#!/usr/bin/env bash
# The whole downstream at the configuration every sweep converged on, on whatever samples are named.
#
#   pad_cores 50   more anchor per contig. Cannot cost a contig -- the window is clamped to the
#                  read's own cores and min_flank tests the flank achieved -- and measured on
#                  hcc1954 it dominates pad 25 at every support level: same precision, +0.009 recall.
#   min_aln 0      a length floor on a SEGMENT breaks split chains and fabricates junctions between
#                  pieces that were never adjacent.
#   shred_rep 0.40 a dense run of small indels is an inversion in unique sequence and a VNTR in
#                  repetitive sequence; measured on the contig's own bases, which need no reference.
#   anchor_rep 0.50 the same question at the call: 27.2% of the calls nothing corroborates sit above
#                  0.4 against 2.9% of the true ones. 0.40 is past the knee -- it costs 33 true calls.
#   normal_dirty 0.25  the "the normal cannot place this locus either" clause, and it now runs on
#                  breakends too, which the same-chromosome guard used to skip.
#
# graph.bin and drop.sigs are reused unchanged; only junctions onward is rebuilt.
set -uo pipefail
cd /home/yuio/svp
SVP=${SVP:-/home/yuio/tmp/svp_new}
REF=/home/datasets/references/genomes/hg38/hg38.fa
PAD=${PAD:-50}

for S in "$@"; do
    G=/home/yuio/svp/outputs/svp/$S
    JX=/home/yuio/tmp/jx/${S}_final
    N=/home/yuio/tmp/jx/$S/normal2M.bam
    [ "$S" = hcc1954 ] && N=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam
    mkdir -p "$JX"
    echo "########## $S  $(date +%T) ##########"

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
    for SP in 3 4 5 6; do
        $SVP call "$JX/ctg.bam" "$N" "$JX/c$SP.vcf" 10 0 50 50 $SP 6 300 0.01 0.40 0.50 0.25 >/dev/null 2>&1
        .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/c$SP.vcf" \
            --min_size 50 --tolerance 500 --out_dir "$JX/m$SP" >/dev/null 2>&1
        awk -F'\t' -v s=$S -v p=$SP 'NR>1 && $1!~/^Minda/{printf "  >> %-9s supp>=%s  TP %5s FN %5s FP %5s  prec %.3f  rec %.3f  F1 %.3f\n",s,p,$2,$3,$4,$5,$6,$7}' \
            "$JX/m$SP/results/overall.tsv"
    done
done
echo "########## bitti $(date +%T) ##########"
