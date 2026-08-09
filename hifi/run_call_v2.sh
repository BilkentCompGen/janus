#!/usr/bin/env bash
# End-to-end re-call of the real pairs with the two measured fixes in `svp call`:
#   - a substitution-mismatch gate (a contig placed where it does not belong is not evidence)
#   - a germline check against the normal READS at the reference locus, which is the question the
#     graph cannot ask: it decides somatic by exact core sequence, and in a repeat the same germline
#     variant does not spell the same core in both samples.
#
# graph / bubbles / junctions are UNCHANGED code, so their outputs are reused rather than recomputed;
# the only new input is a normal alignment, and 2M reads (~11x) is what the measurement used.
set -uo pipefail
cd /home/yuio/svp
REF=/home/datasets/references/genomes/hg38/hg38.fa

for S in hcc1937 hcc1954 h2009 h1437; do
    JX=/home/yuio/tmp/jx/$S
    N=$JX/normal2M.bam
    [ "$S" = hcc1954 ] && N=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam   # already built
    echo "########## $S  $(date +%T) ##########"

    if [ ! -s "$N" ]; then
        echo "  [normal hizalama] $(date +%T)"
        ./svp dump hifi/$S/n.svb | head -8000000 \
          | taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no "$REF" - 2> "$JX/mm2_normal.log" \
          | samtools sort -@ 4 -m 3G -T "$JX/stn" -o "$N" -
        samtools index -@ 8 "$N"
    else
        echo "  [skip] normal hizalama var"
    fi

    ./svp call "$JX/ctg.bam" "$N" "$JX/calls_v2.vcf" 2>&1 | /usr/bin/grep -E "mismatch>|elenen|supp>="
    .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/calls_v2.vcf" \
        --min_size 50 --tolerance 500 --out_dir "$JX/minda_v2" >/dev/null 2>&1
    awk -F'\t' -v s=$S 'NR>1 && $1!~/^Minda/{printf "  >> %-9s TP %5s  FN %5s  FP %6s  prec %.3f  rec %.3f  F1 %.3f\n",s,$2,$3,$4,$5,$6,$7}' \
        "$JX/minda_v2/results/overall.tsv"
done
echo "########## bitti $(date +%T) ##########"
