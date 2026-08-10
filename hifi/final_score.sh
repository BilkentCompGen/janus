#!/usr/bin/env bash
# The last table: every sample called against its DEEP normal, at each support floor, beside the
# Severus output the paper published for the same sample and the same truth.
set -uo pipefail
cd /home/yuio/svp
CAS=/home/yuio/tmp/castle/variant_calls_and_benchmarks

for S in hcc1937 hcc1954 h2009 h1437; do
    JX=/home/yuio/tmp/jx/${S}_final
    [ "$S" = hcc1954 ] && JX=/home/yuio/tmp/jx/hcc1954_pad50
    N=/home/yuio/tmp/jx/$S/normal_deep.bam
    if [ ! -s "$N" ]; then echo "  $S: derin normal yok, atlaniyor"; continue; fi
    d=$(samtools idxstats "$N" | awk '{s+=$3} END{printf "%.1fM", s/1e6}')
    for SP in 3 4 5 6; do
        ./svp call "$JX/ctg.bam" "$N" "$JX/deep$SP.vcf" 10 0 50 50 $SP >/dev/null 2>&1
        .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$JX/deep$SP.vcf" \
            --min_size 50 --tolerance 500 --out_dir "$JX/mdeep$SP" >/dev/null 2>&1
        awk -F'\t' -v s=$S -v p=$SP -v d=$d 'NR>1 && $1!~/^Minda/{printf "  %-9s normal %-6s supp>=%s  TP %5s FN %5s FP %5s  prec %.3f  rec %.3f  F1 %.3f\n",s,d,p,$2,$3,$4,$5,$6,$7}' \
            "$JX/mdeep$SP/results/overall.tsv"
    done
    U=$(echo $S | tr 'a-z' 'A-Z')
    sev=$CAS/$U/$U/severus_HIFI_vaf.vcf
    if [ -s "$sev" ]; then
        .venv/bin/minda truthset --base hifi/$S/truth_clean.vcf --vcfs "$sev" \
            --min_size 50 --tolerance 500 --out_dir /home/yuio/tmp/sev_$S >/dev/null 2>&1
        awk -F'\t' -v s=$S 'NR>1 && $1!~/^Minda/{printf "  %-9s SEVERUS              TP %5s FN %5s FP %5s  prec %.3f  rec %.3f  F1 %.3f\n",s,$2,$3,$4,$5,$6,$7}' \
            /home/yuio/tmp/sev_$S/results/overall.tsv
    fi
    echo
done
echo "bitti $(date +%T)"
