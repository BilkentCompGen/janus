#!/usr/bin/env bash
# One table: every sample, every support floor, svp against the CASTLE Severus output, on both the
# 4-support ensemble the paper reports against and the 2-support one. Severus is one of the callers
# that built both, so its calls are over-represented in each; the pair of numbers is the honest read,
# not either alone.
#
# usage: scoreboard.sh <jx-suffix>     e.g. `scoreboard.sh final` reads /home/yuio/tmp/jx/<S>_final
set -uo pipefail
cd /home/yuio/svp
SUF=${1:-final}
CAS=/home/yuio/tmp/castle/variant_calls_and_benchmarks

row() {  # row <label> <vcf> <sample>
    for t in 4 2; do
        base=hifi/$3/truth_clean.vcf
        [ $t = 2 ] && base=/home/yuio/tmp/ens/${3}_2_1.vcf
        [ -s "$base" ] || continue
        o=/home/yuio/tmp/sb/$3_$1_$t
        .venv/bin/minda truthset --base "$base" --vcfs "$2" --min_size 50 --tolerance 500 \
            --out_dir "$o" >/dev/null 2>&1
        awk -F'\t' -v l="$1" -v s="$3" -v t="$t" 'NR>1 && $1!~/^Minda/{printf "  %-9s %-12s %s-destek  TP %5s FN %5s FP %5s  prec %.3f  rec %.3f  F1 %.3f\n",s,l,t,$2,$3,$4,$5,$6,$7}' \
            "$o/results/overall.tsv" 2>/dev/null
    done
}

mkdir -p /home/yuio/tmp/sb /home/yuio/tmp/ens
for S in hcc1937 hcc1954 h2009 h1437; do
    U=$(echo $S | tr 'a-z' 'A-Z')
    src=$CAS/$U/$U/minda_2_1/None_minda_ensemble.vcf
    [ -s "$src" ] && [ ! -s /home/yuio/tmp/ens/${S}_2_1.vcf ] && \
        /usr/bin/grep -v '^cmd:' "$src" > /home/yuio/tmp/ens/${S}_2_1.vcf
    JX=/home/yuio/tmp/jx/${S}_${SUF}
    for SP in 3 4 5; do
        [ -s "$JX/c$SP.vcf" ] && row "supp>=$SP" "$JX/c$SP.vcf" "$S"
    done
    sev=$CAS/$U/$U/severus_HIFI_vaf.vcf
    [ -s "$sev" ] && row "severus" "$sev" "$S"
    echo
done
