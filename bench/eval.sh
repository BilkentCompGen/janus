#!/bin/bash
# eval.sh <pair> <vcf> <name>: PASS records only, scored into bench/eval/<pair>/<name>{.vcf,.log,/}. HG008: minda truthset mode
# against NIST. CASTLE pairs: (A) minda ensemble with the vcf as the fifth HiFi member of the CASTLE consensus rule (the
# fair rule) and (B) score.py against the CASTLE consensus without Severus, the vcf and Severus HiFi side by side.
p=$1; V=$2; tag=$3; E=/home/yuio/snape/bench/eval/$p; O=$E/$tag; mkdir -p $E
grep -E '^#|PASS' $V > $O.vcf
B=/home/yuio/snape/hifi/castle-calls/variant_calls_and_benchmarks; M=/home/yuio/.local/bin/minda
[ -e $O ] && { T=/home/yuio/trash/$(date +%Y%m%d)-bench-eval; mkdir -p $T; mv $O $T/$p.$tag.$(date +%H%M%S%N); }
if [ $p = hg008 ]; then
  $M truthset --base /home/yuio/snape/hifi/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz --vcfs $O.vcf --min_size 50 --tolerance 500 --out_dir $O > $O.log 2>&1
  awk -F'\t' 'NR==2{printf "%s\t%s\tA_nist\tF1 %.4f\tR %.4f\tP %.4f\n","'$p'","'$tag'",$7,$6,$5}' $O/results/overall.tsv
else
  P=${p^^}; D=$B/$P/$P
  $M ensemble --vcfs $D/severus_HIFI_vaf.vcf $D/savana_$P.haplotagged.classified.somatic_HIFI_vaf.vcf $D/tumor.nanomonsv.result_HIFI_vaf.vcf $D/sniffles2_HIFI_vaf.vcf $O.vcf \
    $D/severus_ont_vaf.vcf $D/savana_$P.haplotagged.classified.somatic_ONT_vaf.vcf $D/tumor.nanomonsv.result_ONT_vaf.vcf $D/sniffles2_ONT_vaf.vcf \
    $D/svaba_merged.vcf $D/${P}_subsample.bam.gripss.filtered_vaf.vcf $D/manta_somaticSV_vaf.vcf \
    --conditions "[[[caller_names[:5], caller_names[5:9], caller_names[9:]], '>=', 2],'&', [caller_names,'>=', 4]]" --min_size 50 --tolerance 500 --out_dir $O > $O.log 2>&1
  awk -F'\t' '$1~/^E_/{printf "%s\t%s\tA_fair\tF1 %.4f\tR %.4f\tP %.4f\n","'$p'","'$tag'",$7,$6,$5}' $O/results/overall.tsv
  for c in "snape $O.vcf" "severus $D/severus_HIFI_vaf.vcf"; do set -- $c; python3 /home/yuio/snape/bench/score.py /home/yuio/snape/hifi/truth/castle_without_severus/$p/None_minda_ensemble.vcf $2 | sed "s/^/$p\t$tag\tB_ws_$1\t/"; done
fi
