#!/usr/bin/env bash
p=${1:-hcc1395}; P=H1395; C=/home/yuio/snape/hifi/castle-calls/variant_calls_and_benchmarks/HCC1395/HCC1395; T=/home/yuio/snape/hifi/truth; M=/home/yuio/.local/bin/minda
SAV_H=$C/savana_$P.haplotagged.classified.somatic_HIFI_vaf.vcf; SAV_O=$C/savana_$P.haplotagged.classified.somatic_ONT_vaf.vcf
SR="$C/svaba_merged.vcf $C/${P}_realigned.bam.gripss.filtered_vaf.vcf $C/manta_somaticSV_vaf.vcf"
mkdir -p $T/castle_without_snape $T/castle_without_severus
$M ensemble --vcfs $C/severus_HIFI_vaf.vcf $SAV_H $C/tumor.nanomonsv.result_HIFI_vaf.vcf $C/sniffles2_HIFI_vaf.vcf \
    $C/severus_ont_vaf.vcf $SAV_O $C/tumor.nanomonsv.result_ONT_vaf.vcf $C/sniffles2_ONT_vaf.vcf $SR \
    --conditions "[[[caller_names[:4], caller_names[4:8], caller_names[8:]], '>=', 2],'&', [caller_names,'>=', 4]]" \
    --min_size 50 --tolerance 500 --out_dir $T/castle_without_snape/$p > $T/castle_without_snape/$p.log 2>&1
$M ensemble --vcfs $SAV_H $C/tumor.nanomonsv.result_HIFI_vaf.vcf $C/sniffles2_HIFI_vaf.vcf \
    $SAV_O $C/tumor.nanomonsv.result_ONT_vaf.vcf $C/sniffles2_ONT_vaf.vcf $SR \
    --conditions "[[[caller_names[:3], caller_names[3:6], caller_names[6:9]], '>=', 2],'&', [caller_names,'>=', 4]]" \
    --min_size 50 --tolerance 500 --out_dir $T/castle_without_severus/$p > $T/castle_without_severus/$p.log 2>&1
for d in castle_without_snape castle_without_severus; do echo "$d/$p: $(grep -vc '^#' $T/$d/$p/None_minda_ensemble.vcf) events"; done
