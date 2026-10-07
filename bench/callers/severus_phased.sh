#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=hkubal/clair3:v2.0.3
dk $IMG /opt/bin/run_clair3.sh --bam_fn=$N --ref_fn=$REF --threads=$P --platform=hifi --model_path=/opt/models/hifi_revio \
   --output=$O/clair3 --enable_phasing --longphase_for_phasing
V=$O/clair3/phased_merge_output.vcf.gz
dk $IMG /opt/bin/longphase haplotag -s $V -b $N -r $REF -t $P --tagSupplementary -o $O/n.hp
dk $IMG /opt/bin/longphase haplotag -s $V -b $T -r $REF -t $P --tagSupplementary -o $O/t.hp
samtools index -@ $P $O/n.hp.bam; samtools index -@ $P $O/t.hp.bam
$L/severus/bin/severus --target-bam $O/t.hp.bam --control-bam $O/n.hp.bam --phasing-vcf $V --out-dir $O/severus --vntr-bed $VNTR -t $P
pass $O/severus/somatic_SVs/severus_somatic.vcf $O/calls.vcf
done_msg
