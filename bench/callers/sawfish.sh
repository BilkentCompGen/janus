#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/sawfish:2.2.1--h9ee0642_0
dk $IMG sawfish discover --threads $P --ref $REF --bam $T --disable-cnv --clobber --output-dir $O/t
dk $IMG sawfish discover --threads $P --ref $REF --bam $N --disable-cnv --clobber --output-dir $O/n
dk $IMG sawfish joint-call --threads $P --ref $REF --sample $O/t --sample $O/n --clobber --output-dir $O/joint
bcftools view -f PASS -i 'GT[0]="alt" && GT[1]="ref" && FMT/AD[1:1]=0 && (INFO/SVLEN="." || abs(INFO/SVLEN)>=50)' $O/joint/genotyped.sv.vcf.gz -Ov -o $O/calls.vcf
done_msg
