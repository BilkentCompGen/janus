#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/longcalld:0.0.11--hbc58adc_0
dk $IMG longcallD call -s -t $P $REF $T > $O/t.mosaic.vcf
dk $IMG longcallD call -t $P $REF $N > $O/n.vcf
python3 $K/subtract.py $O/t.mosaic.vcf $O/n.vcf --somatic-tag --dist 1000 > $O/calls.vcf
done_msg
