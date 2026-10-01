#!/usr/bin/env bash
# longcalld.sh <tumor.bam> <normal.bam> <outdir> <threads>: longcallD 0.0.11 (quay.io/biocontainers/longcalld:0.0.11)
# mosaic mode (-s) on the tumor, germline mode on the normal (both HiFi defaults); calls.vcf = tumor SV records
# (>= 50 bp; longcallD calls only INS/DEL) tagged SOMATIC, PASS, with no normal SV of the same type within 1 kb
# (subtract.py). longcallD has no tumor/normal mode; the subtraction is ours.
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/longcalld:0.0.11--hbc58adc_0
dk $IMG longcallD call -s -t $P $REF $T > $O/t.mosaic.vcf
dk $IMG longcallD call -t $P $REF $N > $O/n.vcf
python3 $K/subtract.py $O/t.mosaic.vcf $O/n.vcf --somatic-tag --dist 1000 > $O/calls.vcf
done_msg
