#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/dysgu:1.8.8--py313h5d164f8_0
dk $IMG dysgu call --mode pacbio-revio --diploid False -p $P -x -o $O/t.vcf $REF $O/tmp_t $T
dk $IMG dysgu call --mode pacbio-revio -p $P -x -o $O/n.vcf $REF $O/tmp_n $N
dk $IMG dysgu filter -p $P --normal-vcf $O/n.vcf -o $O/somatic.vcf $O/t.vcf $N
bcftools view -f PASS -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="TRA" || abs(INFO/SVLEN)>=50' $O/somatic.vcf -Ov -o $O/calls.vcf
done_msg
