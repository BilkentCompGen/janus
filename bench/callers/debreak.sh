#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/debreak:1.3--h9ee0642_0
dk $IMG debreak --bam $T -o $O/t -t $P --min_size 50
dk $IMG debreak --bam $N -o $O/n -t $P --min_size 50 --min_support 2
python3 $K/subtract.py $O/t/debreak.vcf $O/n/debreak.vcf --dist 500 > $O/calls.vcf
done_msg
