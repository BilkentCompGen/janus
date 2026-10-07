#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/cutesv:2.1.4--pyhdfd78af_1
HIFI="--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"
mkdir -p $O/t.work $O/n.work
dk $IMG cuteSV $T $REF $O/t.vcf $O/t.work $HIFI -l 50 -S tumor -t $P
dk $IMG cuteSV $N $REF $O/n.vcf $O/n.work $HIFI -l 50 -s 2 -S normal -t $P
python3 $K/subtract.py $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf
done_msg
