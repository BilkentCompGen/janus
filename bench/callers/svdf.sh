#!/usr/bin/env bash
export GPU=1
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/svdf:41c1607; D=/opt/SVDF/svdf.py
dk $IMG python3 $D call $T --working_dir $O/t --ref $REF --read_type CCS -s 3 -t $P --min_sv_size 50
dk $IMG python3 $D call $N --working_dir $O/n --ref $REF --read_type CCS -s 2 -t $P --min_sv_size 50 --mode sensitive
python3 $K/subtract.py $O/t/variants.vcf $O/n/variants.vcf --dist 1000 > $O/calls.vcf
done_msg
