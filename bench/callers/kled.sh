#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/kled:1.2.11--he5eddf3_0
dk $IMG kled -R $REF --CCS -t $P $T > $O/t.vcf
dk $IMG kled -R $REF --CCS -t $P --F2 $N > $O/n.vcf
python3 $K/subtract.py $O/t.vcf $O/n.vcf --dist 1000 > $O/calls.vcf
done_msg
