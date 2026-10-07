#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/nanovar:1.8.3--py311h26ae33e_1
ln -sf $REF $O/ref.fa; cp $REF.fai $O/ref.fa.fai
dk $IMG nanovar -t $P -x pacbio-ccs -f hg38 $T $O/ref.fa $O/t
dk $IMG nanovar -t $P -x pacbio-ccs -f hg38 -c 1 -s 0 $N $O/ref.fa $O/n
python3 $K/subtract.py $O/t/*.nanovar.pass.vcf $O/n/*.nanovar.pass.vcf --dist 1000 > $O/calls.vcf
done_msg
