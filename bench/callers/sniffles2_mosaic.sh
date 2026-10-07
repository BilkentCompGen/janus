#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
S=$L/sniffles/bin/sniffles
$S --input $T --vcf $O/t.mosaic.vcf --mosaic --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
$S --input $N --vcf $O/n.vcf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite --minsupport 1
python3 $K/subtract.py $O/t.mosaic.vcf $O/n.vcf --dist 1000 > $O/calls.vcf
done_msg
