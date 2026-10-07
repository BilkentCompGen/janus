#!/usr/bin/env bash
export GPU=1
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/svhunter:28d1ea5; H=/opt/SVHunter
C="[$(printf 'chr%s,' {1..22} X)chrY]"
W=$(( P / 4 )); [ $W -lt 1 ] && W=1
for s in t n; do
  b=$T; [ $s = n ] && b=$N
  mkdir -p $O/$s/data $O/$s/pred $O/$s/vcf
  dk $IMG python3 $H/SVHunter.py generate $b $O/$s/data $W "$C"
  dk $IMG python3 $H/SVHunter.py call $H/model_predict.h5 $O/$s/data $b $O/$s/pred $O/$s/vcf $P "$C" 1
  cp $O/$s/vcf/*/*_all.vcf $O/$s.vcf
done
python3 $K/subtract.py $O/t.vcf $O/n.vcf --dist 1000 | awk -F'\t' '/^#/ || !match($8,/SVLEN=[0-9]+/) || substr($8,RSTART+6,RLENGTH-6)+0>=50' > $O/calls.vcf
done_msg
