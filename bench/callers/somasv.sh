#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/somasv:0.0.1-92f1560
dk $IMG python /opt/SomaSV/main.py --tumor-bam $T --normal-bam $N --reference $REF --output-dir $O/somasv \
   --sample-platform HIFI --mode long-read-only --processes $P 2>&1 | grep -v ' - DEBUG - '
v=$O/somasv/final_somatic_variants.vcf
if [ ! -e $v ] && grep -q '^Final somatic variants: 0$' $O/log.txt; then printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > $v; fi
pass $v $O/calls.vcf
done_msg
