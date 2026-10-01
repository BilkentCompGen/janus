#!/usr/bin/env bash
# somasv.sh <tumor.bam> <normal.bam> <outdir> <threads>: SomaSV 0.0.1 (eioyuou/SomaSV @ 92f1560, image
# snape-bench/somasv:0.0.1-92f1560, images/somasv/Dockerfile), long-read-only mode, HIFI platform, no PoN
# (the PoN is a hybrid-mode input) -> calls.vcf = PASS records of final_somatic_variants.vcf.
# Its per-function DEBUG trace is dropped from log.txt.
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/somasv:0.0.1-92f1560
dk $IMG python /opt/SomaSV/main.py --tumor-bam $T --normal-bam $N --reference $REF --output-dir $O/somasv \
   --sample-platform HIFI --mode long-read-only --processes $P 2>&1 | grep -v ' - DEBUG - '
v=$O/somasv/final_somatic_variants.vcf
# SomaSV writes no VCF at all when it has 0 calls: then calls.vcf is a header-only VCF
if [ ! -e $v ] && grep -q '^Final somatic variants: 0$' $O/log.txt; then printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > $v; fi
pass $v $O/calls.vcf
done_msg
