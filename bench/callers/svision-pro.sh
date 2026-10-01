#!/usr/bin/env bash
# svision-pro.sh <tumor.bam> <normal.bam> <outdir> <threads>: SVision-pro 2.5 (image snape-bench/svision-pro:2.5-pysam0.21 =
# quay.io/biocontainers/svision-pro:2.5 + pysam 0.21, images/svision-pro/Dockerfile; the stock pysam 0.9.1 crashes),
# somatic mode, HiFi preset, the 256 LiteUnet model and the hg38 access bed from the v2.5 repo (ref/svision-pro), CPU;
# then the repo's extract_op.py --extract somatic (records tagged _NewCOMP, no POLY) -> calls.vcf
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/svision-pro:2.5-pysam0.21; R=$K/ref/svision-pro
NP=$(( P / 2 )); [ $NP -lt 1 ] && NP=1      # NP prediction processes x 2 cores each (--unet_cpu_num 2)
mkdir -p $O/svp
dk $IMG SVision-pro --target_path $T --base_path $N --genome_path $REF --model_path $R/src/pre_process/model_liteunet_256_8_16_32_32_32.pth \
   --access_path $R/src/pre_process/hg38.access.10M.bed --out_path $O/svp --sample_name tumor --preset hifi --detect_mode somatic \
   --process_num $NP --unet_cpu_num 2 --device cpu
v=$(ls $O/svp/*.vcf | grep -v somatic_s | head -1)
dk $IMG python $R/extract_op.py --input_vcf $v --extract somatic
pass ${v%.vcf}.somatic_s1.vcf $O/calls.vcf
done_msg
