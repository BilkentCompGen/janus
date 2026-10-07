#!/usr/bin/env bash
export GPU=1
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/svision-pro:2.5-gpu; R=$K/ref/svision-pro
NP=$(( P / 2 )); [ $NP -lt 1 ] && NP=1
mkdir -p $O/svp
dk $IMG SVision-pro --target_path $T --base_path $N --genome_path $REF --model_path $R/src/pre_process/model_liteunet_256_8_16_32_32_32.pth \
   --access_path $R/src/pre_process/hg38.access.10M.bed --out_path $O/svp --sample_name tumor --preset hifi --detect_mode somatic \
   --process_num $NP --unet_cpu_num 2 --device gpu --gpu_id 0
v=$(ls $O/svp/*.vcf | grep -v somatic_s | head -1)
dk $IMG python $R/extract_op.py --input_vcf $v --extract somatic
pass ${v%.vcf}.somatic_s1.vcf $O/calls.vcf
done_msg
