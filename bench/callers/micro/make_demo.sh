#!/usr/bin/env bash
set -uo pipefail
M=/home/yuio/snape/bench/callers/micro/demo; S=/home/yuio/snape/bench/callers/ref/somasv/data/demo; IDX=/home/yuio/snape/ref/hg38_parmask.mmi
mkdir -p $M
for s in n t; do
  [ -s $M/$s.bam.bai ] && continue
  [ $s = n ] && { sm=normal; b=$S/HG008_10X_HiFi_chr20_normal_demo_subset.bam; } || { sm=tumor; b=$S/HG008_30X_HiFi_chr20_tumor_demo_subset.bam; }
  samtools fastq -F 0x900 $b 2>/dev/null | gzip -1 > $M/$s.fq.gz
  taskset -c 16-31 bash -c "/home/yuio/.local/bin/minimap2 -ax map-hifi -t 16 -Y --secondary=no -R '@RG\tID:$s\tSM:$sm' $IDX $M/$s.fq.gz 2> $M/$s.minimap2.log \
     | samtools sort -@ 4 -m 1G -T $M/st.$s -o $M/$s.bam - && samtools index $M/$s.bam"
done
