#!/usr/bin/env bash
# make.sh: the micro validation pair: the first 131,000 lines (~32.7k reads, ~0.5 Gb) of hifi/hcc1954/{n,t}.fastq.gz,
# aligned as hifi/align.sh does (minimap2 -ax map-hifi -Y --secondary=no, ref/hg38_parmask.mmi) with a read group
# (SM normal / tumor), sorted and indexed -> micro/{n,t}.bam(.bai). Pinned to cores 16-31.
set -uo pipefail
M=/home/yuio/snape/bench/callers/micro; H=/home/yuio/snape/hifi/hcc1954; IDX=/home/yuio/snape/ref/hg38_parmask.mmi
for s in n t; do
  [ -s $M/$s.bam.bai ] && continue
  sm=$([ $s = n ] && echo normal || echo tumor)
  zcat $H/$s.fastq.gz | head -n 131000 | gzip -1 > $M/$s.fq.gz
  /usr/bin/time -f "align_$s %e s %M KB" taskset -c 16-31 bash -c "/home/yuio/.local/bin/minimap2 -ax map-hifi -t 16 -Y --secondary=no -R '@RG\tID:$s\tSM:$sm' $IDX $M/$s.fq.gz 2> $M/$s.minimap2.log \
     | samtools sort -@ 4 -m 1G -T $M/st.$s -o $M/$s.bam - && samtools index $M/$s.bam" 2>> $M/align.times.txt
done
echo done > $M/make.done
