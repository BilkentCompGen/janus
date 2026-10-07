#!/usr/bin/env bash
set -u
S=/home/yuio/snape; D=$S/bench/runs/sim/source; F=$S/hifi/hg009/n.fastq.gz; LIM=$((64 * 3100000000)); mkdir -p $D
[ -s $D/n.bam.bai ] && exit 0
for f in $D/n.bam.tmp $D/st.*; do [ -e "$f" ] && mv "$f" /home/yuio/trash/2026-10-04-partial/sim.$(basename $f).$(date +%H%M%S); done
pigz -dc $F | awk -v lim=$LIM 'NR % 4 == 2 { s += length($0) } { print } NR % 4 == 0 && s >= lim { exit }' \
    | minimap2 -ax map-hifi -t 16 -Y --secondary=no -R '@RG\tID:simsrc\tSM:simsrc' $S/ref/hg38_parmask.mmi - 2> $D/minimap2.log \
    | samtools sort -@ 4 -m 2G -T $D/st -o $D/n.bam.tmp - && mv $D/n.bam.tmp $D/n.bam && samtools index -@ 8 $D/n.bam
echo "$(date '+%F %T') source alignment done: $(samtools view -c -@ 8 -F 0x900 $D/n.bam) reads" >> $S/bench/runs/night/queue.log
