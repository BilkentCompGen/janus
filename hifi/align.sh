#!/usr/bin/env bash
# align.sh <pair> [threads]: hifi/<pair>/{n,t}.fastq.gz -> hifi/<pair>/{n,t}.bam (+ .bai), the BAMs the other callers
# need. minimap2 -ax map-hifi -Y --secondary=no against ref/hg38_parmask.mmi, samtools sort + index; normal and tumor run
# side by side with <threads> each (default 16). The BAM must hold every read of its fastq.gz (primary records), or it is
# left as <s>.bam.bad. Per-sample wall time, CPU time and peak memory go to hifi/<pair>/align.steps.txt. Run it in tmux.
set -uo pipefail
p=$1; T=${2:-16}; H=/home/yuio/snape/hifi; D=$H/$p; IDX=/home/yuio/snape/ref/hg38_parmask.mmi
[ -s $D/n.fastq.gz ] && [ -s $D/t.fastq.gz ] || { echo "$p: fastq.gz missing"; exit 1; }
one() {   # $1 n|t
    local s=$1 want got
    [ -s $D/$s.bam.bai ] && return 0
    /usr/bin/time -f "align_$s %e s %M KB %U user %S sys" bash -c "minimap2 -ax map-hifi -t $T -Y --secondary=no $IDX $D/$s.fastq.gz 2> $D/$s.minimap2.log \
        | samtools sort -@ 4 -m 2G -T $D/st.$s -o $D/$s.bam.tmp - && mv $D/$s.bam.tmp $D/$s.bam && samtools index -@ 8 $D/$s.bam" 2>> $D/align.steps.txt
    want=$(( $(pigz -dc $D/$s.fastq.gz | wc -l) / 4 )); got=$(samtools view -c -@ 8 -F 0x900 $D/$s.bam)
    [ "$want" = "$got" ] && echo "$p $s: $got reads aligned" || { mv $D/$s.bam $D/$s.bam.bad; echo "[FATAL] $p $s: $got reads in the BAM, $want in the fastq.gz"; }
}
one n & one t & wait
