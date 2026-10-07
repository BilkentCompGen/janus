#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/svdss:2.1.1--he17396a_0
FMD=$K/ref/svdss/hg38_parmask.fmd
[ -s $K/ref/svdss/index.done ] || { echo "missing $FMD: run $K/ref/svdss/index.sh"; exit 1; }
run() {
  dk $IMG SVDSS smooth --reference $REF --bam $1 --threads $P > $O/$2.smoothed.bam
  samtools index -@ $P $O/$2.smoothed.bam
  dk $IMG SVDSS search --index $FMD --bam $O/$2.smoothed.bam --threads $P > $O/$2.sfs.txt
  dk $IMG SVDSS call --reference $REF --bam $O/$2.smoothed.bam --sfs $O/$2.sfs.txt --threads $P --min-sv-length 50 --min-cluster-weight $3 > $O/$2.vcf
}
run $T t 4
run $N n 2
python3 $K/subtract.py $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf
done_msg
