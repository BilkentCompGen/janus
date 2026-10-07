#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
$L/severus/bin/severus --target-bam $T --control-bam $N --out-dir $O/severus --vntr-bed $VNTR -t $P
pass $O/severus/somatic_SVs/severus_somatic.vcf $O/calls.vcf
done_msg
