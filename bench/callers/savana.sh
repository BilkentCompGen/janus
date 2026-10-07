#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
$L/savana/bin/savana --tumour $T --normal $N --ref $REF --outdir $O/savana --threads $P --pb --sample s --overwrite
pass $O/savana/s.classified.somatic.vcf $O/calls.vcf
done_msg
