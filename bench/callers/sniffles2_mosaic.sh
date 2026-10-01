#!/usr/bin/env bash
# sniffles2_mosaic.sh <tumor.bam> <normal.bam> <outdir> <threads>: Sniffles2 2.8.0 (venv lib/callers/sniffles) --mosaic
# on the tumor (single-sample low-VAF mode), germline calling on the normal with --minsupport 1; calls.vcf = the PASS
# mosaic tumor SVs with no normal SV of the same type within 1 kb (subtract.py). Sniffles2 has no tumor/normal mode;
# the subtraction is ours (the Severus paper reports that --mosaic lowered Sniffles2's F1 vs. the multisample run).
source /home/yuio/snape/bench/callers/common.sh
S=$L/sniffles/bin/sniffles
$S --input $T --vcf $O/t.mosaic.vcf --mosaic --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
$S --input $N --vcf $O/n.vcf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite --minsupport 1
python3 $K/subtract.py $O/t.mosaic.vcf $O/n.vcf --dist 1000 > $O/calls.vcf
done_msg
