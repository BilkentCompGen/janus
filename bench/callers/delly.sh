#!/usr/bin/env bash
# delly.sh <tumor.bam> <normal.bam> <outdir> <threads>: Delly 2.6.0 long-read mode (quay.io/biocontainers/delly:2.6.0),
# joint discovery on tumor + normal (-y pb, delly's hg38 exclude list), then delly filter -f somatic with the
# tumor/control sample sheet -> calls.vcf = PASS records
source /home/yuio/snape/bench/callers/common.sh
IMG=quay.io/biocontainers/delly:2.6.0--h3752d28_0
sm() { samtools view -H $1 | awk -F'\t' '/^@RG/{for(i=2;i<=NF;i++) if($i~/^SM:/){print substr($i,4); exit}}'; }
ST=$(sm $T); SN=$(sm $N); [ -n "$ST" ] && [ -n "$SN" ] && [ "$ST" != "$SN" ] || { echo "tumor/normal BAMs need distinct @RG SM tags"; exit 1; }
cp $K/ref/human.hg38.excl.tsv $O/excl.tsv
printf "%s\ttumor\n%s\tcontrol\n" "$ST" "$SN" > $O/samples.tsv
dk $IMG delly lr -y pb -g $REF -x $O/excl.tsv -o $O/joint.bcf $T $N
dk $IMG delly filter -f somatic -o $O/somatic.bcf -s $O/samples.tsv $O/joint.bcf
bcftools view -f PASS $O/somatic.bcf -Ov -o $O/calls.vcf
done_msg
