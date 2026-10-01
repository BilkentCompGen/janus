#!/usr/bin/env bash
# minisv.sh <tumor.bam> <normal.bam> <outdir> <threads>: minisv r29 (image snape-bench/minisv:r29, images/minisv/Dockerfile),
# the README's tumor-normal recipe WITHOUT the normal's own assembly (no hifiasm: > 40 GB): the reads come back out of
# the BAMs (primary records); tumor aligned to GRCh38 (ref/hg38_parmask) + T2T-CHM13v2.0 (UCSC hs1; minimap2 -cx map-hifi -s50
# --ds) + the CHM13 minigraph pangenome (CHM13-464, minigraph -cxlr); normal aligned to GRCh38; minisv e / extract with
# the hs38 centromere mask, merge, keep TUMOR-only candidates (grep TUMOR | grep -v NORMAL), genvcf -> calls.vcf
# (GRCh38 coordinates). The tumor FASTQ (t.fq.gz) stays in the outdir.
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/minisv:r29; R=$K/ref; CEN=/opt/minisv/data/hs38.cen-mask.bed; H38=/home/yuio/snape/ref/hg38_parmask.mmi
[ -s $R/hs1.map-hifi.mmi ] || dk $IMG minimap2 -x map-hifi -t $P -d $R/hs1.map-hifi.mmi $R/hs1.fa.gz
dk $IMG bash -c "set -euo pipefail
samtools fastq -F 0x900 -@ 4 $T | gzip -1 > $O/t.fq.gz
minimap2 -cx map-hifi -s50 --ds -t $P $H38 $O/t.fq.gz | gzip -1 > $O/t.hg38l.paf.gz
minimap2 -cx map-hifi -s50 --ds -t $P $R/hs1.map-hifi.mmi $O/t.fq.gz | gzip -1 > $O/t.chm13l.paf.gz
minigraph -cxlr -t $P $R/CHM13-464.gfa.gz $O/t.fq.gz | gzip -1 > $O/t.chm13g.paf.gz
samtools fastq -F 0x900 -@ 4 $N | minimap2 -cx map-hifi -s50 --ds -t $P $H38 - | gzip -1 > $O/n.hg38l.paf.gz
minisv.js e -n TUMOR -b $CEN $O/t.hg38l.paf.gz $O/t.chm13l.paf.gz $O/t.chm13g.paf.gz | bash > $O/t.rsv
minisv.js extract -n NORMAL -b $CEN $O/n.hg38l.paf.gz > $O/n.rsv
cat $O/t.rsv $O/n.rsv | sort -k1,1 -k2,2n -S4G -T $O | minisv.js merge - > $O/merged.msv
(grep TUMOR $O/merged.msv | grep -v NORMAL || true) > $O/paired.msv
minisv.js genvcf $O/paired.msv > $O/calls.vcf"
done_msg
