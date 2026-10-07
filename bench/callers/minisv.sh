#!/usr/bin/env bash
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
