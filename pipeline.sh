#!/bin/bash
# pipeline.sh <dir> <normal.fastq.gz> <tumor.fastq.gz> [threads]: pack (the reads into svbs, both samples side by side),
# extract, count, graph, prune, select, place, local alignment, lift,
# the pieces' unaligned stretches aligned on their own (clip, short-read mode, patch), and call into <dir>. Each finished step leaves done.<step>, so a rerun continues where it stopped.
D=$1; N=$(realpath $2); T=$(realpath $3); P=${4:-16}; B=/home/yuio/snape/build/snape; R=/home/yuio/snape/ref
mkdir -p $D; cd $D
step() { n=$1; shift; [ -e done.$n ] && return; /usr/bin/time -f "$n %e s %M KB %U user %S sys" "$@" > log.$n 2>&1 || exit 1; tail -1 log.$n >> steps.txt; touch done.$n; }
step pack_n $B pack -i $N -o n.svb -t $((P / 2)) & step pack_t $B pack -i $T -o t.svb -t $((P / 2)) & wait
[ -e done.pack_n ] && [ -e done.pack_t ] || exit 1
for s in n t; do step extract_$s $B extract -i $s.svb -o $s.svb.exc -t $P; done
for s in n t; do step count_$s $B count -i $s.svb -o $s.cnt -t $P; done
step graph $B graph -i n.cnt,t.cnt -o graph.gph
step prune $B prune -i graph.gph -r n.svb,t.svb -o prune -t $P
step select $B select -i graph.gph -d prune.sites -r n.svb,t.svb -o ev -t $P
step place $B place -i ev.fa -r $R/hg38_parmask.cidx -d $R/hg38_parmask.fa -o ev -t $P
step local minimap2 -cx map-hifi --secondary=no -t $P -o local.paf ev.win.fa ev.on.fa
step lift $B lift -i local.paf -d $R/hg38_parmask.fa -o aln.paf
step clip $B clip -i ev -r aln.paf -o ev.gaps.fa
step gaps minimap2 -cx sr --secondary=no -t $P -o ev.gaps.raw.paf $R/hg38_parmask.sr.mmi ev.gaps.fa
step patch $B patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf
step call $B call -i ev -r n.svb,t.svb -d $R/hg38_parmask.fa -g graph.gph -o calls.vcf -v $R/vntr.trf.bed -t $P
