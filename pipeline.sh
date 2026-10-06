#!/bin/bash
D=$1; N=$(realpath $2); T=$(realpath $3); P=${4:-16}; B=/home/yuio/snape/build/snape; R=/home/yuio/snape/ref
mkdir -p $D; cd $D
step() { n=$1; shift; [ -e done.$n ] && return; /usr/bin/time -f "$n %e s %M KB %U user %S sys" "$@" > log.$n 2>&1 || exit 1; tail -1 log.$n >> steps.txt; touch done.$n; }
for s in n t; do f=$N; [ $s = t ] && f=$T; step pack_$s $B pack -i $f -o $s.svb -t $P ${BASES:+-b $BASES}; done
for s in n t; do step extract_$s $B extract -i $s.svb -o $s.svb.exc -t $P; done
X=$(( $(python3 -c "import struct; print(struct.unpack_from('<Q', open('t.svb', 'rb').read(32), 16)[0])") / 3100000000 ))
MN=3; MT=3; [ $X -lt 25 ] && MT=2; [ $X -lt 13 ] && MN=2
[ -n "$MINREADS" ] && MN=$MINREADS && MT=$MINREADS
step count_n $B count -i n.svb -o n.cnt -t $P -m $MN
step count_t $B count -i t.svb -o t.cnt -t $P -m $MT
step graph $B graph -i n.cnt,t.cnt -o graph.gph -m $MT
step prune $B prune -i graph.gph -r n.svb,t.svb -o prune -t $P -m $MT
step select $B select -i graph.gph -d prune.sites -r n.svb,t.svb -o ev -t $P
step place $B place -i ev.fa -r $R/${REF:-hg38_parmask}.cidx -d $R/${REF:-hg38_parmask}.fa -o ev -t $P
step local minimap2 -cx map-hifi --secondary=no -t $P -o local.paf ev.win.fa ev.on.fa
step lift $B lift -i local.paf -d $R/${REF:-hg38_parmask}.fa -o aln.paf
step clip $B clip -i ev -r aln.paf -o ev.gaps.fa
step gaps minimap2 -cx sr --secondary=no -t $P -o ev.gaps.raw.paf $R/${REF:-hg38_parmask}.sr.mmi ev.gaps.fa
step patch $B patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf
step call $B call -i ev -r n.svb,t.svb -d $R/${REF:-hg38_parmask}.fa -g graph.gph -o calls.vcf -v $R/$( [ -n "$REF" ] && echo $REF.trf.bed || echo vntr.trf.bed ) -t $P -m $MT
