#!/usr/bin/env bash
V=${1:-2}; S=/home/yuio/snape; B=$S/build/snape; R=$S/ref; P=8; C="nice -n 15 taskset -c 24-31"
for p in hcc1954 hcc1395; do
    D=$S/bench/runs/trace/${p}_20x; W=$S/bench/runs/trace/${p}_20x_mm$V; mkdir -p $W; cd $W
    for f in n.svb n.svb.idx n.svb.exc t.svb t.svb.idx t.svb.exc n.cnt; do ln -sf $D/$f $W/$f; done
    step() { n=$1; shift; [ -e done.$n ] && return 0; /usr/bin/time -f "$n %e s %M KB" "$@" > log.$n 2>&1 || { echo "$p $n FAILED" >> $W/exp.log; return 1; }; tail -1 log.$n >> steps.txt; touch done.$n; }
    step count_t $C $B count -i t.svb -o t.cnt -t $P -m $V || continue
    step graph $C $B graph -i n.cnt,t.cnt -o graph.gph -m $V || continue
    step prune $C $B prune -i graph.gph -r n.svb,t.svb -o prune -t $P -m $V || continue
    step select $C $B select -i graph.gph -d prune.sites -r n.svb,t.svb -o ev -t $P || continue
    step place $C $B place -i ev.fa -r $R/hg38_parmask.cidx -d $R/hg38_parmask.fa -o ev -t $P || continue
    step local $C minimap2 -cx map-hifi --secondary=no -t $P -o local.paf ev.win.fa ev.on.fa || continue
    step lift $C $B lift -i local.paf -d $R/hg38_parmask.fa -o aln.paf || continue
    step clip $C $B clip -i ev -r aln.paf -o ev.gaps.fa || continue
    step gaps $C minimap2 -cx sr --secondary=no -t $P -o ev.gaps.raw.paf $R/hg38_parmask.sr.mmi ev.gaps.fa || continue
    step patch $C $B patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf || continue
    step call3 $C $B call -i ev -r n.svb,t.svb -d $R/hg38_parmask.fa -g graph.gph -o calls_m3.vcf -v $R/vntr.trf.bed -t $P -m 3
    step callv $C $B call -i ev -r n.svb,t.svb -d $R/hg38_parmask.fa -g graph.gph -o calls_mv.vcf -v $R/vntr.trf.bed -t $P -m $V
    echo "$p done $(date '+%T')" >> $W/exp.log
done
