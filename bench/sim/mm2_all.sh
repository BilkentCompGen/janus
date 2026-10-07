#!/usr/bin/env bash
S=/home/yuio/snape; B=$S/build/snape; R=$S/ref; P=16; G=3100000000; L=$S/bench/runs/night/queue.log
log() { echo "$(date '+%F %T') sim: mm2 $*" >> $L; }
for p in h2009 hcc1937 h1437; do
    W=$S/bench/runs/mm2/${p}_20x; mkdir -p $W; cd $W
    step() { n=$1; shift; [ -e done.$n ] && return 0; nice /usr/bin/time -f "$n %e s %M KB" "$@" > log.$n 2>&1 || { log "$p $n FAILED"; return 1; }; tail -1 log.$n >> steps.txt; touch done.$n; }
    [ -s m2/calls.vcf ] && [ -s m3/calls.vcf ] && continue
    for s in n t; do
        step pack_$s $B pack -i $S/hifi/$p/$s.bam -o $W/$s.svb -t $P -b $((20 * G)) || continue 2
        step extract_$s $B extract -i $W/$s.svb -o $W/$s.svb.exc -t $P || continue 2
    done
    step count_n $B count -i $W/n.svb -o $W/n.cnt -t $P -m 3 || continue
    for V in 3 2; do
        mkdir -p $W/m$V; cd $W/m$V
        for f in n.svb n.svb.idx n.svb.exc t.svb t.svb.idx t.svb.exc n.cnt; do ln -sf $W/$f $f; done
        step count_t $B count -i t.svb -o t.cnt -t $P -m $V || continue 2
        step graph $B graph -i n.cnt,t.cnt -o graph.gph -m $V || continue 2
        step prune $B prune -i graph.gph -r n.svb,t.svb -o prune -t $P -m $V || continue 2
        step select $B select -i graph.gph -d prune.sites -r n.svb,t.svb -o ev -t $P || continue 2
        step place $B place -i ev.fa -r $R/hg38_parmask.cidx -d $R/hg38_parmask.fa -o ev -t $P || continue 2
        step local minimap2 -cx map-hifi --secondary=no -t $P -o local.paf ev.win.fa ev.on.fa || continue 2
        step lift $B lift -i local.paf -d $R/hg38_parmask.fa -o aln.paf || continue 2
        step clip $B clip -i ev -r aln.paf -o ev.gaps.fa || continue 2
        step gaps minimap2 -cx sr --secondary=no -t $P -o ev.gaps.raw.paf $R/hg38_parmask.sr.mmi ev.gaps.fa || continue 2
        step patch $B patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf || continue 2
        step call $B call -i ev -r n.svb,t.svb -d $R/hg38_parmask.fa -g graph.gph -o calls.vcf -v $R/vntr.trf.bed -t $P -m $V || continue 2
        rm -f t.cnt graph.gph ev.fa ev.win.fa ev.on.fa ev.off.fa *.paf
        cd $W
    done
    rm -f $W/n.svb $W/t.svb $W/n.svb.exc $W/t.svb.exc $W/n.cnt
    log "${p} 20x done: $(grep -vc '^#' m3/calls.vcf) records at MM=3, $(grep -vc '^#' m2/calls.vcf) at MM=2"
done
cd $S/bench && OMP_NUM_THREADS=8 nice $S/.venv-ml/bin/python sim/mm2_verdict.py > runs/mm2/verdict.txt 2>&1 && log "verdict done (bench/runs/mm2/verdict.txt)" || log "verdict FAILED"
