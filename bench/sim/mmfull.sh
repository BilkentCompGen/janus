#!/usr/bin/env bash
S=/home/yuio/snape; B=$S/build/snape; R=$S/ref; P=16; G=3100000000; X=$S/bench/runs/mmfull; C="taskset -c 0-15"
L=$S/bench/runs/night/queue.log; TP=/home/yuio/trash/2026-10-04-partial
log() { echo "$(date '+%F %T') mmfull: $*" >> $L; }
step() { n=$1; shift; [ -e done.$n ] && return 0; /usr/bin/time -f "$n %e s %M KB" "$@" > log.$n 2>&1 || { log "$(pwd | sed "s#$X/##") $n FAILED"; return 1; }; tail -1 log.$n >> steps.txt; touch done.$n; }
variant() {
    local d=$1 src=$2 ncnt=$3 V=$4 cm=$5 f
    mkdir -p $d; cd $d
    [ -s calls.vcf ] && [ -e done.call ] && return 0
    for f in n.svb n.svb.idx n.svb.exc t.svb t.svb.idx t.svb.exc; do ln -sf $src/$f $f; done
    ln -sf $src/$ncnt n.cnt
    step count_t $C $B count -i t.svb -o t.cnt -t $P -m $V || return 1
    step graph $C $B graph -i n.cnt,t.cnt -o graph.gph -m $V || return 1
    step prune $C $B prune -i graph.gph -r n.svb,t.svb -o prune -t $P -m $V || return 1
    step select $C $B select -i graph.gph -d prune.sites -r n.svb,t.svb -o ev -t $P || return 1
    step place $C $B place -i ev.fa -r $R/hg38_parmask.cidx -d $R/hg38_parmask.fa -o ev -t $P || return 1
    step local $C minimap2 -cx map-hifi --secondary=no -t $P -o local.paf ev.win.fa ev.on.fa || return 1
    step lift $C $B lift -i local.paf -d $R/hg38_parmask.fa -o aln.paf || return 1
    step clip $C $B clip -i ev -r aln.paf -o ev.gaps.fa || return 1
    step gaps $C minimap2 -cx sr --secondary=no -t $P -o ev.gaps.raw.paf $R/hg38_parmask.sr.mmi ev.gaps.fa || return 1
    step patch $C $B patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf || return 1
    for m in $cm; do
        o=calls.vcf; [ "$m" != "${cm%% *}" ] && o=calls_m$m.vcf
        n=call; [ $o != calls.vcf ] && n=call$m
        step $n $C $B call -i ev -r n.svb,t.svb -d $R/hg38_parmask.fa -g graph.gph -o $o -v $R/vntr.trf.bed -t $P -m $m || return 1
    done
    rm -f t.cnt graph.gph ev.fa ev.win.fa ev.on.fa ev.off.fa local.paf aln.paf ev.paf ev.gaps.raw.paf
    log "$(echo $d | sed "s#$X/##"): $(awk '{s += $2; if ($4 > m) m = $4} END {printf "%.0f s, peak %.1f GB", s, m / 1048576}' steps.txt), $(grep -vc '^#' calls.vcf) records"
}
subset() {
    local W=$1 s
    mkdir -p $W/x10; cd $W/x10
    for s in n t; do
        [ -s $s.svb.idx ] || python3 $S/bench/svb_mix.py $W/x10/$s.svb $W/$s.svb:0:$((10 * G)) > log.mix_$s 2>&1 || return 1
        step extract_$s $C $B extract -i $W/x10/$s.svb -o $W/x10/$s.svb.exc -t $P || return 1
    done
    step count_n2 $C $B count -i n.svb -o n.cnt2 -t $P -m 2 || return 1
    step count_n3 $C $B count -i n.svb -o n.cnt3 -t $P -m 3 || return 1
}
mkdir -p $X $TP
log "start (benchmark queue stopped; every run alone on CPUs 0-15, 16 threads)"
for p in hcc1954 hcc1395 h2009 hcc1937 h1437 hg008 colo829; do
    W=$X/$p; mkdir -p $W; cd $W
    [ -e $W/done ] && continue
    ext=bam; [ -s $S/hifi/$p/t.fastq.gz ] && ext=fastq.gz
    for s in n t; do
        step pack_$s $C $B pack -i $S/hifi/$p/$s.$ext -o $W/$s.svb -t $P || continue 2
        step extract_$s $C $B extract -i $W/$s.svb -o $W/$s.svb.exc -t $P || continue 2
    done
    step count_n3 $C $B count -i n.svb -o n.cnt3 -t $P -m 3 || continue
    [ $ext = bam ] && { variant $W/full_m3 $W n.cnt3 3 "3" || continue; }
    variant $W/full_m2 $W n.cnt3 2 "2" || continue
    if [ $ext = bam ]; then
        subset $W || continue
        variant $W/x10_n2t2 $W/x10 n.cnt2 2 "2" || continue
        variant $W/x10_n3t2 $W/x10 n.cnt3 2 "2" || continue
        case $p in hcc1954|hcc1395) variant $W/x10_n2t1 $W/x10 n.cnt2 1 "2 1" ;; esac
        rm -f $W/x10/*.svb $W/x10/*.exc $W/x10/*.idx $W/x10/n.cnt2 $W/x10/n.cnt3
    fi
    rm -f $W/n.svb $W/t.svb $W/n.svb.exc $W/t.svb.exc $W/n.svb.idx $W/t.svb.idx $W/n.cnt3
    touch $W/done; log "$p finished"
done
cd $S/bench && OMP_NUM_THREADS=8 $S/.venv-ml/bin/python sim/mmfull_verdict.py > $X/verdict.txt 2>&1 && log "verdict done (bench/runs/mmfull/verdict.txt)" || log "verdict FAILED (see bench/runs/mmfull/verdict.txt)"
for t in cutesv kled dysgu svision-pro svdss delly; do
    [ -d $S/bench/runs/competitors/$t/hcc1954 ] && mv $S/bench/runs/competitors/$t/hcc1954 $TP/$t.hcc1954.cachepeak.$(date +%H%M%S)
done
for t in debreak sawfish longcalld; do [ -e $S/bench/runs/competitors/$t/hcc1954/attempts ] && [ ! -s $S/bench/runs/competitors/$t/hcc1954/run.json ] && mv $S/bench/runs/competitors/$t/hcc1954 $TP/$t.hcc1954.interrupted.$(date +%H%M%S); done
systemctl --user reset-failed snape-night 2> /dev/null
systemd-run --user --unit=snape-night -p Restart=on-failure -p RestartSec=60 -p OOMPolicy=continue --setenv=PATH="$PATH" --setenv=HOME="$HOME" bash $S/bench/night.sh
log "benchmark queue restarted"
