#!/usr/bin/env bash
S=/home/yuio/snape; M=$S/bench/runs/sim; B=$M/source/n.bam; G=3100000000; L=$S/bench/runs/night/queue.log
log() { echo "$(date '+%F %T') sim: $*" >> $L; }
until [ -s $B.bai ]; do sleep 120; done
for spec in "mixed none 1" "bnd BND 6" "del DEL 2" "ins INS 3" "dup DUP 4" "inv INV 5"; do
    set -- $spec; D=$M/$1; mkdir -p $D
    [ -e $D/done ] && continue
    if [ ! -e $D/implant.done ]; then
        log "$1: implanting (dominant type $2)"
        nice python3 $S/bench/sim/implant.py --bam $B --out $D/sim --events 2500 --dominant $2 --seed $3 > $D/implant.log 2>&1 && touch $D/implant.done || { log "$1: implant FAILED (see $D/implant.log)"; continue; }
        log "$1: $(tail -2 $D/implant.log | tr '\n' ' ')"
    fi
    for x in 30 20 10; do
        W=$D/x$x; [ -e $W/done.call ] && [ -s $W/calls.vcf ] && continue
        mkdir -p $W; m=3; [ $x = 10 ] && m=2
        BASES=$((x * G)) MINREADS=$m nice bash $S/pipeline.sh $W $D/sim.n.bam $D/sim.t.bam 16 > $W/run.out 2>&1 || { log "$1 ${x}x: snape FAILED"; continue; }
        for f in $W/*.svb $W/*.exc $W/*.idx $W/*.cnt $W/graph.gph $W/ev.*fa $W/*.paf; do [ -e $f ] && rm -f $f; done
        log "$1 ${x}x: $(grep -vc '^#' $W/calls.vcf) records"
    done
    [ -s $D/x30/calls.vcf ] && [ -s $D/x20/calls.vcf ] && [ -s $D/x10/calls.vcf ] && { rm -f $D/sim.n.bam $D/sim.t.bam; touch $D/done; }
done
log "all sets done"
