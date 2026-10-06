#!/usr/bin/env bash
S=/home/yuio/snape; H=$S/bench; X=${1:?experiment: full or coverage}; shift
G=3100000000
DEV=$(basename $(findmnt -no SOURCE /home) | sed 's/p[0-9]*$//'); HW=$(ls -d /sys/block/$DEV/device/hwmon* | head -1)
LEVELS="5 10 15 20 30"
CASTLE="hcc1954 h2009 hcc1937 h1437 hcc1395"
minreads() { [ $1 -lt 13 ] && echo '"n2 t2"' || { [ $1 -lt 25 ] && echo '"n3 t2"' || echo '"n3 t3"'; }; }
beta() { [ $1 -le 5 ] && echo 0.5 || { [ $1 -le 10 ] && echo 0.75 || echo 1; }; }
ALL="hcc1954 h2009 hcc1937 h1437 hg008 hcc1395 colo829 hg009 h209 h2126"
PAIRS="${*:-$ALL}"
E=$H/runs/$X; Q=$E/queue.log; mkdir -p $E
src() { [ -s $S/hifi/$1/$2.bam ] && echo $S/hifi/$1/$2.bam || echo $S/hifi/$1/$2.fastq.gz; }
log() { echo "$(date '+%F %T') $*" >> $Q; }

run() {
    local D=$1 p=$2 b=$3 n=$4 th=${5:-16} m="" f=1
    [ -n "$n" ] && { m=$(minreads $n); f=$(beta $n); }
    [ -s $D/run.json ] && { log "[$(basename $D)] done before"; return; }
    mkdir -p $D/work $D/monitor $D/logs; log "[$(basename $D)] start"
    for c in $(src $p n) $(src $p t) $D/work/n.svb $D/work/t.svb; do [ -e $c ] && dd if=$c iflag=nocache count=0 status=none; done
    sync; sleep 30
    printf '{"experiment": "%s", "pair": "%s", "coverage_x": %s, "bases_per_x": %s, "bases_limit": %s, "min_reads": %s, "verdict_f_beta": %s, "threads": %s, "intermediates": "%s"}\n' \
        $X $p "${n:-null}" $G "${b:-null}" "${m:-\"n3 t3\"}" $f $th deleted > $D/monitor/params.json
    pidstat -h -r -u -d -C 'snape|minimap2' 2 >> $D/monitor/pidstat.txt 2>&1 & m1=$!
    iostat -dxm -t $DEV 5 >> $D/monitor/iostat.txt 2>&1 & m2=$!
    vmstat -t 5 >> $D/monitor/vmstat.txt 2>&1 & m3=$!
    ( while :; do echo "$(date +%s) $(cat $HW/temp1_input) $(cat $HW/temp3_input) $(sensors 2>/dev/null | awk '/Tctl/{print $2}')"; sleep 5; done ) >> $D/monitor/temp.txt & m4=$!
    [ -s $D/start.txt ] || [ -s $D/monitor/start.txt ] || date +%s.%N > $D/start.txt
    BASES=$b /usr/bin/time -f "total %e s wall %M KB peak %U s user %S s sys %I fs-in-blocks %O fs-out-blocks %x exit" \
        taskset -c 0-$((th - 1)) bash $S/pipeline.sh $D/work $(src $p n) $(src $p t) $th 2> $D/total.txt
    kill $m1 $m2 $m3 $m4 2>/dev/null; wait $m1 $m2 $m3 $m4 2>/dev/null
    [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ] || { log "[$(basename $D)] FAILED: $(tail -1 $D/total.txt)"; return; }
    python3 $H/measure.py $D >> $Q 2>&1 || { log "[$(basename $D)] measure.py failed, intermediates left in place"; return; }
    mv $D/work/calls.vcf $D/; mv $D/work/log.* $D/logs/
    rm -rf $D/work; log "[$(basename $D)] done, intermediates deleted"
}

depth() { python3 -c "import json; r=json.load(open('$H/runs/full/$1/run.json'))['reads']; print(int(min(r['n']['bases'], r['t']['bases']) // $G))"; }

if [ $X = full ]; then
    for p in $PAIRS; do run $E/$p $p "" ""; done
elif [ $X = groups ]; then
    for p in ${*:-hcc1954 h2009 hcc1937 h1437 hg008}; do run $E/$p $p "" ""; done
elif [ $X = threads ]; then
    for p in ${*:-hcc1954}; do for th in 1 2 4 8 12 16; do run $E/${p}_${th}t $p "" "" $th; done; done
elif [ $X = coverage ]; then
    for p in ${*:-$CASTLE}; do
        M=$E/mix/$p; mkdir -p $M
        for s in n t; do [ -s $M/$s.svb.idx ] || taskset -c 0-15 $S/build/snape pack -i $(src $p $s) -o $M/$s.svb -t 16 > $M/pack_$s.log 2>&1; done
        for n in $LEVELS; do
            [ $n -le $(depth $p) ] || continue
            D=$E/${p}_${n}x; [ -s $D/run.json ] && continue
            mkdir -p $D/work
            for s in n t; do python3 $H/svb_mix.py $D/work/$s.svb $M/$s.svb:0:$((n * G)) >> $Q 2>&1 || { log "[${p}_${n}x] subset FAILED"; continue 2; }; done
            touch $D/work/done.pack_n $D/work/done.pack_t
            run $D $p "" $n
        done
        rm -rf $M; log "[$p] full packed reads deleted"
    done
elif [ $X = purity ]; then
    for p in ${*:-hcc1954 hcc1395}; do
        M=$E/mix/$p; mkdir -p $M
        for s in n t; do [ -s $M/$s.svb.idx ] || taskset -c 0-15 $S/build/snape pack -i $(src $p $s) -o $M/$s.svb -t 16 > $M/pack_$s.log 2>&1; done
        nb=$(python3 -c "import struct; print(struct.unpack_from('<Q', open('$M/n.svb', 'rb').read(32), 16)[0])")
        nm=$((nb - 24 * G - G / 10)); [ $nm -gt $((30 * G)) ] && nm=$((30 * G))
        [ -s $M/nm.svb.idx ] || python3 $H/svb_mix.py $M/nm.svb $M/n.svb:0:$nm >> $Q 2>&1 || { log "[$p] matched normal mix FAILED"; continue; }
        for q in 20 40 60 80 100; do
            [ -s $M/t$q.svb.idx ] || python3 $H/svb_mix.py $M/t$q.svb $M/t.svb:0:$((q * 30 * G / 100)) $([ $q -lt 100 ] && echo $M/n.svb:$nm:$(((100 - q) * 30 * G / 100))) >> $Q 2>&1 || { log "[${p}_${q}p] mix FAILED, run skipped"; continue; }
            D=$E/${p}_${q}p; mkdir -p $D/work
            for f in svb svb.idx; do ln -sf $M/nm.$f $D/work/n.$f; ln -sf $M/t$q.$f $D/work/t.$f; done
            touch $D/work/done.pack_n $D/work/done.pack_t
            run $D $p "" ""
        done
        rm -rf $M; log "[$p] purity mixes deleted"
    done
else
    echo "run.sh: unknown experiment $X" >&2; exit 1
fi
log "queue done"
