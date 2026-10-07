#!/usr/bin/env bash
S=/home/yuio/snape; H=$S/bench; R=$H/runs; N=$R/night; K=$H/callers; G=3100000000
V=$S/.venv-ml/bin/python; Q=$N/queue.log; TR=/home/yuio/trash/2026-10-04-fastq; TP=/home/yuio/trash/2026-10-04-partial
CASTLE="hcc1954 h2009 hcc1937 h1437 hcc1395"
mkdir -p $N/running $R/competitors $R/align $TP
echo 1000 > /proc/self/oom_score_adj
log() { echo "$(date '+%F %T') $*" >> $Q; }
note() { echo "$(date '+%F %T') $*" >> $N/INCIDENTS.md; log "INCIDENT: $*"; }
free_gb() { df --output=avail -BG /home | tail -1 | tr -dc 0-9; }
wait_disk() { local w=0; while [ $(free_gb) -lt 60 ]; do [ $w = 0 ] && note "disk: only $(free_gb) GB free; queue waits (fastq files in ~/trash still occupy space)"; w=1; sleep 600; done; }
bam_ok() { [ -s $S/hifi/$1/$2.bam.bai ] && [ -s $R/align/$1/$2/run.json ] && grep -q "\[$1\] $2: [0-9]* reads aligned" $R/competitors/queue.log; }

align_lane() {
    local tries=0 miss p s f
    while :; do
        miss=0; for p in $CASTLE; do for s in n t; do bam_ok $p $s || miss=1; done; done
        [ $miss = 0 ] && break
        if pgrep -f "align_fast.sh" > /dev/null; then sleep 60; continue; fi
        tries=$((tries + 1)); [ $tries -gt 3 ] && { note "alignments: BAMs still missing after 3 restarts; gave up"; return; }
        note "alignments: the alignment job was not running while BAMs are missing; (re)started, attempt $tries"
        for p in $CASTLE; do
            for f in $S/hifi/$p/*.bam.tmp $S/hifi/$p/*.bam.st.*; do [ -e "$f" ] && mv "$f" $TP/$p.$(basename $f).$tries; done
            for s in n t; do [ -d $R/align/$p/$s ] && [ ! -s $R/align/$p/$s/run.json ] && mv $R/align/$p/$s $TP/align_${p}_${s}.$tries; done
        done
        bash $H/align_fast.sh >> $N/align.out 2>&1
    done
    log "alignments: all ten BAMs present and read counts verified"
}

trash_fq() {
    local p=$1
    [ -e $S/hifi/$p/n.fastq.gz ] || [ -e $S/hifi/$p/t.fastq.gz ] || return
    bam_ok $p n && bam_ok $p t || return
    mkdir -p $TR/$p; mv $S/hifi/$p/n.fastq.gz $S/hifi/$p/t.fastq.gz $TR/$p/ && log "fastq.gz of $p moved to $TR/$p"
}

trace() {
    local p=$1 D=$R/trace/${1}_20x
    [ -e $D/done.call ] && return
    [ -s $S/hifi/$p/t.fastq.gz ] || { note "trace $p 20x: fastq.gz missing, not run"; return; }
    mkdir -p $D; log "trace $p 20x start (intermediates kept)"
    BASES=$((20 * G)) MINREADS=3 taskset -c 0-15 bash $S/pipeline.sh $D $S/hifi/$p/n.fastq.gz $S/hifi/$p/t.fastq.gz 16 > $D/run.out 2>&1 \
        && log "trace $p 20x done" || note "trace $p 20x FAILED (see $D/run.out)"
}

purity() {
    local q
    if [ ! -e $N/purity.moved ]; then
        mkdir -p $TP/purity-hcc1395-first
        for q in 20 40 60 80 100; do [ -e $R/purity/hcc1395_${q}p ] && mv $R/purity/hcc1395_${q}p $TP/purity-hcc1395-first/; done
        [ -e $R/purity/mix/hcc1395 ] && mv $R/purity/mix/hcc1395 $TP/purity-hcc1395-first/mix
        touch $N/purity.moved
    fi
    if [ ! -s $R/purity/hcc1395_100p/run.json ] || [ ! -s $R/purity/hcc1395_20p/run.json ]; then
        [ -s $S/hifi/hcc1395/t.fastq.gz ] || { note "purity hcc1395: fastq.gz missing, not run"; return; }
        log "purity hcc1395 rerun start"; bash $H/run.sh purity hcc1395
        ( cd $H && OMP_NUM_THREADS=8 $V purity_verdict.py > verdict/purity_stdout.txt 2>&1 ) && log "purity verdict done" || note "purity verdict FAILED"
        [ -s $R/purity/hcc1395_20p/run.json ] || note "purity hcc1395 20% still has no result after the rerun"
    fi
}

hs578t() {
    [ -s $R/full/hs578t/run.json ] && return
    [ -s $S/hifi/hs578t/t.fastq.gz ] || { note "hs578t: fastq.gz missing"; return; }
    log "hs578t full run start"; bash $H/run.sh full hs578t
    [ -s $R/full/hs578t/calls.vcf ] && ( cd $H && OMP_NUM_THREADS=8 $V verdict.py > verdict/stdout.txt 2>&1 ) && log "verdict rerun with hs578t: $(grep hs578t $H/verdict/log.txt)" || note "hs578t full run or verdict FAILED"
}

snape_lane() {
    trash_fq h2009
    trace hcc1954; [ -e $R/trace/hcc1954_20x/done.call ] && trash_fq hcc1954
    purity
    trace hcc1395
    hs578t
    log "snape follow-ups finished"
}

one() {
    local t=$1 p=$2 cpus=$3 mem=$4 kind=$5 D=$R/competitors/$1/$2 W=$S/hifi/$2 a sl c pid f
    [ -s $D/run.json ] && return 0
    [ -f $K/$t.sh ] || return 0
    [ -s $W/t.bam.bai ] && [ -s $W/n.bam.bai ] || return 0
    a=$(cat $D/attempts 2> /dev/null || echo 0); [ $a -ge 2 ] && return 0
    wait_disk
    [ $a -ge 1 ] && note "$t/$p: an earlier attempt left no result (killed, timed out or the queue was restarted); second and last attempt"
    mkdir -p $TP; [ -e $D/out ] && mv $D/out $TP/$t.$p.out.$(date +%H%M%S)
    mkdir -p $D/out; echo $((a + 1)) > $D/attempts; chmod +x $K/$t.sh
    echo "{\"experiment\": \"competitor\", \"tool\": \"$t\", \"pair\": \"$p\", \"coverage_x\": null, \"cpus\": \"$cpus\", \"threads\": 16, \"memory_cap\": \"$mem\", \"run\": \"$kind\"}" > $D/params.json
    sl="sb_$(echo $t | tr -c 'a-z0-9\n' '_')_$p.slice"
    log "[$t/$p] start ($kind, cpus $cpus)"
    CPUS=$cpus DMEM=$mem CGPARENT=$sl setsid /usr/bin/time -f "%e %M %U %S %x" -o $D/time.txt timeout ${TMO:-172800} taskset -c $cpus $K/$t.sh $W/t.bam $W/n.bam $D/out 16 > $D/out/stdout.txt 2>&1 &
    pid=$!; echo $pid > $N/running/$t.$p.pgid; c=/sys/fs/cgroup/$sl
    ( m=0; while kill -0 $pid 2> /dev/null; do x=$(awk '/^anon /{a=$2} /^shmem /{s=$2} END{print a+s}' $c/memory.stat 2> /dev/null); [ "${x:-0}" -gt $m ] && { m=$x; echo $m > $D/anon_peak.txt; }; sleep 2; done ) &
    wait $pid; : > $N/running/$t.$p.pgid
    [ -d $c ] && { grep -E "^(usage|user|system)_usec" $c/cpu.stat; echo "memory_peak $(cat $c/memory.peak)"; [ -s $D/anon_peak.txt ] && echo "anon_peak $(cat $D/anon_peak.txt)"; } > $D/cgroup.txt
    [ -d $c ] && rmdir $c 2> /dev/null
    if python3 $H/measure_tool.py $D >> $Q 2>&1; then
        for f in calls.vcf log.txt stdout.txt; do [ -e $D/out/$f ] && mv $D/out/$f $D/; done
        rm -rf $D/out
        python3 -c "import json,sys; j=json.load(open('$D/run.json')); sys.exit(0 if j['total']['exit'] == 0 and j['calls'] is not None else 1)" \
            || note "$t/$p: the tool failed (exit $(python3 -c "import json; print(json.load(open('$D/run.json'))['total']['exit'])"), see $D/log.txt)"
    else
        note "$t/$p: no measurement (the job was killed before it finished)"
    fi
}

tools() { grep -v '^#' $N/$1.txt 2> /dev/null; }
lane() { local t p; for t in $(tools $2); do for p in $1; do one $t $p 0-31 60g shared; done; done; }
others() { lane "h2009 hcc1937" $1 & lane "h1437 hcc1395" $1 & wait; }
clean() { local t; for t in $(tools $1); do one $t hcc1954 0-15 100g clean; done; }
gpu() { local t p; for t in $(tools gpu); do for p in $CASTLE; do one $t $p 0-15 150g alone-gpu; done; done; }
state() { ls $R/competitors/*/*/run.json $R/competitors/*/*/attempts 2> /dev/null | wc -l; }

log "queue (re)started"
align_lane & A=$!
snape_lane
wait $A
for p in $CASTLE; do trash_fq $p; done
log "phase: severus (unphased), one at a time"
for p in $CASTLE; do one severus $p 0-15 100g clean; done
while :; do
    b=$(state)
    if [ ! -e $N/hold_clean ]; then
        log "phase: fast tools, hcc1954 alone"; clean fast
        log "phase: slow tools, hcc1954 alone"; clean slow
    fi
    log "phase: fast tools, other pairs, two at a time"; others fast
    log "phase: slow tools, other pairs, two at a time"; others slow
    [ -e $N/hold_clean ] || { log "phase: GPU tools, one at a time"; gpu; }
    [ $(state) = $b ] && { log "nothing left to run; sleeping"; sleep 900; }
done
