#!/usr/bin/env bash
# run.sh <experiment> [pair ...]: snape end to end, one run after the other with no pause, from hifi/<pair>/{n,t}.fastq.gz
# to calls.vcf (pipeline.sh, 16 threads, pinned to cores 0-15), each run measured by measure.py. Experiments:
#   full      one run per pair on all its reads                  -> runs/full/<pair>/
#   coverage  per pair, runs at 5x, 10x, ... up to 60x (as deep as both samples go; 1x = G bases): pack keeps only the
#             first N x G bases of each sample (whole reads, file order) -> runs/coverage/<pair>_<N>x/
#             runs go depth by depth (5x of every pair, then 10x, ...) so early results cover all pairs
# Per run folder: calls.vcf, run.json (every measurement), steps.tsv, files.tsv, logs/, monitor/ (raw samples and
# params.json). The intermediates are measured, then moved to ~/trash (full) or deleted (coverage, approved by the user:
# ~7 TB would not fit). A run with run.json is skipped; a full run waits while less than 400 GB is free.
S=/home/yuio/snape; H=$S/bench; X=${1:?experiment: full or coverage}; shift
G=3100000000                       # bases per 1x: the hg38 length, rounded
LEVELS="5 10 15 20 25 30 35 40 45 50 55 60"
ALL="hcc1954 h2009 hcc1937 h1437 hg008 hcc1395 colo829 hg009 h209 h2126"
PAIRS="${*:-$ALL}"
E=$H/runs/$X; Q=$E/queue.log; mkdir -p $E
log() { echo "$(date '+%F %T') $*" >> $Q; }

# one run: run <folder> <pair> <bases or empty>
run() {
    local D=$1 p=$2 b=$3
    [ -s $D/run.json ] && { log "[$(basename $D)] done before"; return; }
    if [ $X = full ]; then
        while [ $(df --output=avail -B1G $S | tail -1) -lt 400 ]; do log "[$p] waiting: under 400 GB free (empty ~/trash)"; sleep 600; done
    fi
    mkdir -p $D/work $D/monitor $D/logs; log "[$(basename $D)] start"
    printf '{"experiment": "%s", "pair": "%s", "coverage_x": %s, "bases_per_x": %s, "bases_limit": %s, "intermediates": "%s"}\n' \
        $X $p "${b:+$((b / G))}" $G "${b:-null}" $([ $X = coverage ] && echo deleted || echo trash) | sed 's/: ,/: null,/' > $D/monitor/params.json
    pidstat -h -r -u -d -C 'snape|minimap2' 2 >> $D/monitor/pidstat.txt 2>&1 & m1=$!
    iostat -dxm -t nvme0n1 5 >> $D/monitor/iostat.txt 2>&1 & m2=$!
    vmstat -t 5 >> $D/monitor/vmstat.txt 2>&1 & m3=$!
    ( while :; do echo "$(date +%s) $(cat /sys/class/hwmon/hwmon1/temp1_input) $(cat /sys/class/hwmon/hwmon1/temp3_input) $(sensors 2>/dev/null | awk '/Tctl/{print $2}')"; sleep 5; done ) >> $D/monitor/temp.txt & m4=$!
    [ -s $D/start.txt ] || [ -s $D/monitor/start.txt ] || date +%s.%N > $D/start.txt
    BASES=$b /usr/bin/time -f "total %e s wall %M KB peak %U s user %S s sys %I fs-in-blocks %O fs-out-blocks %x exit" \
        taskset -c 0-15 bash $S/pipeline.sh $D/work $S/hifi/$p/n.fastq.gz $S/hifi/$p/t.fastq.gz 16 2> $D/total.txt
    kill $m1 $m2 $m3 $m4 2>/dev/null; wait $m1 $m2 $m3 $m4 2>/dev/null
    [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ] || { log "[$(basename $D)] FAILED: $(tail -1 $D/total.txt)"; return; }
    python3 $H/measure.py $D >> $Q 2>&1 || { log "[$(basename $D)] measure.py failed, intermediates left in place"; return; }
    mv $D/work/calls.vcf $D/; mv $D/work/log.* $D/logs/
    if [ $X = coverage ]; then
        rm -rf $D/work; log "[$(basename $D)] done, intermediates deleted"
    else
        local T=/home/yuio/trash/$(date +%Y%m%d)-runs/$p; mkdir -p $(dirname $T); mv $D/work $T; log "[$p] done, intermediates in $T"
    fi
}

# depth of a pair, in whole x, from the full run's pack report (the shallower of its two samples)
depth() { python3 -c "import json; r=json.load(open('$H/runs/full/$1/run.json'))['reads']; print(int(min(r['n']['bases'], r['t']['bases']) // $G))"; }

if [ $X = full ]; then
    for p in $PAIRS; do run $E/$p $p ""; done
elif [ $X = coverage ]; then
    for n in $LEVELS; do
        for p in $PAIRS; do
            [ $n -le $(depth $p) ] && run $E/${p}_${n}x $p $((n * G))
        done
    done
else
    echo "run.sh: unknown experiment $X" >&2; exit 1
fi
log "queue done"
