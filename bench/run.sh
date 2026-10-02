#!/usr/bin/env bash
# run.sh <experiment> [pair ...]: snape end to end, one run after the other with no pause, from hifi/<pair>/{n,t}.fastq.gz
# to calls.vcf (pipeline.sh, 16 threads, pinned to cores 0-15), each run measured by measure.py. Experiments:
#   full      one run per pair on all its reads                  -> runs/full/<pair>/
#   coverage  per pair, runs at 1x, 3x, 5x, 10x, 15x, ... up to 60x (as deep as both samples go; 1x = G bases): pack
#             keeps only the first N x G bases of each sample (whole reads, file order); up to 10x the support
#             threshold is 2 reads instead of 3 (see minreads below) -> runs/coverage/<pair>_<N>x/
#   groups    the 4 CASTLE pairs + HG008 on all reads with select's per-distance-group pieces (2026-10-02), measured like
#             full; intermediates deleted after recording (approved by the user)       -> runs/groups/<pair>/
#   threads   one pair (default hcc1954) on all its reads with 1, 2, 4, 8, 16 and 32 threads, each pinned to as many
#             cores (0..N-1): how well the pipeline parallelizes -> runs/threads/<pair>_<N>t/ (intermediates to trash)
#             runs go depth by depth (5x of every pair, then 10x, ...) so early results cover all pairs
# Per run folder: calls.vcf, run.json (every measurement), steps.tsv, files.tsv, logs/, monitor/ (raw samples and
# params.json). The intermediates are measured, then deleted (approved by the user, 2026-10-02). A run with run.json is
# skipped.
S=/home/yuio/snape; H=$S/bench; X=${1:?experiment: full or coverage}; shift
G=3100000000                       # bases per 1x: the hg38 length, rounded
LEVELS="1 3 5 10 15 20 25 30 35 40 45 50 55 60"
# per depth, chosen with the user (2026-10-01): support reads (snape -m: a core's occurrences to be a node, a tumor's
# reads on a novel edge, a site, a record) and the F-beta the verdict's cut is chosen for (precision first when shallow)
minreads() { [ $1 -le 10 ] && echo 2 || echo 3; }
beta() { [ $1 -le 5 ] && echo 0.5 || { [ $1 -le 10 ] && echo 0.75 || echo 1; }; }
ALL="hcc1954 h2009 hcc1937 h1437 hg008 hcc1395 colo829 hg009 h209 h2126"
PAIRS="${*:-$ALL}"
E=$H/runs/$X; Q=$E/queue.log; mkdir -p $E
log() { echo "$(date '+%F %T') $*" >> $Q; }

# one run: run <folder> <pair> <bases or empty> <depth or empty> [threads]
run() {
    local D=$1 p=$2 b=$3 n=$4 th=${5:-16} m="" f=1
    [ -n "$n" ] && { m=$(minreads $n); f=$(beta $n); }
    [ -s $D/run.json ] && { log "[$(basename $D)] done before"; return; }
    mkdir -p $D/work $D/monitor $D/logs; log "[$(basename $D)] start"
    printf '{"experiment": "%s", "pair": "%s", "coverage_x": %s, "bases_per_x": %s, "bases_limit": %s, "min_reads": %s, "verdict_f_beta": %s, "threads": %s, "intermediates": "%s"}\n' \
        $X $p "${n:-null}" $G "${b:-null}" "${m:-3}" $f $th deleted > $D/monitor/params.json
    pidstat -h -r -u -d -C 'snape|minimap2' 2 >> $D/monitor/pidstat.txt 2>&1 & m1=$!
    iostat -dxm -t nvme0n1 5 >> $D/monitor/iostat.txt 2>&1 & m2=$!
    vmstat -t 5 >> $D/monitor/vmstat.txt 2>&1 & m3=$!
    ( while :; do echo "$(date +%s) $(cat /sys/class/hwmon/hwmon1/temp1_input) $(cat /sys/class/hwmon/hwmon1/temp3_input) $(sensors 2>/dev/null | awk '/Tctl/{print $2}')"; sleep 5; done ) >> $D/monitor/temp.txt & m4=$!
    [ -s $D/start.txt ] || [ -s $D/monitor/start.txt ] || date +%s.%N > $D/start.txt
    BASES=$b MINREADS=$m /usr/bin/time -f "total %e s wall %M KB peak %U s user %S s sys %I fs-in-blocks %O fs-out-blocks %x exit" \
        taskset -c 0-$((th - 1)) bash $S/pipeline.sh $D/work $S/hifi/$p/n.fastq.gz $S/hifi/$p/t.fastq.gz $th 2> $D/total.txt
    kill $m1 $m2 $m3 $m4 2>/dev/null; wait $m1 $m2 $m3 $m4 2>/dev/null
    [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ] || { log "[$(basename $D)] FAILED: $(tail -1 $D/total.txt)"; return; }
    python3 $H/measure.py $D >> $Q 2>&1 || { log "[$(basename $D)] measure.py failed, intermediates left in place"; return; }
    mv $D/work/calls.vcf $D/; mv $D/work/log.* $D/logs/
    rm -rf $D/work; log "[$(basename $D)] done, intermediates deleted"
}

# depth of a pair, in whole x, from the full run's pack report (the shallower of its two samples)
depth() { python3 -c "import json; r=json.load(open('$H/runs/full/$1/run.json'))['reads']; print(int(min(r['n']['bases'], r['t']['bases']) // $G))"; }

if [ $X = full ]; then
    for p in $PAIRS; do run $E/$p $p "" ""; done
elif [ $X = groups ]; then
    for p in ${*:-hcc1954 h2009 hcc1937 h1437 hg008}; do run $E/$p $p "" ""; done
elif [ $X = threads ]; then
    for p in ${*:-hcc1954}; do for th in 1 2 4 8 16 32; do run $E/${p}_${th}t $p "" "" $th; done; done
elif [ $X = coverage ]; then
    for n in $LEVELS; do
        for p in $PAIRS; do
            [ $n -le $(depth $p) ] && run $E/${p}_${n}x $p $((n * G)) $n
        done
    done
else
    echo "run.sh: unknown experiment $X" >&2; exit 1
fi
log "queue done"
