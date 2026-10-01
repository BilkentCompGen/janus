#!/usr/bin/env bash
# run.sh [pair ...]: snape end to end on each pair, one pair after the other with no pause, from hifi/<pair>/{n,t}.fastq.gz
# to calls.vcf (pipeline.sh unchanged, 16 threads, pinned to cores 0-15). Per pair, in bench/runs/<pair>/:
#   calls.vcf            the calls (all the verdict and the evaluation need)
#   run.json             every measurement in one file (see measure.py); steps.tsv and files.tsv are the same as tables
#   logs/                each step's own output
#   monitor/             the raw samples: pidstat (2 s, per snape/minimap2 process), iostat, vmstat, temperatures (5 s)
# The intermediates (svbs, counts, graph, alignments) are measured, then moved to ~/trash. Trash is on the same drive, so
# a pair waits while less than 400 GB is free (empty ~/trash to let it go on). Finished pairs (run.json present) are skipped.
S=/home/yuio/snape; H=$S/bench; Q=$H/runs/queue.log
PAIRS="${*:-h2009 hcc1937 h1437 hg008 hcc1395 colo829 hg009 h209 h2126}"
log() { echo "$(date '+%F %T') $*" >> $Q; }
for p in $PAIRS; do
    D=$H/runs/$p; [ -s $D/run.json ] && { log "[$p] done before"; continue; }
    while [ $(df --output=avail -B1G $S | tail -1) -lt 400 ]; do log "[$p] waiting: under 400 GB free (empty ~/trash)"; sleep 600; done
    mkdir -p $D/work $D/monitor $D/logs; log "[$p] start"
    pidstat -h -r -u -d -C 'snape|minimap2' 2 >> $D/monitor/pidstat.txt 2>&1 & m1=$!
    iostat -dxm -t nvme0n1 5 >> $D/monitor/iostat.txt 2>&1 & m2=$!
    vmstat -t 5 >> $D/monitor/vmstat.txt 2>&1 & m3=$!
    ( while :; do echo "$(date +%s) $(cat /sys/class/hwmon/hwmon1/temp1_input) $(cat /sys/class/hwmon/hwmon1/temp3_input) $(sensors 2>/dev/null | awk '/Tctl/{print $2}')"; sleep 5; done ) >> $D/monitor/temp.txt & m4=$!
    [ -s $D/start.txt ] || date +%s.%N > $D/start.txt
    /usr/bin/time -f "total %e s wall %M KB peak %U s user %S s sys %I fs-in-blocks %O fs-out-blocks %x exit" \
        taskset -c 0-15 bash $S/pipeline.sh $D/work $S/hifi/$p/n.fastq.gz $S/hifi/$p/t.fastq.gz 16 2> $D/total.txt
    kill $m1 $m2 $m3 $m4 2>/dev/null; wait $m1 $m2 $m3 $m4 2>/dev/null
    [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ] || { log "[$p] FAILED: $(tail -1 $D/total.txt)"; continue; }
    python3 $H/measure.py $D >> $Q 2>&1 || { log "[$p] measure.py failed, intermediates left in place"; continue; }
    mv $D/work/calls.vcf $D/; mv $D/work/log.* $D/logs/
    T=/home/yuio/trash/$(date +%Y%m%d)-runs/$p; mkdir -p $(dirname $T); mv $D/work $T
    log "[$p] done, intermediates in $T"
done
log "queue done"
