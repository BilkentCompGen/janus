#!/usr/bin/env bash
H=/home/yuio/snape/bench; L=$H/runs/all.log; Q=$H/runs/competitors/queue.log; W=$H/runs/competitors/work
log() { echo "$(date '+%F %T') $*" >> $L; }
until grep -q "\[h2009\] n: .* reads aligned" $Q; do sleep 20; done
for pid in $(pgrep -f "bash align_queue.sh"); do kill -- -$(ps -o pgid= $pid | tr -d ' ') 2>/dev/null; done
sleep 5
T=/home/yuio/trash/2026-10-03-align-partial; mkdir -p $T
for f in $W/h2009/t.bam.tmp $W/h2009/t.bam.st.*; do [ -e "$f" ] && mv "$f" $T/; done
[ -e $H/runs/align/h2009/t ] && [ ! -s $H/runs/align/h2009/t/run.json ] && mv $H/runs/align/h2009/t $T/align_h2009_t
log "alignments: remaining ones with 32 threads"; bash $H/align_fast.sh
log "all done"
