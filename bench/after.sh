#!/usr/bin/env bash
H=/home/yuio/snape/bench; L=$H/runs/all.log; V=/home/yuio/snape/.venv-ml/bin/python
log() { echo "$(date '+%F %T') $*" >> $L; }
until grep -q "all done" $L; do sleep 300; done
log "after: coverage 40x"; bash $H/run.sh coverage
log "after: coverage verdict"; ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 $V cov_verdict_matched.py > verdict/coverage_stdout.txt 2>&1 )
T=/home/yuio/trash/2026-10-03-purity-hcc1395-first; mkdir -p $T
for q in 20 40 60 80 100; do [ -e $H/runs/purity/hcc1395_${q}p ] && mv $H/runs/purity/hcc1395_${q}p $T/; done
[ -e $H/runs/purity/mix/hcc1395 ] && mv $H/runs/purity/mix/hcc1395 $T/mix
log "after: purity hcc1395"; bash $H/run.sh purity hcc1395
log "after: purity verdict"; ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 $V purity_verdict.py > verdict/purity_stdout.txt 2>&1 )
log "after: full hs578t"; bash $H/run.sh full hs578t
log "after: done"
