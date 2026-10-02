#!/usr/bin/env bash
# all.sh: every benchmark of the final snape (tag snape-final), one after the other on cores 0-15, progress in
# runs/all.log: the full runs of the ten pairs, the verdict, minda's scoring (CASTLE consensus rule / NIST), the coverage
# series, the thread scaling (HCC1954, 1-32 threads), then the competitors. Each stage resumes where it stopped.
H=/home/yuio/snape/bench; L=$H/runs/all.log; mkdir -p $H/runs
log() { echo "$(date '+%F %T') $*" >> $L; }
log "full runs"; bash $H/run.sh full
log "verdict"; ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 /home/yuio/snape/.venv-ml/bin/python verdict.py > verdict/stdout.txt 2>&1 )
log "minda"; for p in hcc1954 h2009 hcc1937 h1437 hg008; do taskset -c 0-15 bash $H/eval.sh $p $H/verdict/$p.vcf snape >> $H/runs/minda.txt 2>&1; done
log "coverage"; bash $H/run.sh coverage
log "threads"; bash $H/run.sh threads
log "competitors"; bash $H/compete.sh
log "all done"
