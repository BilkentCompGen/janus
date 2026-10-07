#!/usr/bin/env bash
H=/home/yuio/snape/bench; L=$H/runs/all.log; V=/home/yuio/snape/.venv-ml/bin/python
log() { echo "$(date '+%F %T') $*" >> $L; }
log "germline branch: full runs of the six truth pairs"; bash $H/run.sh full hcc1954 h2009 hcc1937 h1437 hcc1395 hg008
log "germline branch: verdict (57 features)"; ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 $V verdict.py > verdict/stdout.txt 2>&1 )
log "germline branch: minda"; : > $H/runs/minda.txt; for p in hcc1954 h2009 hcc1937 h1437 hcc1395 hg008; do taskset -c 0-15 bash $H/eval.sh $p $H/verdict/$p.vcf snape >> $H/runs/minda.txt 2>&1; done
log "germline branch rerun done"
