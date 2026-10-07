#!/usr/bin/env bash
H=/home/yuio/snape/bench; L=$H/runs/all.log; V=/home/yuio/snape/.venv-ml/bin/python; mkdir -p $H/runs $H/verdict
log() { echo "$(date '+%F %T') $*" >> $L; }
log "build"; make -C /home/yuio/snape -j16 > $H/runs/build.log 2>&1 || { log "build FAILED"; exit 1; }
log "full runs"; bash $H/run.sh full
log "verdict"; ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 $V verdict.py > verdict/stdout.txt 2>&1 )
log "minda"; for p in hcc1954 h2009 hcc1937 h1437 hcc1395 hg008; do taskset -c 0-15 bash $H/eval.sh $p $H/verdict/$p.vcf snape >> $H/runs/minda.txt 2>&1; done
log "coverage"; bash $H/run.sh coverage
log "coverage verdict"; ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 $V cov_verdict_matched.py > verdict/coverage_stdout.txt 2>&1 )
log "purity"; bash $H/run.sh purity
log "purity verdict"; ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 $V purity_verdict.py > verdict/purity_stdout.txt 2>&1 )
log "threads"; bash $H/run.sh threads
log "download hs578t"; bash /home/yuio/snape/hifi/fetch.sh hs578t > $H/runs/fetch_hs578t.log 2>&1
log "full hs578t"; [ -s /home/yuio/snape/hifi/hs578t/t.fastq.gz ] && bash $H/run.sh full hs578t
log "competitors"; bash $H/compete.sh
log "all done"
