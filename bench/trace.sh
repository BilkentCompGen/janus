#!/usr/bin/env bash
S=/home/yuio/snape; H=$S/bench; L=$H/runs/all.log; G=3100000000
log() { echo "$(date '+%F %T') $*" >> $L; }
until grep -q "after: done" $L; do sleep 300; done
for p in hcc1954 hcc1395; do
    D=$H/runs/trace/${p}_20x; mkdir -p $D
    for f in $S/hifi/$p/n.fastq.gz $S/hifi/$p/t.fastq.gz; do dd if=$f iflag=nocache count=0 status=none; done
    log "trace: $p 20x start (intermediates kept)"
    BASES=$((20 * G)) MINREADS=3 taskset -c 0-15 bash $S/pipeline.sh $D $S/hifi/$p/n.fastq.gz $S/hifi/$p/t.fastq.gz 16 > $D/run.out 2>&1 \
        && log "trace: $p 20x done" || log "trace: $p 20x FAILED"
done
log "trace: done"
