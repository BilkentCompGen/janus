#!/usr/bin/env bash
H=/home/yuio/snape/bench; L=$H/runs/all.log
log() { echo "$(date '+%F %T') $*" >> $L; }
log "alignments only (competitors on hold)"; bash $H/align_only.sh
log "all done"
