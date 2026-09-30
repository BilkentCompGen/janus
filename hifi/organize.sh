#!/usr/bin/env bash
# organize.sh: runs beside the fastq conversions (tmux session "organize"), every 5 minutes:
#  - moves each finished raw/<pair>_<n|t>.fastq.gz into its pair folder as <pair>/<n|t>.fastq.gz (no link is left: the
#    realignment that read from raw/ was stopped by the user on 2026-10-01), and unlinks the raw/ links earlier versions left;
#  - exits once all ten fastq.gz of the five pairs are in place.
H=/home/yuio/snape/hifi; RAW=$H/raw; L=$H/_aln2/organize.log
PAIRS="hcc1954 hg008 h2009 hcc1937 h1437"
log() { echo "$(date '+%F %T') $*" >> "$L"; }
while :; do
    for f in "$RAW"/*.fastq.gz; do
        [ -L "$f" ] && unlink "$f" && { log "unlinked raw/$(basename "$f")"; continue; }
        [ -f "$f" ] || continue
        b=$(basename "$f" .fastq.gz); p=${b%_*}; s=${b##*_}
        [ -d "$H/$p" ] && [[ $s == n || $s == t ]] || continue
        mv "$f" "$H/$p/$s.fastq.gz" && log "moved $b.fastq.gz to $p/$s.fastq.gz"
    done
    n=0; for p in $PAIRS; do for s in n t; do [ -s "$H/$p/$s.fastq.gz" ] && n=$((n + 1)); done; done
    [ $n -eq 10 ] && { log "all ten fastq.gz in place; done"; exit 0; }
    sleep 300
done
