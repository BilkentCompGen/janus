#!/usr/bin/env bash
# purge_sources.sh (tmux session "purge"): every 2 minutes deletes the source files (.sra, unaligned .bam) that the fetch
# scripts have moved to ~/trash/reads-converted/ after their fastq.gz passed both checks (download md5/ETag, and the
# fastq.gz read count equal to the source's, with a gzip integrity test). The user approved deleting such verified
# sources on 2026-10-01. Files that failed a check go to reads-converted/bad/ and are never touched here.
D=/home/yuio/trash/reads-converted; L=/home/yuio/snape/hifi/_aln2/purge.log
while :; do
    for f in "$D"/*; do
        [ -f "$f" ] || continue
        s=$(du -BG --apparent-size "$f" | cut -f1)
        rm -f -- "$f" && echo "$(date '+%F %T') deleted verified source $(basename "$f") ($s)" >> "$L"
    done
    sleep 120
done
