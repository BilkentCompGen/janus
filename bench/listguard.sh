#!/usr/bin/env bash
N=/home/yuio/snape/bench/runs/night
while true; do
    ok=" $(grep -hv '^#' $N/fast.txt $N/slow.txt $N/gpu.txt 2>/dev/null | tr '\n' ' ') "
    for f in $N/running/*.pgid; do
        [ -e $f ] || continue
        t=$(basename $f | cut -d. -f1); p=$(basename $f | cut -d. -f2); g=$(cat $f)
        kill -0 -- -$g 2>/dev/null || continue
        case "$ok" in *" $t "*) continue;; esac
        for c in $(docker ps -q); do docker inspect -f '{{range .Mounts}}{{.Source}} {{end}}' $c | grep -q "/competitors/$t/$p" && docker kill $c > /dev/null 2>&1; done
        kill -TERM -- -$g 2>/dev/null
        echo "$(date '+%F %T') LIST GUARD: $t/$p is not in the tool lists (lane still on an old list); killed" >> $N/INCIDENTS.md
    done
    sleep 60
done
