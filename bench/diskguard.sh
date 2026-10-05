#!/usr/bin/env bash
N=/home/yuio/snape/bench/runs/night; held=0
while true; do
    a=$(df --output=avail -BG /home | tail -1 | tr -dc 0-9)
    if [ $held = 0 ] && [ $a -lt 60 ]; then
        mine=$(docker ps -q --filter status=running); [ -n "$mine" ] && docker pause $mine > /dev/null
        systemctl --user freeze snape-night.service snape-seq.service snape-rerun.service 2>/dev/null
        echo "$(date '+%F %T') DISK GUARD: $a GB free on /home; queue, experiments and containers held until 120 GB are free" >> $N/INCIDENTS.md
        held=1
    elif [ $held = 1 ] && [ $a -ge 120 ]; then
        [ -n "$mine" ] && docker unpause $mine > /dev/null
        systemctl --user thaw snape-night.service snape-seq.service snape-rerun.service 2>/dev/null
        echo "$(date '+%F %T') DISK GUARD: $a GB free again, everything resumed" >> $N/INCIDENTS.md
        held=0
    fi
    sleep 60
done
