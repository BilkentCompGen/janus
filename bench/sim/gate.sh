#!/usr/bin/env bash
L=/home/yuio/snape/bench/runs/night/queue.log; cur=thaw
while :; do
    ph=$(grep -E "phase:|nothing left|queue \(re\)started" $L | tail -1)
    case "$ph" in *alone*|*"one at a time"*) want=freeze ;; *) want=thaw ;; esac
    for u in snape-sim.service snape-simalign.service; do systemctl --user is-active -q $u && systemctl --user $want $u 2> /dev/null; done
    if [ $want != $cur ]; then
        echo "$(date '+%F %T') sim: $want (queue: ${ph#* * })" >> $L; cur=$want
    fi
    sleep 20
done
