#!/usr/bin/env bash
N=/home/yuio/snape/bench/runs/night; mkdir -p $N/running
while :; do
    a=$(awk '/MemAvailable/{print int($2 / 1048576)}' /proc/meminfo)
    if [ "$a" -lt 14 ]; then
        {
            echo "$(date '+%F %T') memguard: only $a GB of memory available; killing the running competitor jobs"
            ps -eo rss,pid,args --sort=-rss | head -8 | cut -c1-170 | sed 's/^/    /'
            docker ps --format '    container {{.ID}} {{.Image}} {{.RunningFor}}'
            for f in $N/running/*.pgid; do [ -s "$f" ] && echo "    job $(basename $f .pgid)"; done
        } >> $N/INCIDENTS.md
        c=$(docker ps -q); [ -n "$c" ] && docker kill $c > /dev/null 2>&1
        for f in $N/running/*.pgid; do [ -s "$f" ] && pkill -9 -s $(cat $f) 2> /dev/null; done
        sleep 90
    fi
    sleep 4
done
