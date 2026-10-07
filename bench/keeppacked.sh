#!/usr/bin/env bash
R=/home/yuio/snape/bench/runs; K=/home/yuio/snape/hifi
while true; do
    left=0
    for p in hcc1954 h2009 hcc1937 h1437 hcc1395 hg008; do
        [ -s $K/$p/t.svb.exc ] || [ -e $K/$p/svb.skipped ] && continue
        W=$R/full/$p/work
        if [ -e $W/done.extract_t ] && [ ! -e $K/$p/svb.skipped ]; then
            free=$(df --output=avail -BG /home | tail -1 | tr -dc 0-9)
            if [ $free -lt 300 ]; then
                touch $K/$p/svb.skipped
                echo "$(date '+%F %T') packed reads of $p not kept: only $free GB free" >> $R/night/queue.log
                continue
            fi
            mkdir -p $K/$p
            for f in n.svb n.svb.idx n.svb.exc t.svb t.svb.idx t.svb.exc; do ln -f $W/$f $K/$p/$f; done
            echo "$(date '+%F %T') kept packed reads of $p in $K/$p" >> $R/night/queue.log
        else
            left=1
        fi
    done
    [ $left = 0 ] && break
    sleep 20
done
