#!/bin/bash
S=/home/yuio/snape; T=/home/yuio/trash/2026-10-03-clean; mkdir -p $T/bench $T/hifi
mv $S/tmp $S/scratch $S/build $T/
for d in runs eval verdict __pycache__; do [ -e $S/bench/$d ] && mv $S/bench/$d $T/bench/; done
[ -e /tmp/claude-1000/-home-yuio-snape ] && mv /tmp/claude-1000/-home-yuio-snape $T/claude_scratch
mkdir -p $S/bench/runs $S/bench/verdict
du -sh $T; df -h /home | tail -1
