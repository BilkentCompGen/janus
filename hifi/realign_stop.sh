#!/bin/bash
# stops realign.sh and its disk guard (tmux sessions align, diskguard) with the running minimap2 / samtools sort
tmux kill-session -t align 2>/dev/null; tmux kill-session -t diskguard 2>/dev/null
pkill -CONT -f "minimap2 -ax map-hifi"; pkill -f "minimap2 -ax map-hifi"; pkill -f "samtools sort -@ 4 -m 2G -T /home/yuio/snape/hifi/_aln2"
pkill -f "bash realign.sh"; pkill -f "diskguard.sh"; sleep 3; pkill -9 -f "minimap2 -ax map-hifi"; true
