#!/bin/bash
# stops fetch_reads.sh with its transfers and conversions; downloads resume from the .aria2 control files, an unfinished
# conversion restarts on the next fetch_launch.sh
pkill -f "bash fetch_reads.sh"; pkill -x aria2c; pkill -f "fasterq-dump"; pkill -f "samtools fastq"; sleep 3
pkill -9 -x aria2c; pkill -9 -f fasterq-dump; tmux kill-session -t reads 2>/dev/null; true
