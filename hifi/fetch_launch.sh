#!/bin/bash
# starts (or resumes) fetch_reads.sh in a detached tmux session "reads", so it survives the terminal and Claude session.
# Watch: tmux attach -t reads (detach with Ctrl-b d). Safe to run again: it does nothing while the session is alive.
tmux has-session -t reads 2>/dev/null && { echo "tmux session 'reads' is already running"; exit 0; }
mkdir -p /home/yuio/snape/hifi/raw
tmux new-session -d -s reads -c /home/yuio/snape/hifi "bash fetch_reads.sh 2>&1 | tee -a raw/fetch.log"
echo "started in tmux session 'reads'"
