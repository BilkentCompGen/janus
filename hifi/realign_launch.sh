#!/bin/bash
# starts (or resumes) realign.sh in a detached tmux session "align"; does nothing while the session is alive
tmux has-session -t align 2>/dev/null && { echo "tmux session 'align' is already running"; exit 0; }
mkdir -p /home/yuio/snape/hifi/_aln2
tmux new-session -d -s align -c /home/yuio/snape/hifi "bash realign.sh 2>&1 | tee -a _aln2/realign.log"
echo "started in tmux session 'align'"
