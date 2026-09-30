#!/bin/bash
# diskguard.sh: every 30 s, pauses (SIGSTOP) the realignment (minimap2 + samtools sort of realign.sh) when free space on
# /home falls below LOW GB, and resumes it (SIGCONT) above HIGH GB, so a full disk never breaks it. Logs to _aln2/diskguard.log.
LOW=30; HIGH=60; st=run; L=/home/yuio/snape/hifi/_aln2/diskguard.log
pids() { pgrep -f "minimap2 -ax map-hifi .*hifi/raw/"; pgrep -f "samtools sort .*_aln2/"; }
while :; do
  f=$(df -BG --output=avail /home | tail -1 | tr -dc 0-9)
  if [ $st = run ] && [ $f -lt $LOW ]; then kill -STOP $(pids) 2>/dev/null; st=paused; echo "$(date '+%F %T') free ${f}G: alignment paused" >> $L; fi
  if [ $st = paused ] && [ $f -gt $HIGH ]; then kill -CONT $(pids) 2>/dev/null; st=run; echo "$(date '+%F %T') free ${f}G: alignment resumed" >> $L; fi
  sleep 30
done
