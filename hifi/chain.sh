#!/usr/bin/env bash
# Runs the pipeline for one LCP level and scores it. LV comes from the environment (default 5) and is
# passed straight through, so `LV=4 bash chain.sh` builds g4 and scores sev_ctg_L4.
#
# This exists as a file rather than an inline `bash -c` because launching detached work through
# `bash -c '...'` is refused here, and the run is far too long to hold a foreground shell.
set -uo pipefail
cd /home/yuio/svp/colo829
LV=${LV:-5}
for f in t.svb n.svb; do
    [ -s "$f" ] || { echo "[FATAL] $f yok -- zincir durdu"; exit 1; }
done
LV=$LV bash run.sh && LV=$LV bash score.sh
