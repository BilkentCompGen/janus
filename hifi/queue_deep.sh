#!/usr/bin/env bash
# Deep normals for the other three, queued behind whatever is already running. Waits on the two
# specific artefacts rather than on "is any minimap2 running", because that test blocks on this
# script's own alignments and deadlocks the queue.
#
# 2M extra reads each rather than 6M: the depth response measured on hcc1954 (3x -> 290 FP,
# 5.5x -> 172, 11x -> 136) is steepest at the low end, so doubling to ~22x takes most of what is
# left, at a third of the wall clock.
set -uo pipefail
cd /home/yuio/svp
FINAL=/tmp/claude-1000/-home-claude/03a7c320-c0bf-47b5-8a15-2f7c6c653738/tasks/b9op76izn.output

until grep -q "########## bitti" "$FINAL" 2>/dev/null; do sleep 30; done
until [ -s /home/yuio/tmp/jx/hcc1954/normal_deep.bam.bai ]; do sleep 30; done
echo "[kuyruk] onceki isler bitti, derin normaller basliyor $(date +%T)"

for S in hcc1937 h2009 h1437; do
    EXTRA=2000000 bash hifi/deep_normal.sh $S
done
echo "[kuyruk] hepsi bitti $(date +%T)"
