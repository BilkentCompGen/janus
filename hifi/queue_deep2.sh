#!/usr/bin/env bash
# The remaining deep normals, queued behind h1437's. Measured on hcc1954, going from 11x to 44x took
# the false calls from 136 to 116 and F1 from 0.867 to 0.874 -- past Severus's 0.871 on that sample.
# h2009 and hcc1937 have a wider precision gap than hcc1954 did, so the same question asked with the
# same evidence should be worth at least as much there.
set -uo pipefail
cd /home/yuio/svp
until [ -s /home/yuio/tmp/jx/h1437/normal_deep.bam.bai ]; do sleep 30; done
for S in h2009 hcc1937; do
    EXTRA=4000000 bash hifi/deep_normal.sh $S
done
echo "[kuyruk] derin normaller bitti $(date +%T)"
