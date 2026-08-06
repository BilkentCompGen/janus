#!/usr/bin/env bash
# Wait for the VAF sweep to release the cores, then run the compression ladder.
while pgrep -f vaf_test.sh > /dev/null 2>&1; do sleep 60; done
exec bash /home/yuio/svp/hifi/hpc_test.sh
