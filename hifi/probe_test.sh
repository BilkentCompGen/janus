#!/usr/bin/env bash
# Are the "somatic" cores under the false positives actually absent from the normal reads?
# Two sets from the same sample, same graph, same normal scan: the contigs that align with more than
# 2% divergence (which is where the false positives come from) and a size-matched control of contigs
# that align cleanly everywhere.
set -uo pipefail
cd /home/yuio/svp
S=hcc1954
JX=/home/yuio/tmp/jx/$S
for SET in dirty clean; do
    echo "########## $SET  $(date +%T) ##########"
    taskset -c 0-15 tools/normal_probe outputs/svp/$S/graph.bin "$JX/$SET.fa" \
        <(./svp dump hifi/$S/n.svb) "$JX/probe_$SET.tsv" 16 2>&1 | /usr/bin/grep -vE "^\[probe\]   [0-9]+ okuma"
done
echo "########## bitti $(date +%T) ##########"
