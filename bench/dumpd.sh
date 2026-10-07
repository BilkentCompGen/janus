#!/usr/bin/env bash
p=$1; D=$2; S=/home/yuio/snape; B=$S/build/snape; P=16; mkdir -p $D; cd $D
for s in n t; do [ -s $s.svb ] || taskset -c 0-15 $B pack -i $S/hifi/$p/$s.fastq.gz -o $s.svb -t $P; done
for s in n t; do [ -s $s.svb.exc ] || taskset -c 0-15 $B extract -i $s.svb -o $s.svb.exc -t $P; done
for s in n t; do [ -s $s.cnt ] || taskset -c 0-15 $B count -i $s.svb -o $s.cnt -t $P; done
[ -s graph.gph ] || taskset -c 0-15 $B graph -i n.cnt,t.cnt -o graph.gph
[ -s prune.sites ] || taskset -c 0-15 $B prune -i graph.gph -r n.svb,t.svb -o prune -t $P
SNAPE_DUMP_D=$D/d.tsv taskset -c 0-15 $B select -i graph.gph -d prune.sites -r n.svb,t.svb -o ev -t $P
