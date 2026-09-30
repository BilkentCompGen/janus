#!/bin/bash
# run.sh <run dir> <snape binary> [first seed] [select root with <pair>/ev.{tsv,paf,fa}]: call on the five pairs' select outputs, held-out verdict, fair evaluation.
# Writes <run dir>/<pair>/{calls,verdict}.vcf, evaluation files in bench/out/, and one summary line to <run dir>/summary.txt.
d=$(realpath -m $1); B=$2; s=${3:-0}; A=${4:-/home/yuio/tmp/snape-castle/accepted3/castle-bench}; H=/home/yuio/snape/bench; mkdir -p $d
for p in hcc1954 hcc1937 hg008 h1437 h2009; do mkdir -p $d/$p; (cd $A/$p && nice -n 5 $B call -i ev.paf -r ev.tsv -d /home/datasets/references/genomes/hg38/hg38.fa -o $d/$p/calls.vcf -v /home/yuio/snape/ref/vntr.trf.bed > $d/$p/call.log 2>&1) & done; wait
OMP_NUM_THREADS=8 nice /home/yuio/snape/.venv-ml/bin/python $H/verdict.py $d $s > $d/verdict.log 2>&1
t=$(basename $d)_s$s
for p in hcc1954 hcc1937 hg008 h1437 h2009; do nice $H/eval.sh $p $d/$p/verdict.vcf $t > $d/$p/eval.txt 2>&1 & done; wait
cat $d/*/eval.txt | grep -E 'A_fair|A_nist' | awk -v t=$t '{f[$1]=$5} END {printf "%s hcc1954 %.4f hg008 %.4f h2009 %.4f h1437 %.4f hcc1937 %.4f mean %.4f\n", t, f["hcc1954"], f["hg008"], f["h2009"], f["h1437"], f["hcc1937"], (f["hcc1954"]+f["hg008"]+f["h2009"]+f["h1437"]+f["hcc1937"])/5}' | tee -a $d/summary.txt
