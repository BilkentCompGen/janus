#!/usr/bin/env bash
cd /home/yuio/snape
D=bench/runs/trace/hcc1395_20x
python3 bench/lowcov.py hcc1395 $D bench/verdict/coverage_matched/hcc1395_20x.vcf > $D/lowcov.out 2>&1
echo "exit $?" >> $D/lowcov.out
