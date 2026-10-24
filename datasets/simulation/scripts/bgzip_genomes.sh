#!/bin/bash
# Compress the simulation's genome FASTAs in place with bgzip. Only files >= 100 MB, so the small
# provenance files (sim.vcf, sim.novel_insertions.fa, *.bed) stay readable by name. bgzip removes
# the original once the .gz is written, so peak extra disk is one file.
set -u
f=$1
[ -s "$f" ] || exit 0
b=$(stat -c%s "$f")
bgzip -@ 8 "$f" || { echo "[FAIL] $f"; exit 1; }
a=$(stat -c%s "$f.gz")
echo "[OK] $f  $((b/1073741824)) -> $((a/1073741824)) GiB"
