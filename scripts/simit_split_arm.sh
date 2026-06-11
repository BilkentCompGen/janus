#!/usr/bin/env bash
# Sim-it per-color-gate (split) arm vs the max-keep baseline.
# Resumes from the existing wgs_simit pre_prune.bin (skips the 21-min pool/sort/merge),
# re-prunes with --prune-mode split (mm=10), runs the full pipeline, benchmarks vs truth.
set -uo pipefail
cd /home/yuio/svp
SRC=outputs/svp/wgs_simit
DST=outputs/svp/simit_split
mkdir -p "$DST"
[ -s "$DST/pre_prune.bin" ] || { echo "seeding pre_prune.bin"; cp "$SRC/pre_prune.bin" "$DST/pre_prune.bin"; }
echo "=== split arm: svp run (resume from pre_prune, --prune-mode split, mm=10) ==="
/usr/bin/time -v ./svp run \
  simulation/reads/normal_simit/hifi.fastq.gz \
  simulation/reads/wgs_simit/hifi.fastq.gz \
  "$DST" --ref simulation/ref.fasta --min-mult 10 --prune-mode split --threads 16 > "$DST/run.log" 2>&1
echo "run exit=$?  calls=$(grep -vc '^#' "$DST/calls.vcf" 2>/dev/null || echo NA)"
grep -E 'prune min_mult|kept=|Maximum resident' "$DST/run.log" | tail -3
echo "=== benchmark split arm vs truth ==="
bash scripts/svp_benchmark.sh "$DST" > "$DST/benchmark.out" 2>&1
echo "DONE — see $DST/benchmark.out"
