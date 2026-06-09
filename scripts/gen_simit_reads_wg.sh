#!/usr/bin/env bash
# gen_simit_reads_wg.sh — whole-genome HiFi reads for the somatic benchmark via Sim-it.
#   normal = simulation/ref.fasta        -> reads/normal_simit/hifi.fastq.gz
#   tumor  = simulation/truth/wgs/out.fa -> reads/wgs_simit/hifi.fastq.gz
# Resumable: a sample is skipped if its hifi.fastq.gz already exists.
# Usage: gen_simit_reads_wg.sh [depth=30] [threads=16]
set -euo pipefail
ROOT=/home/yuio/svp
DEPTH="${1:-30}"
THREADS="${2:-16}"
REF="$ROOT/simulation/ref.fasta"
TUMOR_FA="$ROOT/simulation/truth/wgs/out.fa"
RN="$ROOT/simulation/reads/normal_simit"
RT="$ROOT/simulation/reads/wgs_simit"

echo "[gen $(date +%F_%T)] START depth=$DEPTH threads=$THREADS"
if [ -s "$RN/hifi.fastq.gz" ]; then
  echo "[gen $(date +%T)] normal reads already exist -> skip"
else
  echo "[gen $(date +%T)] normal reads from ref.fasta"
  bash "$ROOT/scripts/sim_reads_simit.sh" -i "$REF" -o "$RN" -n normal -d "$DEPTH" -t "$THREADS"
fi
if [ -s "$RT/hifi.fastq.gz" ]; then
  echo "[gen $(date +%T)] tumor reads already exist -> skip"
else
  echo "[gen $(date +%T)] tumor reads from truth/wgs/out.fa"
  bash "$ROOT/scripts/sim_reads_simit.sh" -i "$TUMOR_FA" -o "$RT" -n tumor -d "$DEPTH" -t "$THREADS"
fi
echo "[gen $(date +%F_%T)] DONE"
echo "  normal: $RN/hifi.fastq.gz"
echo "  tumor : $RT/hifi.fastq.gz"
