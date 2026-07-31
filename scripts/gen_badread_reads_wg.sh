#!/usr/bin/env bash
# gen_badread_reads_wg.sh — whole-genome HiFi reads for the somatic benchmark via Badread.
#   normal = simulation/ref.fasta            -> simulation/reads/normal/hifi.fastq
#   tumor  = simulation/truth/<id>/out.fa    -> simulation/reads/<id>/hifi.fastq  (id default: hard)
# Regular-difficulty reads: Q30 / stdev 3 on the pacbio2021 model (set in sim_reads_badread.sh).
# Uses all cores by default (Badread is single-threaded; the script shards contigs across threads).
# Resumable: a sample is skipped if its hifi.fastq already exists.
# Usage: gen_badread_reads_wg.sh [sim_id=hard] [depth=30] [threads=$(nproc)]
set -euo pipefail
ROOT=/home/yuio/svp
ID="${1:-hard}"
DEPTH="${2:-30}"
THREADS="${3:-$(nproc)}"
REF="$ROOT/simulation/ref.fasta"
TUMOR_FA="$ROOT/simulation/truth/$ID/out.fa"
RN="$ROOT/simulation/reads/normal"
RT="$ROOT/simulation/reads/$ID"

echo "[gen $(date +%F_%T)] START id=$ID depth=$DEPTH threads=$THREADS"
[ -f "$TUMOR_FA" ] || { echo "[gen] tumor genome missing: $TUMOR_FA (run 'make -C simulation $ID')" >&2; exit 2; }
if [ -s "$RN/hifi.fastq" ]; then
  echo "[gen $(date +%T)] normal reads already exist -> skip"
else
  echo "[gen $(date +%T)] normal reads from ref.fasta"
  bash "$ROOT/scripts/sim_reads_badread.sh" -i "$REF" -o "$RN" -d "$DEPTH" -t "$THREADS" -s 1
fi
if [ -s "$RT/hifi.fastq" ]; then
  echo "[gen $(date +%T)] tumor reads already exist -> skip"
else
  echo "[gen $(date +%T)] tumor reads from truth/$ID/out.fa"
  bash "$ROOT/scripts/sim_reads_badread.sh" -i "$TUMOR_FA" -o "$RT" -d "$DEPTH" -t "$THREADS" -s 1001
fi
echo "[gen $(date +%F_%T)] DONE"
echo "  normal: $RN/hifi.fastq"
echo "  tumor : $RT/hifi.fastq"
