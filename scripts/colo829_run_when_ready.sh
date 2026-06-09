#!/usr/bin/env bash
# colo829_run_when_ready.sh — auto-continuation: wait for the COLO829 download+convert to
# finish, then run the svp caller and the COLO829 truvari benchmark, unattended.
#
# Waits for BOTH FASTQs to exist AND the download screen's DONE marker (so we don't start
# on a half-converted file). Resumable: svp run skips completed stages; benchmark re-runs cheap.
# Usage: scripts/colo829_run_when_ready.sh [min_mult=10] [threads=16]
set -euo pipefail
ROOT=/home/yuio/svp
MM="${1:-10}"
THREADS="${2:-16}"
DL_LOG="$ROOT/colo829/download.log"
TUMOR="$ROOT/colo829/tumor/COLO829_tumor_hifi.fastq.gz"
NORMAL="$ROOT/colo829/normal/COLO829BL_normal_hifi.fastq.gz"
OUT="$ROOT/outputs/svp/colo829"
log(){ echo "[colo829-run $(date +%F_%T)] $*"; }

log "WAIT: for reads + download DONE marker (poll 120s)"
while :; do
    if [ -s "$TUMOR" ] && [ -s "$NORMAL" ] && grep -q "] DONE" "$DL_LOG" 2>/dev/null; then
        log "reads present and download DONE -> proceed"; break
    fi
    sleep 120
done

# guard: files must be stable (not still being written) — two size reads 30s apart match
stable(){ local a b; a=$(stat -c %s "$1"); sleep 30; b=$(stat -c %s "$1"); [ "$a" = "$b" ]; }
while ! { stable "$TUMOR" && stable "$NORMAL"; }; do log "reads still growing; wait"; sleep 30; done

log "reads: tumor=$(du -h "$TUMOR"|cut -f1) normal=$(du -h "$NORMAL"|cut -f1)"
mkdir -p "$OUT"

log "svp run (mm=$MM, $THREADS threads) -> $OUT  [normal=COLO829-BL, tumor=COLO829]"
t0=$SECONDS
"$ROOT/svp" run "$NORMAL" "$TUMOR" "$OUT" --ref "$ROOT/simulation/ref.fasta" --min-mult "$MM" --threads "$THREADS"
log "svp run done in $((SECONDS-t0))s; calls.vcf=$(grep -vc '^#' "$OUT/calls.vcf" 2>/dev/null || echo NA) records"

log "benchmark vs COLO829 truth"
bash "$ROOT/scripts/colo829_benchmark.sh" "$OUT" 2>&1 | tee "$OUT/benchmark.out"
log "DONE — results in $OUT (truvari_geq50/, benchmark.out)"
