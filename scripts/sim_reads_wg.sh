#!/usr/bin/env bash
# sim_reads_wg.sh — whole-genome modern PacBio-HiFi reads via tools/sim_hifi.
#
# pbsim cannot produce Q30 in a single pass (CLR models cap ~Q14, measured 0.958)
# and its multipass->ccs HiFi path costs ~1160 CPU-hours/genome (~3 days for a
# tumor+normal pair). sim_hifi models the Q30 HiFi *output* directly (sample
# fragments + inject errors at 0.1% with the Sequel sub:ins:del=22:45:33 profile
# and homopolymer deletion bias), so there is no consensus stage: each primary
# chromosome is an independent job and they run fully in parallel, and every
# injected error is written to an exact ground-truth file.
#
# Usage:  sim_reads_wg.sh <input.fa> <out_dir>
# Env:  DEPTH=30 LEN_MEAN=15000 LEN_SD=5000 ERR=0.001 RATIO=22:45:33
#       HP_DEL_BIAS=3 PAR=24 SEED=42

set -euo pipefail

INPUT="${1:?usage: sim_reads_wg.sh <input.fa> <out_dir>}"
OUTDIR="${2:?usage: sim_reads_wg.sh <input.fa> <out_dir>}"

DEPTH="${DEPTH:-30}"; LEN_MEAN="${LEN_MEAN:-15000}"; LEN_SD="${LEN_SD:-5000}"
ERR="${ERR:-0.001}"; RATIO="${RATIO:-22:45:33}"; HP_DEL_BIAS="${HP_DEL_BIAS:-3}"
PAR="${PAR:-24}"; SEED="${SEED:-42}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SIM="$ROOT/tools/sim_hifi"
SAMTOOLS="${SAMTOOLS:-samtools}"
[ -x "$SIM" ] || { echo "[wg] sim_hifi missing: $SIM (run: make tools/sim_hifi)" >&2; exit 2; }
[ -f "$INPUT" ] || { echo "[wg] input not found: $INPUT" >&2; exit 2; }

mkdir -p "$OUTDIR" "$OUTDIR/errors"
OUTDIR="$(cd "$OUTDIR" && pwd)"
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
[ -f "$INPUT.fai" ] || "$SAMTOOLS" faidx "$INPUT"

mapfile -t CHROMS < <(awk '$1 ~ /^([0-9]+|X|Y)(_hap[AB])?$/ {print $2"\t"$1}' "$INPUT.fai" | sort -k1 -nr | cut -f2)
[ "${#CHROMS[@]}" -gt 0 ] || { echo "[wg] no primary chromosomes in $INPUT.fai" >&2; exit 3; }
echo "[wg] $(date '+%F %T') input=$INPUT chroms=${#CHROMS[@]} depth=${DEPTH}x err=$ERR par=$PAR"

sim_one() {
    local chrom="$1" reg="$OUTDIR/$chrom.region.fa" idx
    "$SAMTOOLS" faidx "$INPUT" "$chrom" > "$reg"
    # deterministic per-chrom seed so the whole run is reproducible
    idx=$(printf '%s' "$chrom" | cksum | cut -d' ' -f1)
    "$SIM" "$reg" "$OUTDIR/$chrom.fq.gz" --errors "$OUTDIR/errors/$chrom.errors.tsv.gz" \
        --depth "$DEPTH" --len-mean "$LEN_MEAN" --len-sd "$LEN_SD" \
        --err "$ERR" --ratio "$RATIO" --hp-del-bias "$HP_DEL_BIAS" \
        --seed "$((SEED + idx))" > "$OUTDIR/$chrom.sim.log" 2>&1
    : > "$reg"
    echo "[wg] done: $chrom $(date '+%T')"
}

for chrom in "${CHROMS[@]}"; do
    sim_one "$chrom" &
    while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do wait -n; done
done
wait
echo "[wg] === concatenating FASTQ ==="
: > "$OUTDIR/hifi.fastq.gz"
for chrom in "${CHROMS[@]}"; do
    f="$OUTDIR/$chrom.fq.gz"
    [ -f "$f" ] && { cat "$f" >> "$OUTDIR/hifi.fastq.gz"; : > "$f"; } || echo "[wg] WARN missing fq for $chrom" >&2
done
N=$(zcat "$OUTDIR/hifi.fastq.gz" | awk 'NR%4==2{n++; t+=length($0)} END{print n"\t"t}')
echo "[wg] $(date '+%F %T') DONE $OUTDIR/hifi.fastq.gz  reads/total_bp = $N  (errors: $OUTDIR/errors/*.errors.tsv.gz)"
