#!/usr/bin/env bash
# Whole-genome 30x PacBio HiFi via pbsim3 multi-pass CLR + REAL ccs consensus.
# This is the realistic path (correlated, context-dependent HiFi errors from the
# actual ccs polisher), as opposed to the fast direct injector tools/sim_hifi.
#
#   Phase 1: per primary chromosome, multi-pass pbsim (parallel, capped) -> subreads
#            BAM; the MAF is redirected to /dev/null via a pre-created symlink so it is
#            never written to disk (ccs needs only the BAM; post-consensus errors aren't
#            recoverable from the subread MAF anyway), and .ref / region are truncated.
#   Phase 2: ccs each subreads BAM (all cores) -> HiFi FASTQ, appended; BAM truncated.
#   Output:  <out_dir>/hifi.fastq.gz  (+ <out_dir>/ccs_report.txt)
#
# pbsim/ccs flags mirror scripts/sim_reads.sh (the proven single-region pipeline).
# usage: sim_reads_wg_ccs.sh <input.fa> <out_dir>
# env:   DEPTH=30 LEN_MEAN=15000 LEN_SD=5000 PASS_NUM=10 DIFF_RATIO=22:45:33
#        METHOD=qshmm|errhmm  QSHMM=<model>  ERRHMM=<model>  PBSIM_PAR=4 CCS_THREADS=0

set -euo pipefail
INPUT="${1:?usage: sim_reads_wg_ccs.sh <input.fa> <out_dir>}"
OUTDIR="${2:?usage: sim_reads_wg_ccs.sh <input.fa> <out_dir>}"

DEPTH="${DEPTH:-30}"; LEN_MEAN="${LEN_MEAN:-15000}"; LEN_SD="${LEN_SD:-5000}"
PASS_NUM="${PASS_NUM:-10}"; DIFF_RATIO="${DIFF_RATIO:-22:45:33}"
METHOD="${METHOD:-qshmm}"; PBSIM_PAR="${PBSIM_PAR:-16}"; CCS_THREADS="${CCS_THREADS:-0}"; CCS_PAR="${CCS_PAR:-2}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PB="$ROOT/lib/pbsim3/src/pbsim"; CCS="$ROOT/lib/ccs/ccs"; S="${SAMTOOLS:-samtools}"; SHIM="$ROOT/scripts/fakebin"
QSHMM="${QSHMM:-$ROOT/lib/pbsim3/data/QSHMM-RSII.model}"
ERRHMM="${ERRHMM:-$ROOT/lib/pbsim3/data/ERRHMM-SEQUEL.model}"
[ -x "$PB" ] && [ -x "$CCS" ] && [ -f "$INPUT" ] || { echo "[ccs-wg] missing pbsim/ccs/input" >&2; exit 2; }
chmod +x "$SHIM/gzip" 2>/dev/null || true   # gzip->cat shim neutralises pbsim's discarded-MAF compression (which throttles pbsim ~3.7x via the popen pipe)

mkdir -p "$OUTDIR"; OUTDIR="$(cd "$OUTDIR" && pwd)"
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
[ -f "$INPUT.fai" ] || "$S" faidx "$INPUT"
mapfile -t CHROMS < <(awk '$1 ~ /^([0-9]+|X|Y)(_hap[AB])?$/ {print $2"\t"$1}' "$INPUT.fai" | sort -k1 -nr | cut -f2)
[ "${#CHROMS[@]}" -gt 0 ] || { echo "[ccs-wg] no primary chromosomes" >&2; exit 3; }
echo "[ccs-wg] $(date '+%F %T') input=$INPUT chroms=${#CHROMS[@]} depth=${DEPTH}x method=$METHOD pass=$PASS_NUM par=$PBSIM_PAR"

pbsim_one() {
    local chrom="$1" reg="$OUTDIR/$chrom.region.fa"
    "$S" faidx "$INPUT" "$chrom" > "$reg"
    ln -sf /dev/null "$OUTDIR/${chrom}_0001.maf.gz"   # discard ground-truth MAF: pbsim follows the symlink, 0 bytes hit disk
    if [ "$METHOD" = errhmm ]; then
        ( cd "$OUTDIR" && PATH="$SHIM:$PATH" "$PB" --strategy wgs --method errhmm --errhmm "$ERRHMM" --depth "$DEPTH" \
            --length-mean "$LEN_MEAN" --length-sd "$LEN_SD" --pass-num "$PASS_NUM" \
            --genome "$chrom.region.fa" --prefix "$chrom" > "$chrom.pbsim.log" 2>&1 )
    else
        ( cd "$OUTDIR" && PATH="$SHIM:$PATH" "$PB" --strategy wgs --method qshmm --qshmm "$QSHMM" --depth "$DEPTH" \
            --length-mean "$LEN_MEAN" --length-sd "$LEN_SD" --pass-num "$PASS_NUM" \
            --difference-ratio "$DIFF_RATIO" --genome "$chrom.region.fa" --prefix "$chrom" > "$chrom.pbsim.log" 2>&1 )
    fi
    for f in "$OUTDIR/${chrom}_"*.maf.gz "$OUTDIR/${chrom}_"*.ref "$reg"; do [ -f "$f" ] && : > "$f"; done
    echo "[ccs-wg] pbsim done: $chrom $(date '+%T')"
}

echo "[ccs-wg] === phase 1: multipass pbsim (parallel x$PBSIM_PAR) ==="
for chrom in "${CHROMS[@]}"; do
    pbsim_one "$chrom" &
    while [ "$(jobs -rp | wc -l)" -ge "$PBSIM_PAR" ]; do wait -n; done
done
wait
ccs_one() {
    local chrom="$1" bam="$OUTDIR/${chrom}_0001.bam"
    [ -f "$bam" ] || { echo "[ccs-wg] WARN missing BAM $chrom" >&2; return; }
    "$CCS" -j "$CCS_THREADS" --report-file "$OUTDIR/${chrom}.ccs_report.txt" "$bam" "$OUTDIR/${chrom}.hifi.fq.gz" >/dev/null 2>&1
    for f in "$bam" "$OUTDIR/${chrom}.hifi.zmw_metrics.json.gz"; do [ -f "$f" ] && : > "$f"; done   # free the subread BAM as soon as consumed
    echo "[ccs-wg] ccs done: $chrom $(date '+%T')"
}
echo "[ccs-wg] === phase 2: ccs (parallel x$CCS_PAR, -j $CCS_THREADS each; one ccs plateaus at ~20/32 cores, so 2 fill the box) ==="
: > "$OUTDIR/hifi.fastq.gz"; : > "$OUTDIR/ccs_report.txt"
for chrom in "${CHROMS[@]}"; do
    ccs_one "$chrom" &
    while [ "$(jobs -rp | wc -l)" -ge "$CCS_PAR" ]; do wait -n; done
done
wait
for chrom in "${CHROMS[@]}"; do   # concat sequentially: concurrent appends to one file would interleave/corrupt
    fq="$OUTDIR/${chrom}.hifi.fq.gz"
    [ -s "$fq" ] && { cat "$fq" >> "$OUTDIR/hifi.fastq.gz"; : > "$fq"; }
    [ -f "$OUTDIR/${chrom}.ccs_report.txt" ] && { echo "### $chrom"; cat "$OUTDIR/${chrom}.ccs_report.txt"; } >> "$OUTDIR/ccs_report.txt"
done
N=$(zcat "$OUTDIR/hifi.fastq.gz" | awk 'NR%4==2{n++; t+=length($0)} END{print n"\t"t}')
echo "[ccs-wg] $(date '+%F %T') DONE $OUTDIR/hifi.fastq.gz  reads/total_bp = $N"
