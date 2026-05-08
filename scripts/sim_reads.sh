#!/usr/bin/env bash
# sim_reads.sh — pbsim3 multi-pass CLR + ccs HiFi consensus.
#
# Usage:
#   sim_reads.sh -i <input.fa> -o <out_dir> [-r <region>] [-d depth]
#                [-l len_mean] [-s len_sd] [-a accuracy] [-p pass_num]
#                [-x diff_ratio] [-t threads] [--keep-bam]
#
# Defaults match realistic PacBio Sequel HiFi at 30x:
#   depth=30  len_mean=15000  len_sd=5000  accuracy=0.85
#   pass_num=10  diff_ratio=22:45:33  threads=0 (auto)
#
# Output layout:
#   <out_dir>/hifi.fastq.gz       — final HiFi reads
#   <out_dir>/pbsim.maf.gz        — pbsim ground-truth alignments
#   <out_dir>/ccs_report.txt      — ccs filter report
#   <out_dir>/pbsim.log           — pbsim console log
#   <out_dir>/subreads.bam        — kept only with --keep-bam

set -euo pipefail

INPUT=""
OUTDIR=""
REGION=""
DEPTH=30
LEN_MEAN=15000
LEN_SD=5000
ACCURACY=0.85
PASS_NUM=10
DIFF_RATIO="22:45:33"
THREADS=0
KEEP_BAM=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        -i) INPUT="$2"; shift 2;;
        -o) OUTDIR="$2"; shift 2;;
        -r) REGION="$2"; shift 2;;
        -d) DEPTH="$2"; shift 2;;
        -l) LEN_MEAN="$2"; shift 2;;
        -s) LEN_SD="$2"; shift 2;;
        -a) ACCURACY="$2"; shift 2;;
        -p) PASS_NUM="$2"; shift 2;;
        -x) DIFF_RATIO="$2"; shift 2;;
        -t) THREADS="$2"; shift 2;;
        --keep-bam) KEEP_BAM=1; shift;;
        -h|--help) sed -n '2,16p' "$0"; exit 0;;
        *) echo "[sim_reads] unknown arg: $1" >&2; exit 2;;
    esac
done

[ -n "$INPUT" ]  || { echo "[sim_reads] -i <input.fa> required" >&2; exit 2; }
[ -n "$OUTDIR" ] || { echo "[sim_reads] -o <out_dir> required" >&2; exit 2; }
[ -f "$INPUT" ]  || { echo "[sim_reads] input not found: $INPUT" >&2; exit 2; }

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PBSIM="$ROOT/lib/pbsim3/src/pbsim"
PBSIM_MODEL="$ROOT/lib/pbsim3/data/QSHMM-RSII.model"
CCS="$ROOT/lib/ccs/ccs"
SAMTOOLS="${SAMTOOLS:-samtools}"

[ -x "$PBSIM" ] || { echo "[sim_reads] pbsim missing: $PBSIM (run 'make install-deps')" >&2; exit 2; }
[ -x "$CCS" ]   || { echo "[sim_reads] ccs missing: $CCS (run 'make install-deps')" >&2; exit 2; }
[ -f "$PBSIM_MODEL" ] || { echo "[sim_reads] pbsim model missing: $PBSIM_MODEL" >&2; exit 2; }

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"
INPUT_ABS="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"

# 1. extract region (or copy whole file).
GENOME="$OUTDIR/region.fa"
if [ -n "$REGION" ]; then
    [ -f "$INPUT_ABS.fai" ] || "$SAMTOOLS" faidx "$INPUT_ABS"
    echo "[sim_reads] extracting region '$REGION' from $INPUT_ABS"
    "$SAMTOOLS" faidx "$INPUT_ABS" "$REGION" > "$GENOME"
    [ -s "$GENOME" ] || { echo "[sim_reads] empty extract — region '$REGION' not in $INPUT_ABS" >&2; exit 3; }
else
    echo "[sim_reads] no region given — using whole FASTA (heavy)"
    cp "$INPUT_ABS" "$GENOME"
fi

# 2. pbsim multi-pass — produces sd_NNNN.bam (one per FASTA record).
echo "[sim_reads] pbsim3: depth=$DEPTH len=$LEN_MEAN±$LEN_SD acc=$ACCURACY passes=$PASS_NUM diff=$DIFF_RATIO"
( cd "$OUTDIR" && \
    "$PBSIM" --strategy wgs --method qshmm \
        --qshmm "$PBSIM_MODEL" \
        --depth "$DEPTH" \
        --length-mean "$LEN_MEAN" --length-sd "$LEN_SD" \
        --accuracy-mean "$ACCURACY" \
        --pass-num "$PASS_NUM" \
        --difference-ratio "$DIFF_RATIO" \
        --genome region.fa \
        --prefix sd \
        > pbsim.log 2>&1 )

# 3. ccs each subreads BAM into HiFi FASTQ; concat into hifi.fastq.gz.
shopt -s nullglob
BAMS=( "$OUTDIR"/sd_*.bam )
[ "${#BAMS[@]}" -gt 0 ] || { echo "[sim_reads] pbsim produced no BAM — see $OUTDIR/pbsim.log" >&2; exit 4; }

echo "[sim_reads] ccs: ${#BAMS[@]} subreads BAM → HiFi"
: > "$OUTDIR/hifi.fastq.gz"   # truncate
for bam in "${BAMS[@]}"; do
    stem="${bam%.bam}.hifi"
    "$CCS" -j "$THREADS" --report-file "$OUTDIR/ccs_report.txt" "$bam" "$stem.fastq.gz"
    cat "$stem.fastq.gz" >> "$OUTDIR/hifi.fastq.gz"
    rm -f "$stem".fastq.gz "$stem".zmw_metrics.json.gz
done

# 4. tidy intermediates.
mv -f "$OUTDIR"/sd_0001.maf.gz "$OUTDIR/pbsim.maf.gz" 2>/dev/null || true
rm -f "$OUTDIR"/sd_*.maf.gz "$OUTDIR"/sd_*.ref "$GENOME"
if [ "$KEEP_BAM" -eq 1 ]; then
    if [ "${#BAMS[@]}" -eq 1 ]; then
        mv -f "${BAMS[0]}" "$OUTDIR/subreads.bam"
    fi
else
    rm -f "$OUTDIR"/sd_*.bam
fi

# 5. summary line.
N=$(zcat "$OUTDIR/hifi.fastq.gz" | awk 'NR%4==2{n++; t+=length($0)} END{print n"\t"t}')
READS=$(echo "$N" | cut -f1)
TOTAL=$(echo "$N" | cut -f2)
MEAN=$(awk -v r="$READS" -v t="$TOTAL" 'BEGIN{ if(r>0) printf "%.0f",t/r; else print 0 }')
echo "[sim_reads] done: $OUTDIR/hifi.fastq.gz reads=$READS total_bp=$TOTAL mean_len=$MEAN"
