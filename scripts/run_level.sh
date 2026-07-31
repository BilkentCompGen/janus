#!/bin/bash
# Build an svp graph at a given LCP level and run the full trace analysis (identical metrics
# to the L5/L4 comparison). Reusable per-level driver.
#   usage: run_level.sh <level> <graphdir> <tracedir> [SVP_NOCACHE]
# Box rule: svp builds ALWAYS taskset -c 0-15 + --threads 16.
cd /home/yuio/svp
LEVEL=$1; GDIR=$2; TDIR=$3; NOCACHE=${4:-}
mkdir -p "$GDIR" "$TDIR"
LOG="outputs/svp/l${LEVEL}_run.log"
{
echo "[L$LEVEL] $(date '+%H:%M:%S') build start (cache=$([ -n "$NOCACHE" ] && echo off || echo on))"
env ${NOCACHE:+SVP_NOCACHE=1} taskset -c 0-15 /usr/bin/time -v ./svp graph \
    simulation/reads/normal/hifi.fastq simulation/reads/hard/hifi.fastq "$GDIR" \
    --level "$LEVEL" --threads 16
if [ ! -f "$GDIR/graph.bin" ]; then echo "[L$LEVEL] BUILD FAILED"; exit 1; fi
echo "[L$LEVEL] build done; graph.bin=$(du -h "$GDIR/graph.bin"|cut -f1)"
echo "[L$LEVEL] === trace run ==="
SVP_TRACE_LEVEL=$LEVEL TRACE_THREADS=16 ./tools/trace run "$GDIR/graph.bin" \
    simulation/ref.fasta simulation/truth/hard/out.fa simulation/truth/hard/out.vcf "$TDIR" 2>&1 \
  | grep -E 'LCP level|recall|nodes:|marking|compact'
echo "[L$LEVEL] === diag ==="
SVP_TRACE_LEVEL=$LEVEL ./tools/trace diag "$GDIR/graph.bin" simulation/ref.fasta \
    simulation/truth/hard/out.vcf "$TDIR/variants.tsv" "$TDIR" 2>&1 \
  | grep -E 'MISSED|FOUND|attribution|undercov|nonnovel|no_context|no_novel|low_cov'
echo "[L$LEVEL] === level_report ==="
python3 scripts/level_report.py "$GDIR/graph.bin" "$TDIR" "LEVEL $LEVEL"
echo "[L$LEVEL] DONE $(date '+%H:%M:%S')"
} >> "$LOG" 2>&1
