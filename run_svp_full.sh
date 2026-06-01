#!/bin/bash
# Full svp pipeline (lcp + graph + bubbles) on WG normal+wgs, 16 threads.
# Records /usr/bin/time -v for the full pipeline.
set -e
cd /home/yuio/svp
OUT=outputs/svp/wgs
LOG=outputs/svp/full_$(date +%Y%m%d_%H%M%S).log
exec >>"$LOG" 2>&1
echo "=========================================="
echo "svp full pipeline (lcp + graph + bubbles)"
echo "  normal: simulation/reads/normal/hifi.fastq.gz"
echo "  tumor : simulation/reads/wgs/hifi.fastq.gz"
echo "  threads: 16  level: 5  flank: 10  min-mult: 3"
echo "  outdir: $OUT"
echo "  start: $(date '+%F %T')"
echo "=========================================="
mkdir -p $OUT
/usr/bin/time -v ./svp run \
  simulation/reads/normal/hifi.fastq.gz \
  simulation/reads/wgs/hifi.fastq.gz \
  $OUT \
  --level 5 --flank 10 --min-mult 3 --threads 16
echo "=========================================="
echo "SVP-COMPLETE $(date '+%F %T')"
ls -la $OUT/
