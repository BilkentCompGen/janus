#!/usr/bin/env bash
# sim_reads_badread.sh — realistic PacBio HiFi reads via Badread (replaces Sim-it).
#
# Reads-only: the variants are already baked into the input genome (rule 1) —
#   tumor  = simulation/truth/<id>/out.fa  (insilicoSV: ref + somatic SNP/indel/SV)
#   normal = simulation/ref.fasta          (the bare reference)
# Badread adds NO variants; it only samples reads and injects a realistic HiFi
# error/qscore profile. The truth VCF (out.vcf) stays the ground truth.
#
# Profile: Q30 mean / stdev 3 (--identity 30,3, normal qscore distribution) over the
# pacbio2021 model — Badread's *recent* PacBio Sequel IIe HiFi model. (pacbio2016 is
# the old RS II model and is deliberately NOT used.) HiFi reads are clean CCS
# consensus, so adapters/junk/random/chimeras/glitches are all off.
#
# Badread is single-threaded and runs at a fixed bp/s, so for a whole-genome reference
# we split the contigs into ~equal chunks (a whole chromosome cannot be split by Badread
# itself and would otherwise monopolize one core for hours) and greedy bin-pack the
# chunks into THREADS balanced bins, running one Badread per bin in parallel at the
# target depth, then concatenate. Per-region depth == whole-genome depth.
#
# Usage: sim_reads_badread.sh -i <genome.fa> -o <out_dir> [-d depth=30] [-t threads=16] [-s seed=1]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BADREAD="$ROOT/lib/Badread/venv/bin/badread"
SAMTOOLS="${SAMTOOLS:-samtools}"

DEPTH=30
THREADS=16
SEED=1
IDENTITY="30,3"
MODEL=pacbio2021
LENGTH="15000,4000"
INPUT=""
OUTDIR=""
while [ $# -gt 0 ]; do
  case "$1" in
    -i) INPUT="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    -d) DEPTH="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -s) SEED="$2"; shift 2;;
    *) echo "[badread] unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$INPUT" ] && [ -n "$OUTDIR" ] || { echo "[badread] need -i <genome.fa> and -o <out_dir>" >&2; exit 2; }
[ -f "$INPUT" ] || { echo "[badread] input not found: $INPUT" >&2; exit 2; }
[ -x "$BADREAD" ] || { echo "[badread] badread missing: $BADREAD (run 'make install-deps')" >&2; exit 2; }
[ "$THREADS" -ge 1 ] 2>/dev/null || THREADS=$(nproc)   # caller may pass 0 ("auto"); Badread is 1-threaded so we shard

mkdir -p "$OUTDIR"
WORK="$OUTDIR/badread_work"
mkdir -p "$WORK"
[ -f "$INPUT.fai" ] || "$SAMTOOLS" faidx "$INPUT"

# Split contigs into ~equal chunks and greedy bin-pack them into THREADS bins balanced
# by bp -> shard.<i>.regions (newline-separated samtools region strings). Chunking is
# what keeps every core busy: without it the largest chromosome's whole-contig shard
# would run single-threaded for hours while the other cores sit idle.
rm -f "$WORK"/shard.*.regions
python3 - "$INPUT.fai" "$THREADS" "$WORK" <<'PY'
import sys, math
fai, T, work = sys.argv[1], int(sys.argv[2]), sys.argv[3]
contigs = []; total = 0
for line in open(fai):
    f = line.split('\t'); contigs.append((f[0], int(f[1]))); total += int(f[1])
chunk = max(1, math.ceil(total / (T * 8)))   # ~8 chunks per thread -> tight packing
units = []                                    # (region_string, bp)
for name, ln in contigs:
    s = 1
    while s <= ln:
        e = min(s + chunk - 1, ln)
        units.append((f"{name}:{s}-{e}", e - s + 1)); s = e + 1
units.sort(key=lambda u: -u[1])
load = [0]*T; bins = [[] for _ in range(T)]
for reg, bp in units:                         # largest chunk -> currently-lightest bin
    i = min(range(T), key=lambda b: load[b]); bins[i].append(reg); load[i] += bp
for i in range(T):
    open(f"{work}/shard.{i}.regions", "w").write("".join(r + "\n" for r in bins[i]))
print(f"[badread] {len(contigs)} contigs -> {len(units)} chunks -> {T} bins; "
      f"bin bp max={max(load):,} min={min(load):,}")
PY

echo "[badread] $OUTDIR: depth=${DEPTH}x identity=$IDENTITY model=$MODEL length=$LENGTH threads=$THREADS seed=$SEED"
t0=$SECONDS
pids=()
for i in $(seq 0 $((THREADS-1))); do
  REGIONS="$WORK/shard.$i.regions"
  [ -s "$REGIONS" ] || continue
  (
    "$SAMTOOLS" faidx "$INPUT" $(cat "$REGIONS") > "$WORK/shard.$i.fa"
    "$BADREAD" simulate --reference "$WORK/shard.$i.fa" --quantity "${DEPTH}x" \
      --error_model "$MODEL" --qscore_model "$MODEL" \
      --identity "$IDENTITY" --length "$LENGTH" \
      --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 \
      --start_adapter 0,0 --end_adapter 0,0 --seed $((SEED+i)) \
      2> "$WORK/shard.$i.log" | gzip > "$WORK/shard.$i.fastq.gz"
    rm -f "$WORK/shard.$i.fa"
  ) &
  pids+=("$!")
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ "$fail" -eq 0 ] || { echo "[badread] a shard failed — see $WORK/shard.*.log" >&2; exit 4; }
echo "[badread] simulation done in $((SECONDS-t0))s"

# concatenate shard fastqs (gzip members concat -> valid gz stream)
cat "$WORK"/shard.*.fastq.gz > "$OUTDIR/hifi.fastq.gz.tmp"
mv -f "$OUTDIR/hifi.fastq.gz.tmp" "$OUTDIR/hifi.fastq.gz"
rm -f "$WORK"/shard.*.fastq.gz "$WORK"/shard.*.regions

GLEN=$(awk '{b+=$2} END{print b}' "$INPUT.fai")
read -r RN RBP < <(zcat "$OUTDIR/hifi.fastq.gz" | awk 'NR%4==2{n++; t+=length($0)} END{print n, t}')
echo "[badread] reads=$RN read_bases=$RBP genome=$GLEN coverage=$(awk "BEGIN{printf \"%.2f\",$RBP/$GLEN}")x"
echo "[badread] wrote $OUTDIR/hifi.fastq.gz ($(du -h "$OUTDIR/hifi.fastq.gz" | cut -f1))"
