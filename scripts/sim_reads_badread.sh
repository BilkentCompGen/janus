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
FORMAT=fastq          # fastq -> one hifi.fastq (default, unchanged); svb -> 2-bit shard.<i>.svb dir
while [ $# -gt 0 ]; do
  case "$1" in
    -i) INPUT="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    -d) DEPTH="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -s) SEED="$2"; shift 2;;
    -f) FORMAT="$2"; shift 2;;
    *) echo "[badread] unknown arg: $1" >&2; exit 2;;
  esac
done
[ "$FORMAT" = fastq ] || [ "$FORMAT" = svb ] || { echo "[badread] -f must be fastq|svb" >&2; exit 2; }
SVB_PACK="$ROOT/snape pack"
[ "$FORMAT" != svb ] || [ -x "$SVB_PACK" ] || { echo "[badread] svb_pack missing: $SVB_PACK (run 'make')" >&2; exit 2; }
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
    # svb: pack Badread's FASTQ straight to a 2-bit shard (no FASTQ on disk).
    # fastq: written uncompressed -- snape reads plain FASTQ only. sink chosen by FORMAT.
    if [ "$FORMAT" = svb ]; then SINK=("$SVB_PACK" "$OUTDIR/shard.$i.svb"); else SINK=(cat); fi
    "$BADREAD" simulate --reference "$WORK/shard.$i.fa" --quantity "${DEPTH}x" \
      --error_model "$MODEL" --qscore_model "$MODEL" \
      --identity "$IDENTITY" --length "$LENGTH" \
      --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 \
      --start_adapter 0,0 --end_adapter 0,0 --seed $((SEED+i)) \
      2> "$WORK/shard.$i.log" | "${SINK[@]}" > "$WORK/shard.$i.out" 2>> "$WORK/shard.$i.log"
    if [ "$FORMAT" != svb ]; then mv -f "$WORK/shard.$i.out" "$WORK/shard.$i.fastq"; else rm -f "$WORK/shard.$i.out"; fi
    rm -f "$WORK/shard.$i.fa"
  ) &
  pids+=("$!")
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ "$fail" -eq 0 ] || { echo "[badread] a shard failed — see $WORK/shard.*.log" >&2; exit 4; }
echo "[badread] simulation done in $((SECONDS-t0))s"

GLEN=$(awk '{b+=$2} END{print b}' "$INPUT.fai")
if [ "$FORMAT" = svb ]; then
  # shards ARE the output ($OUTDIR/shard.<i>.svb); no concat. snape reads the dir in parallel.
  rm -f "$WORK"/shard.*.regions
  read -r RN RBP < <(python3 - "$OUTDIR" <<'PY'
import sys, glob, struct
n=t=0
for f in sorted(glob.glob(sys.argv[1]+"/shard.*.svb")):
    h=open(f,'rb').read(32); assert h[:4]==b'SVB1', f
    nr,tb=struct.unpack_from('<QQ',h,8); n+=nr; t+=tb
print(n, t)
PY
)
  echo "[badread] reads=$RN read_bases=$RBP genome=$GLEN coverage=$(awk "BEGIN{printf \"%.2f\",$RBP/$GLEN}")x"
  echo "[badread] wrote $(ls "$OUTDIR"/shard.*.svb | wc -l) svb shards to $OUTDIR ($(du -sh "$OUTDIR" | cut -f1))"
  printf 'reads=%s read_bases=%s coverage=%.2fx shards=%s\n' "$RN" "$RBP" "$(awk "BEGIN{printf \"%.2f\",$RBP/$GLEN}")" "$(ls "$OUTDIR"/shard.*.svb | wc -l)" > "$OUTDIR/svb.done"
else
  # concatenate shard fastqs -- uncompressed, because snape reads plain FASTQ only
  cat "$WORK"/shard.*.fastq > "$OUTDIR/hifi.fastq.tmp"
  mv -f "$OUTDIR/hifi.fastq.tmp" "$OUTDIR/hifi.fastq"
  rm -f "$WORK"/shard.*.fastq "$WORK"/shard.*.regions
  read -r RN RBP < <(awk 'NR%4==2{n++; t+=length($0)} END{print n, t}' "$OUTDIR/hifi.fastq")
  echo "[badread] reads=$RN read_bases=$RBP genome=$GLEN coverage=$(awk "BEGIN{printf \"%.2f\",$RBP/$GLEN}")x"
  echo "[badread] wrote $OUTDIR/hifi.fastq ($(du -h "$OUTDIR/hifi.fastq" | cut -f1))"
fi
