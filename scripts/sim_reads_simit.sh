#!/usr/bin/env bash
# sim_reads_simit.sh — realistic PacBio Revio HiFi reads via Sim-it (replaces pbsim3+ccs).
#
# Reads-only: the variants are already baked into the input genome
#   - tumor  = simulation/truth/wgs/out.fa  (insilicoSV: ref + somatic SNP/indel/SV)
#   - normal = simulation/ref.fasta         (the bare reference)
# so all SV counts are set to 0 here (Sim-it adds NO new SVs; out.vcf stays the truth).
#
# HiFi params from the current PacBio Revio standard (median read accuracy >=Q30 / 99.9%,
# read length ~15 kb, N50 >15 kb) + Sim-it's context-dependent HiFi error profile.
#
# Sim-it simulates a diploid: it writes HAP1, HAP2 and the merged HAP12. With a haploid
# input + 0 SVs the two haplotypes are independent samplings of the same sequence, so
# HAP12 == 2x per-haplotype depth. Measured: HAP12_cov = 1.863 * depth  =>  depth 16 ~= 30x.
#
# Usage: sim_reads_simit.sh -i <genome.fa> -o <out_dir> -n <project> [-d depth=16] [-t threads=16]
set -euo pipefail

SIMIT_DIR=/home/yuio/svp/lib/Sim-it
PROFILE="$SIMIT_DIR/error_profile_PB_Sequel_CCS_hifi.txt"
PERL5LIB_LOCAL="$HOME/perl5/lib/perl5"

DEPTH=16
THREADS=16
MEDIAN=15000
LENRANGE="1000-50000"
ACCURACY="99.9%"
PROJECT="reads"
INPUT=""
OUTDIR=""
while [ $# -gt 0 ]; do
  case "$1" in
    -i) INPUT="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    -n) PROJECT="$2"; shift 2;;
    -d) DEPTH="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    *) echo "[simit] unknown arg: $1" >&2; exit 2;;
  esac
done
[ -n "$INPUT" ] && [ -n "$OUTDIR" ] || { echo "[simit] need -i <genome.fa> and -o <out_dir>" >&2; exit 2; }
[ -f "$INPUT" ] || { echo "[simit] input not found: $INPUT" >&2; exit 2; }
[ -f "$PROFILE" ] || { echo "[simit] HiFi error profile missing: $PROFILE" >&2; exit 2; }
[ -d "$PERL5LIB_LOCAL" ] || { echo "[simit] Parallel::ForkManager not installed at $PERL5LIB_LOCAL" >&2; exit 2; }

mkdir -p "$OUTDIR"
WORK="$OUTDIR/simit_work"
mkdir -p "$WORK"
CFG="$WORK/config_${PROJECT}.txt"

cat > "$CFG" <<EOF
Project:
-----------------------------
Project name             = $PROJECT
Reference sequence       = $INPUT
Replace ambiguous nts(N) = no
Max threads              = $THREADS


Structural variation:
-----------------------------
VCF input                =
Foreign sequences        =

Deletions                = 0
Length (bp)              = 30-150000

Insertions               = 0
Length (bp)              = 30-100000

Tandem duplications      = 0
Length (bp)              = 50-10000
Copies                   = 1-20

Inversions               = 0
Length (bp)              = 150-1300000

Complex substitutions    = 0
Length (bp)              = 30-100000

Inverted duplications    = 0
Length (bp)              = 150-350000

Heterozygosity           = 0%


Long Read simulation:
-----------------------------
Sequencing depth         = $DEPTH
Median length            = $MEDIAN
Length range             = $LENRANGE
Accuracy                 = $ACCURACY
Error profile            = $PROFILE
EOF

EST=$(awk "BEGIN{printf \"%.0f\", 1.863*$DEPTH}")
echo "[simit] $PROJECT: input=$INPUT depth=$DEPTH (HAP12 ~= ${EST}x) threads=$THREADS acc=$ACCURACY median=$MEDIAN"
t0=$SECONDS
env PERL5LIB="$PERL5LIB_LOCAL" perl "$SIMIT_DIR/Sim-it1.4.5.pl" -c "$CFG" -o "$WORK"
echo "[simit] Sim-it generation done in $((SECONDS-t0))s"

HAP12="$WORK/Long_reads_${PROJECT}_HAP12.fasta"
[ -f "$HAP12" ] || { echo "[simit] expected output missing: $HAP12" >&2; exit 4; }

# coverage report (uses the input genome length)
GLEN=$(awk '!/^>/{b+=length($0)} END{print b}' "$INPUT")
RBP=$(awk '!/^>/{b+=length($0)} END{print b}' "$HAP12")
RN=$(grep -c '^>' "$HAP12" || true)
echo "[simit] reads=$RN  read_bases=$RBP  genome=$GLEN  coverage=$(awk "BEGIN{printf \"%.2f\",$RBP/$GLEN}")x"

# FASTA (single-line reads) -> FASTQ with uniform Q30 quality ('?'); svp ignores quality.
GZ="gzip"
command -v pigz >/dev/null 2>&1 && GZ="pigz -p $THREADS"
echo "[simit] converting HAP12 -> hifi.fastq.gz via $GZ"
perl -ne 'if(/^>/){s/^>/@/;print}else{chomp;print "$_\n+\n", "?" x length($_), "\n"}' "$HAP12" \
  | $GZ > "$OUTDIR/hifi.fastq.gz.tmp"
mv -f "$OUTDIR/hifi.fastq.gz.tmp" "$OUTDIR/hifi.fastq.gz"
echo "[simit] wrote $OUTDIR/hifi.fastq.gz ($(du -h "$OUTDIR/hifi.fastq.gz" | cut -f1))"

# free the ~2x-sized FASTA intermediates (keep config + log_*.txt)
rm -f "$WORK"/Long_reads_*HAP*.fasta
echo "[simit] $PROJECT done."
