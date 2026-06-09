#!/usr/bin/env bash
# colo829_benchmark.sh — truvari benchmark of svp's calls.vcf against the COLO829
# somatic-SV truth (Valle-Inclan 2022, GRCh38-lifted). SV-only truth: no SNP/indel tiers.
#
# Naming already matches: both ref.fasta and the truth use un-prefixed contigs (1..22,X).
# svp can only emit DEL/INS, so the headline is the >=50bp DEL+INS recall/precision;
# DUP/INV/BND in the truth are reported as structurally-out-of-scope FNs for context.
#
# Usage: scripts/colo829_benchmark.sh <svp_out_dir>   (expects <dir>/calls.vcf)
set -euo pipefail
source /home/yuio/svp/.venv/bin/activate

DIR="${1:?usage: $0 <svp_out_dir>}"
CALLS="$DIR/calls.vcf"
TRUTH_RAW=/home/yuio/svp/colo829/truth/truthset_somaticSVs_COLO829_hg38lifted.vcf
REF=/home/yuio/svp/simulation/ref.fasta
[ -f "$CALLS" ] || { echo "FATAL: $CALLS missing"; exit 1; }
[ -f "$TRUTH_RAW" ] || { echo "FATAL: truth $TRUTH_RAW missing"; exit 1; }

prep_vcf() {
    local src="$1" dst="$2"
    if [ -f "$dst" ] && [ -f "$dst.tbi" ] && [ "$dst" -nt "$src" ]; then echo "[skip] $dst prepped"; return; fi
    bcftools sort "$src" -O z -o "$dst"
    tabix -f -p vcf "$dst"
}

prep_vcf "$TRUTH_RAW" "$DIR/colo829_truth.vcf.gz"
prep_vcf "$CALLS"     "$DIR/calls.vcf.gz"

run_bench() {
    local label="$1" out="$DIR/truvari_$1"; shift
    rm -rf "$out"
    truvari bench -b "$DIR/colo829_truth.vcf.gz" -c "$DIR/calls.vcf.gz" -f "$REF" -o "$out" "$@"
    echo "=== $label ==="
    python3 -c "
import json
s=json.load(open('$out/summary.json'))
for k in ('TP-base','TP-comp','FP','FN','precision','recall','f1','base cnt','comp cnt'):
    if k in s: v=s[k]; print(f'  {k:>12}: '+(f'{v:.4f}' if isinstance(v,float) else str(v)))
"
}

# Headline: >=50bp, sequence-agnostic (real reads -> allele bytes won't match truth exactly).
run_bench "geq50" --sizemin 50 --sizefilt 50 --pctsize 0.7 --pctseq 0 --refdist 500 --sizemax 1000000000

echo
echo "=== per-SVTYPE breakdown (geq50) ==="
D="$DIR/truvari_geq50"
for v in tp-base fn fp; do
    f="$D/$v.vcf.gz"; [ -f "$f" ] || f="$D/$v.vcf"
    printf "%-9s " "$v:"
    zcat -f "$f" 2>/dev/null | grep -v '^#' | grep -oE 'SVTYPE=[A-Za-z]+' | sort | uniq -c | tr '\n' ' '; echo
done
echo
echo "=== truth SVTYPE totals (for reference; svp reaches DEL+INS only) ==="
grep -v '^#' "$TRUTH_RAW" | grep -oE 'SVTYPE=[A-Za-z]+' | sort | uniq -c
