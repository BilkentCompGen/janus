#!/usr/bin/env bash
# Three-tier truvari benchmark on svp's calls.vcf: >=50bp SVs, 1-49bp indels, SNPs.
# Usage: scripts/svp_benchmark.sh <dir>   (expects <dir>/calls.vcf)
set -euo pipefail
source /home/yuio/svp/.venv/bin/activate

DIR="${1:?usage: $0 <dir>}"
CALLS="$DIR/calls.vcf"
TRUTH_RAW="/home/yuio/svp/simulation/truth/wgs/out.vcf"
REF="/home/yuio/svp/simulation/ref.fasta"
TRUTH_DIR="$DIR"

[ -f "$CALLS" ] || { echo "FATAL: $CALLS missing"; exit 1; }

prep_vcf() {
    local src="$1" dst="$2"
    if [ -f "$dst" ] && [ -f "$dst.tbi" ] && [ "$dst" -nt "$src" ]; then
        echo "[skip] $dst already prepped"
        return
    fi
    bcftools sort "$src" -O z -o "$dst"
    tabix -f -p vcf "$dst"
}

prep_vcf "$TRUTH_RAW" "$TRUTH_DIR/truth.vcf.gz"
prep_vcf "$CALLS"     "$DIR/calls.vcf.gz"

run_bench() {
    local label="$1" out="$DIR/truvari_$1"
    shift
    rm -rf "$out"
    truvari bench \
        -b "$TRUTH_DIR/truth.vcf.gz" \
        -c "$DIR/calls.vcf.gz" \
        -f "$REF" \
        -o "$out" \
        "$@"
    echo "=== $label ==="
    python3 -c "
import json
with open('$out/summary.json') as f: s = json.load(f)
for k in ('TP-base','TP-comp','FP','FN','precision','recall','f1','base cnt','comp cnt'):
    if k in s:
        v = s[k]
        if isinstance(v, float): print(f'  {k:>12}: {v:.4f}')
        else: print(f'  {k:>12}: {v}')
"
}

run_bench "geq50" --sizemin 50 --sizefilt 50 --pctsize 0.7 --pctseq 0 --refdist 500 --sizemax 1000000000
run_bench "lt50"  --sizemin 1  --sizefilt 1  --pctsize 0.7 --pctseq 0 --refdist 50  --sizemax 49
run_bench "snp"   --sizemin 0  --sizefilt 0  --refdist 0   --sizemax 1
