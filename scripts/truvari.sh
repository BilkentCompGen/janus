#!/usr/bin/env bash
set -euo pipefail

source .venv/bin/activate

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <base.vcf> <comp.vcf> <reference.fa> [output_dir]"
    exit 1
fi

BASE_VCF="$1"
COMP_VCF="$2"
REF="$3"
OUT_DIR="${4:-outputs/truvari}"

mkdir -p "$OUT_DIR"

prep_vcf() {
    local vcf="$1"
    local out="$2"
    bcftools sort "$vcf" | bgzip -c > "$out"
    tabix -p vcf "$out"
}

echo "Preparing base VCF..."
prep_vcf "$BASE_VCF" "$OUT_DIR/base.vcf.gz"

echo "Preparing comp VCF..."
prep_vcf "$COMP_VCF" "$OUT_DIR/comp.vcf.gz"

echo "Running truvari bench..."
rm -rf "$OUT_DIR/results"
truvari bench \
    -b "$OUT_DIR/base.vcf.gz" \
    -c "$OUT_DIR/comp.vcf.gz" \
    -f "$REF" \
    --pctseq 0 \
    --sizemax 1000000000 \
    -o "$OUT_DIR/results"

echo ""
echo "=== Truvari Results ==="
python3 -c "
import json, sys, gzip, re

with open('$OUT_DIR/results/summary.json') as f:
    s = json.load(f)
for k in ['TP-base', 'TP-comp', 'FP', 'FN', 'precision', 'recall', 'f1']:
    if k in s:
        v = s[k]
        print(f'  {k:>20}: {v if v is None else (f\"{v:.4f}\" if isinstance(v, float) else v)}')

def parse_sv_entries(vcf_gz):
    entries = []
    with gzip.open(vcf_gz, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            chrom, pos, info = fields[0], int(fields[1]), fields[7]
            svlen_m = re.search(r'SVLEN=(\d+)', info)
            svtype_m = re.search(r'SVTYPE=(\w+)', info)
            svlen = int(svlen_m.group(1)) if svlen_m else 0
            svtype = svtype_m.group(1) if svtype_m else '?'
            entries.append((chrom, pos, svtype, svlen))
    return entries

def print_entries(label, entries):
    if not entries:
        print(f'\n  {label}: (none)')
        return
    by_type = {}
    for c, p, t, l in entries:
        by_type.setdefault(t, []).append((c, p, l))
    print(f'\n  {label} ({len(entries)} total):')
    for t in sorted(by_type):
        items = sorted(by_type[t], key=lambda x: x[2])
        sizes = [x[2] for x in items]
        print(f'    {t}: {len(items)}  (sizes: {\" \".join(str(s) for s in sizes)})')
    # Show the big ones in detail
    big = [(c, p, t, l) for c, p, t, l in entries if l >= 1000]
    if big:
        print(f'    Large (>=1kb):')
        for c, p, t, l in sorted(big, key=lambda x: -x[3]):
            print(f'      {c}:{p}  {t}  {l}bp')

fn = parse_sv_entries('$OUT_DIR/results/fn.vcf.gz')
fp = parse_sv_entries('$OUT_DIR/results/fp.vcf.gz')
print_entries('FN (missed)', fn)
print_entries('FP (false calls)', fp)

# Check if any FP might be merging multiple nearby FNs
print('\n  Nearby FP/FN pairs (within 500bp):')
found = False
for fc, fp_pos, ft, fl in fp:
    for nc, fn_pos, nt, nl in fn:
        if fc == nc and abs(fp_pos - fn_pos) < 500:
            found = True
            print(f'    FP {fc}:{fp_pos} {ft} {fl}bp  <->  FN {nc}:{fn_pos} {nt} {nl}bp  (dist={abs(fp_pos-fn_pos)})')
if not found:
    print('    (none)')
"
echo ""
echo "Full results: $OUT_DIR/results/summary.json"
