#!/usr/bin/env bash
# COLO829 scoring, against the WHOLE Valle-Inclan truth set: 134 records = 66 mated breakend pairs
# (one junction each) plus two unmated <INS> records, scored on position and length instead.
#
# Minda is the instrument, at the settings the Severus paper used for this exact benchmark: 500 bp
# tolerance, 30 bp minimum SV length, and only calls at VAF >= 10% considered (for every tool). Those
# published numbers -- F1 0.79 PacBio, 0.74 ONT R9 -- are the reference point this run is measured
# against, so the protocol has to match them.
#
# truvari is NOT run. It cannot process this truth set in this environment: it aborts in pysam with
# "Can't write record: [Errno 22] tp-base.vcf" even when scoring the truth against itself, and even
# when cut down to two bracket-free records (python 3.14.6 / pysam 0.23.3 / truvari 5.4.0). It works
# on the simulation VCFs, so this is specific to these records, not to the callers.
#
# The VAF filter is worth watching rather than trusting: in the contig arm there are no
# reference-supporting reads by construction, so every call comes out near VAF 1.0 and the filter
# removes almost nothing. It is applied anyway, to keep the protocol identical to the paper's.
set -euo pipefail
cd /home/yuio/svp/colo829
SVP=/home/yuio/svp
LV=${LV:-5}
TRUTH=COLO829_somaticSV/truthset_somaticSVs_COLO829.vcf
RAW=sev_ctg_L$LV/all_SVs/severus_all.vcf
VCF=sev_ctg_L$LV/severus_vaf10.vcf
OUT=minda_out_L$LV
STAMP=$(date +%Y%m%d-%H%M%S)
TRASH=/home/yuio/trash/colo-score-$STAMP

python3 - "$RAW" "$VCF" <<'PY'
import re, sys
src, out = sys.argv[1], sys.argv[2]
kept = tot = 0
with open(out, 'w') as o:
    for line in open(src):
        if line[0] == '#':
            o.write(line); continue
        tot += 1
        f = line.rstrip('\n').split('\t')
        d = dict(zip(f[8].split(':'), f[9].split(':')))
        try: vaf = float(d.get('VAF', '0'))
        except ValueError: vaf = 0.0
        m = re.search(r'SVLEN=(-?\d+)', f[7])
        svlen = abs(int(m.group(1))) if m else 0
        if vaf >= 0.10 and (svlen >= 30 or 'SVTYPE=BND' in f[7]):
            o.write(line); kept += 1
print(f"VAF>=10% + >=30bp: {tot} -> {kept} cagri")
PY

echo; echo "############ MINDA (500 bp, min 30 bp -- makalenin protokolu) ############"
[ -d "$OUT" ] && { mkdir -p "$TRASH"; mv "$OUT" "$TRASH/"; }
mkdir -p "$OUT"
"$SVP/.venv/bin/minda" truthset --base "$TRUTH" --vcfs "$VCF" \
    --out_dir "$OUT" --tolerance 500 --min_size 30 2>&1 | tail -4
echo; cat "$OUT/results/overall.tsv"

echo; echo "############ junction seviyesi (kendi skorlayicim, +/-500 bp) ############"
python3 "$SVP/scripts/colo_bp_eval.py" "$TRUTH" "$VCF" 500

echo; echo "############ referans: Severus makalesi, ayni truth set + ayni Minda ############"
echo "  PacBio   F1 0.79"
echo "  ONT R9   F1 0.74   (tumor/normal modu, 100x/100x)"
