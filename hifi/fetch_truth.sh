#!/usr/bin/env bash
# fetch_truth.sh: truth sets and published call sets of the new tumor/normal pairs, into hifi/truth/<set>/.
#   colo829: somatic SV truth set, Valle-Inclan et al. 2022 (Zenodo 5019338), md5 from Zenodo
#   hcc1395: SEQC2 SV tables, Talsania et al. 2022, Genome Biology 23:255 (Additional files 3, 4, 5, 15; S4 is the
#            high-confidence consensus call set)
#   hg009:   NIST's Severus somatic SV calls for HG009-T p16 vs HG009-N-WT p4 (no benchmark exists yet)
#   prcgap:  the PRCGAP paper's result tarball (Zenodo 21231495): calls for NCI-H209, NCI-H2126, HCC1954; md5 from Zenodo
set -uo pipefail
T=/home/yuio/snape/hifi/truth; log() { echo "$(date '+%F %T') $*"; }
zen() {   # $1 record, $2 folder: every file of a Zenodo record, md5 checked
    mkdir -p "$T/$2"
    curl -s "https://zenodo.org/api/records/$1" | python3 -c "
import json,sys
for f in json.load(sys.stdin)['files']: print(f['key'], f['checksum'].split(':')[1], f['links']['self'])" |
    while read -r k md5 url; do
        [ -s "$T/$2/$k" ] && [ "$(md5sum "$T/$2/$k" | cut -d' ' -f1)" = "$md5" ] && { log "[$2] $k present"; continue; }
        # Zenodo throttles each connection to ~80 KB/s (2026-10-01): 16 connections of 1 MiB pieces gave 4.4 MB/s
        aria2c -x16 -s16 -k1M --file-allocation=none --console-log-level=warn --summary-interval=0 --max-tries=0 \
               --retry-wait=10 --auto-file-renaming=false -d "$T/$2" -o "$k" "$url" > /dev/null
        [ "$(md5sum "$T/$2/$k" | cut -d' ' -f1)" = "$md5" ] && log "[$2] $k md5 ok" || log "[FATAL] $2 $k md5 mismatch"
    done
}
zen 5019338 colo829
zen 21231495 prcgap
mkdir -p "$T/hcc1395"
for n in 3 4 5 15; do
    f=13059_2022_2816_MOESM${n}_ESM.xlsx
    [ -s "$T/hcc1395/$f" ] || curl -sSfL --retry 5 -o "$T/hcc1395/$f" "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-022-02816-6/MediaObjects/$f" \
        && log "[hcc1395] $f" || log "[FATAL] hcc1395 $f"
done
mkdir -p "$T/hg009"; G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG009/NIST/HG009-T_bulk/20250715p16/analysis/GRCh38/BCM_Revio_HG009T-p16_20260313/HG009-T_vs_HG009-N-WT-p4/structural_variants
for f in HG009T-p16.severus.somatic.vcf.gz HG009T-p16.severus.somatic.vcf.gz.tbi HG009T-p16.breakpoint_clusters.tsv; do
    curl -sSf --retry 5 -o "$T/hg009/$f" "$G/$f" && log "[hg009] $f" || log "[FATAL] hg009 $f"
done
cat > "$T/SOURCES.md" <<'EOF'
# Truth and call sets of the added pairs
- colo829/: somatic SV truth set, Valle-Inclan et al., Cell Genomics 2022; Zenodo 5019338 (truthset_somaticSVs_COLO829_hg38lifted.vcf is the GRCh38 set)
- hcc1395/: SEQC2 SV study, Talsania et al., Genome Biology 23:255 (2022), Additional files 3 (initial), 4 (high-confidence consensus), 5 (validation), 15 (BND summary)
- hg009/: NIST Severus somatic SV calls, HG009-T p16 vs HG009-N-WT p4 (GRCh38); not a benchmark
- prcgap/: PRCGAP results (Sakamoto et al., bioRxiv 10.64898/2026.05.28.728591), Zenodo 21231495; calls for NCI-H209, NCI-H2126, HCC1954 in personalized-reference coordinates
EOF
log "done"
