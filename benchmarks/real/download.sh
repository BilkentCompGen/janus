#!/usr/bin/env bash
# =============================================================================
# download.sh — Download real-world benchmark datasets for SVP
#
# Usage:
#   ./download.sh <dataset> [dataset2 ...]
#   ./download.sh all          # download everything
#   ./download.sh --list       # show available datasets
#
# Datasets (ordered by size / practicality):
#   grch38          GRCh38 reference genome (~850MB compressed)
#   colo829-truth   COLO829 somatic SV truth set VCF from Valle-Inclan et al.
#   hg002-truth     GIAB HG002 Tier 1 germline SV truth set
#   hg002-asm       HG002 T2T diploid assembly (~2GB compressed)
#   colo829-demo    colorSV demo FASTAs — small COLO829 SV-region subset
#   colo829-hifi    COLO829/COLO829BL PacBio Revio HiFi BAMs (~120GB total)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SAMPLES_DIR="$SCRIPT_DIR/samples"
TRUTH_DIR="$SCRIPT_DIR/ground_truth"

mkdir -p "$SAMPLES_DIR" "$TRUTH_DIR"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
info()  { echo ":: $*"; }
warn()  { echo "!! $*" >&2; }
die()   { echo "ERROR: $*" >&2; exit 1; }

download() {
    local url="$1" dest="$2"
    if [[ -f "$dest" ]]; then
        info "Already exists: $dest — skipping"
        return 0
    fi
    info "Downloading $(basename "$dest") ..."
    info "  URL: $url"
    local dir
    dir="$(dirname "$dest")"
    mkdir -p "$dir"
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "$dest" "$url"
    elif command -v curl &>/dev/null; then
        curl -L --progress-bar -o "$dest" "$url"
    else
        die "Neither wget nor curl found"
    fi
    info "  Saved: $dest ($(du -h "$dest" | cut -f1))"
}

require_cmd() {
    command -v "$1" &>/dev/null || die "Required command not found: $1"
}

# ---------------------------------------------------------------------------
# Dataset: GRCh38 reference genome
# ---------------------------------------------------------------------------
download_grch38() {
    info "=== GRCh38 Reference Genome ==="
    info "The 'no-alt analysis set' — standard for variant calling (~850MB gz)"

    local url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    local dest="$SAMPLES_DIR/GRCh38.fa.gz"

    download "$url" "$dest"

    # Decompress if needed
    if [[ -f "$dest" && ! -f "${dest%.gz}" ]]; then
        info "Decompressing GRCh38 reference..."
        gunzip -k "$dest"
    fi
    info "GRCh38 reference ready: ${dest%.gz}"
}

# ---------------------------------------------------------------------------
# Dataset: COLO829 somatic SV truth set (Valle-Inclan et al. 2022)
# ---------------------------------------------------------------------------
download_colo829_truth() {
    info "=== COLO829 Somatic SV Truth Set ==="
    info "68 validated somatic SVs from Valle-Inclan et al. (Cell Genomics 2022)"
    info "Source: Zenodo 4716169 — GRCh38 coordinates"

    local dest_dir="$TRUTH_DIR/colo829"
    mkdir -p "$dest_dir"

    # The Zenodo record has a zip with all files. Download individual files.
    # GRCh38-lifted truth set VCF (68 somatic SVs)
    local zenodo_api="https://zenodo.org/api/records/4716169/files"

    download "$zenodo_api/truthset_somaticSVs_COLO829_hg38lifted.vcf/content" \
             "$dest_dir/truthset_somaticSVs_COLO829_hg38lifted.vcf"

    # Original GRCh37 truth set
    download "$zenodo_api/truthset_somaticSVs_COLO829.vcf/content" \
             "$dest_dir/truthset_somaticSVs_COLO829.vcf"

    info "COLO829 truth set ready in: $dest_dir/"
    info "  Reference: Valle-Inclan et al. 2022, doi:10.1016/j.xgen.2022.100167"
}

# ---------------------------------------------------------------------------
# Dataset: HG002 GIAB germline SV truth set
# ---------------------------------------------------------------------------
download_hg002_truth() {
    info "=== HG002 GIAB Tier 1 SV Truth Set ==="
    info "12,745 sequence-resolved INS+DEL calls >=50bp (GRCh37 coordinates)"

    local dest_dir="$TRUTH_DIR/hg002"
    mkdir -p "$dest_dir"

    local base="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6"

    download "$base/HG002_SVs_Tier1_v0.6.vcf.gz" \
             "$dest_dir/HG002_SVs_Tier1_v0.6.vcf.gz"
    download "$base/HG002_SVs_Tier1_v0.6.vcf.gz.tbi" \
             "$dest_dir/HG002_SVs_Tier1_v0.6.vcf.gz.tbi"
    download "$base/HG002_SVs_Tier1_v0.6.bed" \
             "$dest_dir/HG002_SVs_Tier1_v0.6.bed"
    download "$base/README_SV_v0.6.txt" \
             "$dest_dir/README_SV_v0.6.txt"

    info "HG002 truth set ready in: $dest_dir/"
    info "  NOTE: v0.6 uses GRCh37 coordinates. See README for liftover info."
}

# ---------------------------------------------------------------------------
# Dataset: HG002 T2T diploid assembly
# ---------------------------------------------------------------------------
download_hg002_asm() {
    info "=== HG002 T2T Diploid Assembly ==="
    info "hg002v1.1 — Verkko + PacBio HiFi + ONT R10 ultralong (~2GB gz)"
    info "This is a complete diploid FASTA — directly usable as SVP input."

    download "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz" \
             "$SAMPLES_DIR/hg002v1.1.fasta.gz"

    # Decompress
    if [[ -f "$SAMPLES_DIR/hg002v1.1.fasta.gz" && ! -f "$SAMPLES_DIR/hg002v1.1.fasta" ]]; then
        info "Decompressing HG002 assembly..."
        gunzip -k "$SAMPLES_DIR/hg002v1.1.fasta.gz"
    fi

    info "HG002 assembly ready: $SAMPLES_DIR/hg002v1.1.fasta"
}

# ---------------------------------------------------------------------------
# Dataset: COLO829 PacBio Revio HiFi BAMs (large!)
# ---------------------------------------------------------------------------
download_colo829_hifi() {
    info "=== COLO829 PacBio Revio HiFi BAMs ==="
    info "60X HiFi reads for tumor (COLO829) and normal (COLO829BL)"
    info "WARNING: This download is ~120GB total. Ensure sufficient disk space."
    info ""
    info "These BAMs can be:"
    info "  1. Assembled into FASTA with hifiasm for SVP input"
    info "  2. Used directly by rival callers (Delly, Manta, etc.)"
    info ""

    local df_avail
    df_avail=$(df -BG --output=avail "$SAMPLES_DIR" 2>/dev/null | tail -1 | tr -d ' G' || echo "unknown")
    if [[ "$df_avail" != "unknown" && "$df_avail" -lt 150 ]]; then
        warn "Only ${df_avail}GB available. HiFi BAMs need ~120GB + assembly space."
        warn "Continue? (Ctrl+C to abort, Enter to continue)"
        read -r
    fi

    local hifi_dir="$SAMPLES_DIR/colo829_hifi"
    mkdir -p "$hifi_dir"

    local pacbio_base="https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829"

    # Tumor COLO829 HiFi BAMs
    info "--- Tumor (COLO829) ---"
    for bam in \
        "m84039_230312_025934_s1.hifi_reads.default.bam" \
        "m84039_230328_000836_s3.hifi_reads.default.bam"; do
        download "$pacbio_base/$bam" "$hifi_dir/$bam"
    done

    # Normal COLO829BL HiFi BAMs
    info "--- Normal (COLO829BL) ---"
    for bam in \
        "m84039_230318_031003_s4.hifi_reads.default.bam" \
        "m84039_230324_015837_s2.hifi_reads.default.bam"; do
        download "$pacbio_base/$bam" "$hifi_dir/$bam"
    done

    info "HiFi BAMs downloaded to: $hifi_dir/"
    info "Next step: run ./prepare.sh assemble-colo829 to build FASTA assemblies."
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
usage() {
    cat <<'USAGE'
Usage: ./download.sh <dataset> [dataset2 ...]

Available datasets (ordered by download size):
  colo829-truth   COLO829 somatic SV truth set VCF (~1MB)
  hg002-truth     GIAB HG002 Tier 1 germline SV truth set (~5MB)
  grch38          GRCh38 reference genome (~850MB compressed, ~3GB uncompressed)
  hg002-asm       HG002 T2T diploid assembly (~2GB compressed, ~6GB uncompressed)
  colo829-hifi    COLO829 PacBio Revio HiFi BAMs (~120GB — needs assembly)

Shorthand:
  truth           Download all truth sets (colo829-truth + hg002-truth)
  quick           GRCh38 + HG002 assembly + truth sets (~3GB) — full benchmark ready
  all             Download everything

Examples:
  ./download.sh quick                    # fast start: truth + demo + ref
  ./download.sh grch38 colo829-truth     # just what you need
  ./download.sh colo829-hifi             # full HiFi for assembly
USAGE
}

if [[ $# -eq 0 || "$1" == "--help" || "$1" == "-h" ]]; then
    usage
    exit 0
fi

if [[ "$1" == "--list" ]]; then
    usage
    exit 0
fi

# Expand shorthand
datasets=()
for arg in "$@"; do
    case "$arg" in
        truth)
            datasets+=(colo829-truth hg002-truth)
            ;;
        quick)
            datasets+=(grch38 hg002-truth hg002-asm colo829-truth)
            ;;
        all)
            datasets+=(grch38 colo829-truth hg002-truth hg002-asm colo829-hifi)
            ;;
        *)
            datasets+=("$arg")
            ;;
    esac
done

# Deduplicate
declare -A seen
unique_datasets=()
for ds in "${datasets[@]}"; do
    if [[ -z "${seen[$ds]:-}" ]]; then
        seen[$ds]=1
        unique_datasets+=("$ds")
    fi
done

info "Downloading ${#unique_datasets[@]} dataset(s): ${unique_datasets[*]}"
echo ""

for ds in "${unique_datasets[@]}"; do
    case "$ds" in
        grch38)          download_grch38 ;;
        colo829-truth)   download_colo829_truth ;;
        hg002-truth)     download_hg002_truth ;;
        hg002-asm)       download_hg002_asm ;;
        colo829-hifi)    download_colo829_hifi ;;
        *)               die "Unknown dataset: $ds (run with --list to see options)" ;;
    esac
    echo ""
done

info "All downloads complete."
info "Next steps:"
info "  ./prepare.sh index-ref         # index GRCh38 reference"
info "  ./prepare.sh assemble-colo829  # assemble HiFi BAMs (if downloaded)"
info "  ./run_svp.sh                   # run SVP benchmark"
