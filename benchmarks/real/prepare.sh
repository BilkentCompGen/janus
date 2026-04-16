#!/usr/bin/env bash
# =============================================================================
# prepare.sh — Prepare downloaded real-world data for SVP benchmarking
#
# Usage:
#   ./prepare.sh <step> [step2 ...]
#   ./prepare.sh all
#
# Steps:
#   index-ref           Index GRCh38 reference (samtools faidx)
#   assemble-colo829    Assemble COLO829 tumor+normal from HiFi BAMs (hifiasm)
#   extract-chr         Extract per-chromosome FASTAs for SVP
#   prep-rival-bams     Sort/index BAMs for rival callers
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SVP_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SAMPLES_DIR="$SCRIPT_DIR/samples"
TRUTH_DIR="$SCRIPT_DIR/ground_truth"

THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"

info()  { echo ":: $*"; }
warn()  { echo "!! $*" >&2; }
die()   { echo "ERROR: $*" >&2; exit 1; }

require_cmd() {
    command -v "$1" &>/dev/null || die "Required command not found: $1. Install it first."
}

# ---------------------------------------------------------------------------
# Step: Index the GRCh38 reference
# ---------------------------------------------------------------------------
index_ref() {
    info "=== Indexing GRCh38 Reference ==="

    local ref="$SAMPLES_DIR/GRCh38.fa"
    [[ -f "$ref" ]] || ref="$SAMPLES_DIR/GRCh38.fa.gz"
    [[ -f "$ref" ]] || die "GRCh38 not found. Run: ./download.sh grch38"

    # Decompress if still gzipped
    if [[ "$ref" == *.gz ]]; then
        if [[ ! -f "${ref%.gz}" ]]; then
            info "Decompressing GRCh38..."
            gunzip -k "$ref"
        fi
        ref="${ref%.gz}"
    fi

    require_cmd samtools

    if [[ ! -f "${ref}.fai" ]]; then
        info "Creating FASTA index..."
        samtools faidx "$ref"
    else
        info "Index already exists: ${ref}.fai"
    fi

    info "GRCh38 reference indexed: $ref"
}

# ---------------------------------------------------------------------------
# Step: Assemble COLO829 from HiFi BAMs using hifiasm
# ---------------------------------------------------------------------------
assemble_colo829() {
    info "=== Assembling COLO829 from HiFi BAMs ==="
    info "Using hifiasm for haplotype-resolved assembly"
    info "This will take several hours and needs ~64GB RAM."

    require_cmd hifiasm
    require_cmd samtools

    local hifi_dir="$SAMPLES_DIR/colo829_hifi"
    [[ -d "$hifi_dir" ]] || die "HiFi BAMs not found. Run: ./download.sh colo829-hifi"

    local asm_dir="$SAMPLES_DIR/colo829_assembly"
    mkdir -p "$asm_dir"

    # --- Convert BAMs to FASTQ (hifiasm prefers FASTQ input) ---
    info "Converting HiFi BAMs to FASTQ..."

    # Tumor
    local tumor_fq="$asm_dir/colo829_tumor.fastq.gz"
    if [[ ! -f "$tumor_fq" ]]; then
        info "  Extracting tumor reads..."
        samtools fastq -@ "$THREADS" \
            "$hifi_dir"/m84039_230312_025934_s1.hifi_reads.default.bam \
            "$hifi_dir"/m84039_230328_000836_s3.hifi_reads.default.bam \
            | gzip > "$tumor_fq"
    fi

    # Normal
    local normal_fq="$asm_dir/colo829bl_normal.fastq.gz"
    if [[ ! -f "$normal_fq" ]]; then
        info "  Extracting normal reads..."
        samtools fastq -@ "$THREADS" \
            "$hifi_dir"/m84039_230318_031003_s4.hifi_reads.default.bam \
            "$hifi_dir"/m84039_230324_015837_s2.hifi_reads.default.bam \
            | gzip > "$normal_fq"
    fi

    # --- Assemble tumor ---
    info "Assembling COLO829 tumor..."
    local tumor_prefix="$asm_dir/colo829_tumor"
    if [[ ! -f "${tumor_prefix}.bp.p_ctg.gfa" ]]; then
        hifiasm -o "$tumor_prefix" -t "$THREADS" "$tumor_fq"
    else
        info "  Tumor assembly already exists, skipping."
    fi

    # --- Assemble normal ---
    info "Assembling COLO829BL normal..."
    local normal_prefix="$asm_dir/colo829bl_normal"
    if [[ ! -f "${normal_prefix}.bp.p_ctg.gfa" ]]; then
        hifiasm -o "$normal_prefix" -t "$THREADS" "$normal_fq"
    else
        info "  Normal assembly already exists, skipping."
    fi

    # --- Convert GFA to FASTA ---
    info "Converting assembly GFA to FASTA..."

    for prefix in "$tumor_prefix" "$normal_prefix"; do
        local gfa="${prefix}.bp.p_ctg.gfa"
        local fa="${prefix}.fa"
        if [[ -f "$gfa" && ! -f "$fa" ]]; then
            awk '/^S/{print ">"$2; print $3}' "$gfa" > "$fa"
            info "  Created: $fa ($(du -h "$fa" | cut -f1))"
        fi
    done

    # Create symlinks with convenient names
    ln -sf "colo829_assembly/colo829_tumor.fa" "$SAMPLES_DIR/colo829_tumor.fa" 2>/dev/null || true
    ln -sf "colo829_assembly/colo829bl_normal.fa" "$SAMPLES_DIR/colo829_normal.fa" 2>/dev/null || true

    info "Assembly complete."
    info "  Tumor:  $SAMPLES_DIR/colo829_tumor.fa"
    info "  Normal: $SAMPLES_DIR/colo829_normal.fa"
}

# ---------------------------------------------------------------------------
# Step: Prepare HG002 — rename chromosomes and extract paternal haplotype
# SVP matches chromosomes by bare number ("19"), so we rename:
#   GRCh38: "chr19 AC:..." -> "19"
#   HG002:  "chr19_PATERNAL" -> "19"
# ---------------------------------------------------------------------------
prep_hg002() {
    info "=== Preparing HG002 for SVP ==="

    local ref="$SAMPLES_DIR/GRCh38.fa"
    local asm="$SAMPLES_DIR/hg002v1.1.fasta"

    [[ -f "$ref" ]] || die "GRCh38 not found. Run: ./download.sh grch38"
    [[ -f "$asm" ]] || die "HG002 assembly not found. Run: ./download.sh hg002-asm"

    local ref_out="$SAMPLES_DIR/GRCh38_svp.fa"
    local tumor_out="$SAMPLES_DIR/hg002_paternal_svp.fa"

    # Rename GRCh38 chromosomes: "chr19 AC:..." -> "19"
    if [[ ! -f "$ref_out" ]]; then
        info "  Renaming GRCh38 chromosomes to bare numbers..."
        python3 -c "
import sys
with open('$ref', 'r') as fin, open('$ref_out', 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            name = line[1:].split()[0]  # e.g. chr19
            if name.startswith('chr'):
                bare = name[3:]
            else:
                bare = name
            fout.write('>' + bare + '\n')
        else:
            fout.write(line)
"
        info "  Created: $ref_out"
    else
        info "  GRCh38 renamed FASTA exists, skipping."
    fi

    # Extract HG002 paternal haplotype, rename: "chr19_PATERNAL" -> "19"
    if [[ ! -f "$tumor_out" ]]; then
        info "  Extracting HG002 paternal haplotype (renamed to bare numbers)..."
        python3 -c "
import sys
with open('$asm', 'r') as fin, open('$tumor_out', 'w') as fout:
    writing = False
    for line in fin:
        if line.startswith('>'):
            name = line[1:].strip()
            if '_PATERNAL' in name:
                bare = name.replace('chr', '').replace('_PATERNAL', '')
                fout.write('>' + bare + '\n')
                writing = True
            else:
                writing = False
        elif writing:
            fout.write(line)
"
        info "  Created: $tumor_out"
    else
        info "  HG002 paternal FASTA exists, skipping."
    fi

    info "HG002 preparation complete."
    info "  Reference (SVP-ready): $ref_out"
    info "  Tumor (SVP-ready):     $tumor_out"
    info "  Normal = Reference (germline benchmark)"
}

# ---------------------------------------------------------------------------
# Step: Sort and index BAMs for rival callers
# ---------------------------------------------------------------------------
prep_rival_bams() {
    info "=== Preparing BAMs for Rival Callers ==="

    require_cmd samtools

    local hifi_dir="$SAMPLES_DIR/colo829_hifi"
    [[ -d "$hifi_dir" ]] || die "HiFi BAMs not found. Run: ./download.sh colo829-hifi"

    local ref="$SAMPLES_DIR/GRCh38.fa"
    [[ -f "$ref" ]] || die "GRCh38 not found. Run: ./download.sh grch38"

    local bam_dir="$SAMPLES_DIR/aligned_bams"
    mkdir -p "$bam_dir"

    require_cmd minimap2

    info "Aligning and sorting HiFi reads to GRCh38..."
    info "This will take several hours."

    # Align tumor
    local tumor_bam="$bam_dir/colo829_tumor.sorted.bam"
    if [[ ! -f "$tumor_bam" ]]; then
        info "  Aligning tumor reads..."
        # Merge tumor HiFi BAMs, convert to FASTQ, align, sort
        samtools cat "$hifi_dir"/m84039_230312_*.bam "$hifi_dir"/m84039_230328_*.bam \
            | samtools fastq -@ 4 - \
            | minimap2 -a -x map-hifi -t "$THREADS" "$ref" - \
            | samtools sort -@ 4 -o "$tumor_bam"
        samtools index -@ "$THREADS" "$tumor_bam"
        info "  Tumor aligned: $tumor_bam"
    else
        info "  Tumor BAM exists, skipping."
    fi

    # Align normal
    local normal_bam="$bam_dir/colo829_normal.sorted.bam"
    if [[ ! -f "$normal_bam" ]]; then
        info "  Aligning normal reads..."
        samtools cat "$hifi_dir"/m84039_230318_*.bam "$hifi_dir"/m84039_230324_*.bam \
            | samtools fastq -@ 4 - \
            | minimap2 -a -x map-hifi -t "$THREADS" "$ref" - \
            | samtools sort -@ 4 -o "$normal_bam"
        samtools index -@ "$THREADS" "$normal_bam"
        info "  Normal aligned: $normal_bam"
    else
        info "  Normal BAM exists, skipping."
    fi

    info "Aligned BAMs ready in: $bam_dir/"
    info "  Tumor:  $tumor_bam"
    info "  Normal: $normal_bam"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
usage() {
    cat <<'USAGE'
Usage: ./prepare.sh <step> [step2 ...]

Steps:
  index-ref           Index GRCh38 reference with samtools faidx
  prep-hg002          Rename GRCh38 + extract HG002 paternal haplotype for SVP
  assemble-colo829    Assemble COLO829 tumor+normal from HiFi BAMs (hifiasm)
  prep-rival-bams     Align HiFi reads to GRCh38, sort+index for rival callers

Shorthand:
  all                 Run all steps in order

Environment variables:
  THREADS=N           Number of threads (default: nproc)
USAGE
}

if [[ $# -eq 0 || "$1" == "--help" || "$1" == "-h" ]]; then
    usage
    exit 0
fi

steps=("$@")
if [[ "$1" == "all" ]]; then
    steps=(index-ref prep-hg002 assemble-colo829 prep-rival-bams)
fi

for step in "${steps[@]}"; do
    case "$step" in
        index-ref)         index_ref ;;
        prep-hg002)        prep_hg002 ;;
        assemble-colo829)  assemble_colo829 ;;
        prep-rival-bams)   prep_rival_bams ;;
        *)                 die "Unknown step: $step" ;;
    esac
    echo ""
done

info "Preparation complete."
