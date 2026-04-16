#!/usr/bin/env bash
# =============================================================================
# run_svp.sh — Run SVP on real-world datasets and compare against truth sets
#
# Usage:
#   ./run_svp.sh hg002                 # full HG002 benchmark (all chr 1-22)
#   ./run_svp.sh hg002 --chr 19        # single chromosome
#   ./run_svp.sh colo829               # COLO829 benchmark (needs assembly)
#   ./run_svp.sh colo829 --chr 7       # COLO829 single chromosome
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SVP_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SAMPLES_DIR="$SCRIPT_DIR/samples"
TRUTH_DIR="$SCRIPT_DIR/ground_truth"
RESULTS_DIR="$SCRIPT_DIR/results"

SVP="$SVP_ROOT/svp"
DEPTH="${DEPTH:-5}"
PARALLEL="${PARALLEL:-8}"

info()  { echo ":: $*"; }
warn()  { echo "!! $*" >&2; }
die()   { echo "ERROR: $*" >&2; exit 1; }

[[ -x "$SVP" ]] || die "SVP binary not found at $SVP. Run 'make' first."

mkdir -p "$RESULTS_DIR" "$SVP_ROOT/outputs/vcf" "$SVP_ROOT/outputs/gfa"

# ---------------------------------------------------------------------------
# Run SVP per-chromosome in parallel, merge VCF, report timing
# ---------------------------------------------------------------------------
run_svp_genome() {
    local ref="$1" normal="$2" tumor="$3" label="$4"
    shift 4

    # Parse --chr flag
    local single_chr=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --chr) single_chr="$2"; shift 2 ;;
            *) shift ;;
        esac
    done

    local chromosomes=()
    if [[ -n "$single_chr" ]]; then
        chromosomes=("$single_chr")
    else
        for c in $(seq 1 22); do chromosomes+=("$c"); done
    fi

    info "=== SVP Benchmark: $label ==="
    info "  ref    = $ref"
    info "  normal = $normal"
    info "  tumor  = $tumor"
    info "  depth  = $DEPTH"
    info "  chr    = ${chromosomes[*]}"
    info ""

    local time_log="$RESULTS_DIR/${label}.time.csv"
    echo "chr,time_sec,peak_mem_kb" > "$time_log"

    # Write a self-contained per-chromosome runner script
    local runner="/tmp/svp_runner_${label}.sh"
    cat > "$runner" <<RUNNER
#!/usr/bin/env bash
chr=\$1
out_name="${label}_chr\${chr}"
timefile="/tmp/svp_${label}_chr\${chr}.time"
/usr/bin/time -f "%e %M" -o "\$timefile" \\
    "$SVP" --ref "$ref" --normal "$normal" --tumor "$tumor" \\
           --depth "$DEPTH" --chr "\$chr" --output "\$out_name" 2>/dev/null
if [[ -f "\$timefile" ]]; then
    read -r t m < "\$timefile"
    echo "\${chr},\${t},\${m}" >> "$time_log"
    echo "    chr\${chr}: \${t}s, \$((m / 1024))MB"
fi
RUNNER
    chmod +x "$runner"

    if [[ ${#chromosomes[@]} -eq 1 ]]; then
        bash "$runner" "${chromosomes[0]}"
    else
        info "  Running ${#chromosomes[@]} chromosomes ($PARALLEL parallel)..."
        printf '%s\n' "${chromosomes[@]}" | xargs -P "$PARALLEL" -I{} bash "$runner" {}
    fi

    # Merge per-chromosome VCFs
    local merged="$RESULTS_DIR/${label}.vcf"
    {
        echo "##fileformat=VCFv4.2"
        echo "##source=svp"
        echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">'
        echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
        echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of variant">'
        echo '##INFO=<ID=INSSEQ,Number=1,Type=String,Description="Inserted sequence">'
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        for chr in "${chromosomes[@]}"; do
            local vcf="$SVP_ROOT/outputs/vcf/${label}_chr${chr}.vcf"
            [[ -f "$vcf" ]] && grep -v "^#" "$vcf" 2>/dev/null || true
        done
    } > "$merged"

    # Summary
    local total_time=0 peak_mem=0
    for chr in "${chromosomes[@]}"; do
        local tf="/tmp/svp_${label}_chr${chr}.time"
        if [[ -f "$tf" ]]; then
            read -r t m < "$tf"
            total_time=$(awk "BEGIN{print $total_time + $t}")
            [[ "$m" -gt "$peak_mem" ]] 2>/dev/null && peak_mem=$m
        fi
    done

    local call_count
    call_count=$(grep -vc "^#" "$merged" 2>/dev/null || echo 0)

    info ""
    info "  Results: $merged ($call_count calls)"
    info "  Total CPU time: ${total_time}s"
    info "  Peak memory: $((peak_mem / 1024))MB"
    info "  Timing log: $time_log"

    echo "$merged"
}

# ---------------------------------------------------------------------------
# Compare output against truth set
# ---------------------------------------------------------------------------
compare_results() {
    local svp_vcf="$1" truth_vcf="$2" ref="$3" label="$4"

    info ""
    info "=== Accuracy: $label ==="

    if [[ -f "$SVP_ROOT/scripts/compare_vcf.py" ]]; then
        info "--- Position-tolerant diff (50bp) ---"
        python3 "$SVP_ROOT/scripts/compare_vcf.py" "$truth_vcf" "$svp_vcf" 2>&1 || true
    fi

    if command -v truvari &>/dev/null; then
        info "--- Truvari benchmark ---"
        local truvari_out="$RESULTS_DIR/truvari_${label}"
        bash "$SVP_ROOT/scripts/truvari.sh" "$truth_vcf" "$svp_vcf" "$ref" "$truvari_out" 2>&1 || true
    fi
}

# ---------------------------------------------------------------------------
# HG002 — full germline SV benchmark
# ---------------------------------------------------------------------------
run_hg002() {
    local ref="$SAMPLES_DIR/GRCh38_svp.fa"
    local tumor="$SAMPLES_DIR/hg002_paternal_svp.fa"

    [[ -f "$ref" ]]   || die "Prepared GRCh38 not found. Run: ./prepare.sh prep-hg002"
    [[ -f "$tumor" ]] || die "Prepared HG002 not found. Run: ./prepare.sh prep-hg002"

    info "HG002 germline benchmark: GRCh38 as normal, HG002 paternal as tumor."
    info "Detects germline structural variants relative to the reference."
    info ""

    local svp_vcf
    svp_vcf=$(run_svp_genome "$ref" "$ref" "$tumor" "hg002" "$@")

    # Compare against GIAB truth set
    local truth="$TRUTH_DIR/hg002/HG002_SVs_Tier1_v0.6.vcf.gz"
    if [[ -f "$truth" ]]; then
        warn "NOTE: GIAB v0.6 uses GRCh37 coordinates — results need liftover for exact matching."
        compare_results "$svp_vcf" "$truth" "$ref" "hg002"
    else
        warn "HG002 truth set not found. Run: ./download.sh hg002-truth"
        info "SVP results: $svp_vcf"
    fi
}

# ---------------------------------------------------------------------------
# COLO829 — somatic SV benchmark (requires assembled FASTA)
# ---------------------------------------------------------------------------
run_colo829() {
    local ref="$SAMPLES_DIR/GRCh38.fa"
    local tumor="$SAMPLES_DIR/colo829_tumor.fa"
    local normal="$SAMPLES_DIR/colo829_normal.fa"

    [[ -f "$ref" ]]    || die "GRCh38 not found. Run: ./download.sh grch38"
    [[ -f "$tumor" ]]  || die "COLO829 tumor assembly not found. Run: ./prepare.sh assemble-colo829"
    [[ -f "$normal" ]] || die "COLO829BL normal assembly not found. Run: ./prepare.sh assemble-colo829"
    [[ -f "${ref}.fai" ]] || die "GRCh38 not indexed. Run: ./prepare.sh index-ref"

    local svp_vcf
    svp_vcf=$(run_svp_genome "$ref" "$normal" "$tumor" "colo829" "$@")

    local truth=""
    for f in "$TRUTH_DIR/colo829/truthset_somaticSVs_COLO829_hg38lifted.vcf" \
             "$TRUTH_DIR/colo829/truthset_somaticSVs_COLO829.vcf"; do
        [[ -f "$f" ]] && truth="$f" && break
    done

    if [[ -n "$truth" ]]; then
        compare_results "$svp_vcf" "$truth" "$ref" "colo829"
    else
        warn "COLO829 truth set not found. Run: ./download.sh colo829-truth"
        info "SVP results: $svp_vcf"
    fi
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
usage() {
    cat <<'USAGE'
Usage: ./run_svp.sh <dataset> [--chr N]

Datasets:
  hg002     HG002 T2T assembly — germline SV benchmark (primary)
  colo829   COLO829 assembly — somatic SV benchmark (needs HiFi assembly)

Options:
  --chr N   Run only chromosome N (default: all 1-22)

Environment:
  DEPTH=N      LCP depth (default: 5)
  PARALLEL=N   Parallel chromosome jobs (default: 8)

Examples:
  ./run_svp.sh hg002                # full genome, all 22 chromosomes
  ./run_svp.sh hg002 --chr 19       # single chromosome
  ./run_svp.sh colo829 --chr 7      # COLO829 chromosome 7
USAGE
}

if [[ $# -eq 0 || "$1" == "--help" || "$1" == "-h" ]]; then
    usage
    exit 0
fi

dataset="$1"; shift

case "$dataset" in
    hg002)    run_hg002 "$@" ;;
    colo829)  run_colo829 "$@" ;;
    *)        die "Unknown dataset: $dataset" ;;
esac
