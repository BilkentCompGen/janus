#!/usr/bin/env bash
# Self-assembly benchmark: lasm and rust-mdbg on the same reads, then QUAST against the GT assembly.
#
# Both assemblers are reference-free; the reference enters ONLY in the evaluation, never in the
# assembly. Wall time and peak RSS come from /usr/bin/time -v so the speed claim is measured.
#
# usage: asm_bench.sh <reads.fastq[.gz]> <outdir> <gt.fa> [threads=16] [tools=both]
set -euo pipefail

READS=$1; OUT=$2; GT=$3; THREADS=${4:-16}; TOOLS=${5:-both}
LASM=/home/yuio/svp/lib/lasm/build/lasm
MDBG_DIR=/home/yuio/svp/lib/rust-mdbg
export PATH=$HOME/.local/bin:$HOME/.cargo/bin:$PATH
mkdir -p "$OUT"

t() {   # t <tag> <cmd...>  -- run under /usr/bin/time -v, keep the log
    local tag=$1; shift
    echo "[bench] $tag: $*" >&2
    /usr/bin/time -v -o "$OUT/$tag.time" "$@" > "$OUT/$tag.log" 2>&1
}

gfa2fa() {  # segments of a GFA as FASTA, skipping the ones with no sequence
    awk -v pat="${3:-.}" '$1=="S" && $3!="*" && $2 ~ pat {print ">"$2"\n"$3}' "$1" > "$2"
}

if [[ $TOOLS == both || $TOOLS == lasm ]]; then
    # The bloom filter's projected count defaults to 1e8, which is chromosome-scale. Whole genome needs
    # it scaled or the filter saturates and stops discriminating; LASM_BLOOM is that knob.
    # Cap it at INT_MAX: BlockedAtomicBloomFilter takes `int n`, so 3e9 wraps negative and the
    # allocation throws std::bad_alloc -- which lasm survives, silently running with no filter at all.
    t lasm "$LASM" hifi "$READS" -o "$OUT/lasm.gfa" -l "${LASM_L:-4}" \
        --density "${LASM_DENSITY:-0.05}" -t "$THREADS" -@ 4 --file-reading-mode sequential \
        --projected-element-count "${LASM_BLOOM:-100000000}"
    # M* segments are lasm's strand-resolved layout nodes; C* are not, and carry both strands of the
    # same locus, so scoring them would double-count. M-only is the assembly, C is reported alongside.
    gfa2fa "$OUT/lasm.gfa" "$OUT/lasm.fa"    '^M'
    gfa2fa "$OUT/lasm.gfa" "$OUT/lasm_C.fa"  '^C'
fi

if [[ $TOOLS == both || $TOOLS == mdbg ]]; then
    t mdbg "$MDBG_DIR/target/release/rust-mdbg" "$READS" \
        -k "${MDBG_K:-21}" -l "${MDBG_L:-14}" --density "${MDBG_D:-0.003}" \
        --minabund 2 --threads "$THREADS" --prefix "$OUT/mdbg"
    # magic_simplify is where mdbg turns the minimizer-space graph into base-space contigs
    t mdbg_simplify bash "$MDBG_DIR/utils/magic_simplify" "$OUT/mdbg"
    cp "$OUT/mdbg.msimpl.fa" "$OUT/mdbg.fa"
fi

echo "[bench] assemblies:" >&2
for f in "$OUT"/lasm.fa "$OUT"/mdbg.fa; do [[ -s $f ]] && echo "  $f $(stat -c%s "$f") bytes" >&2; done
