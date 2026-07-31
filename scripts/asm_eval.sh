#!/usr/bin/env bash
# Score self-assemblies against the GT assembly with QUAST.
#
# The assemblers are reference-free; the reference is used HERE only, to say what came out.
# --min-identity 90 / --extensive-mis-size 1000 are QUAST's large-genome defaults; misassemblies are
# what matter for the SV plan, since a misassembly is a false breakpoint.
#
# usage: asm_eval.sh <gt.fa> <outdir> <asm1.fa> [asm2.fa ...]
set -euo pipefail

GT=$1; OUT=$2; shift 2
QUAST=/home/yuio/tmp/asmeval/bin/quast.py
export PATH=$HOME/.local/bin:$PATH
# QUAST 5.2 passes its config through module globals, which loky's default spawn start method does not
# inherit -- the workers then see min_contig=None and die. Forking keeps the globals.
# QUAST 5.2 hands its config to workers through module globals, which joblib's process pool does not
# inherit -- the workers then see min_contig=None and die. Running joblib in-process keeps the globals;
# the aligner subprocesses still get all --threads, which is where the time actually goes.
export JOBLIB_MULTIPROCESSING=0
mkdir -p "$OUT"

LABELS=$(for f in "$@"; do basename "$f" .fa; done | paste -sd,)
"$QUAST" -r "$GT" -o "$OUT" -t "${THREADS:-16}" --labels "$LABELS" \
    --min-identity 90 --extensive-mis-size 1000 --min-contig 500 --no-icarus --no-html \
    ${QUAST_EXTRA:-} "$@" 2>&1 | tail -20

cat "$OUT/report.txt"
