# CLAUDE.md

Guidance for Claude Code working in this repository.

## Project

`svp` is a somatic structural-variant caller. Active branch: **`lcp-cortex-merged`** (this branch).
It builds a **2-colour LCP-core graph directly from reads** (normal vs tumor HiFi FASTQ) and
emits candidate somatic bubbles. No assembly, no alignment, no reference graph — a single
binary. Output is **coordinate-free** today (per-node multiplicity, not genome positions);
mapping back to coordinates / a VCF is the next step.

### Two hard rules — do not violate

1. **The caller consumes READS ONLY** — two FASTQ inputs (normal, tumor). An assembly/`.fa`
   is **never** a caller input. Assemblies *are* allowed as a ground-truth / analysis oracle
   in `tools/`: extract their LCP cores and match read-derived cores by signature. LCP is
   deterministic, so equal signature = exact match — no alignment needed.
2. **Whole-genome only.** Nothing chr19-specific in the caller. Sim scope is chosen by the
   simulation config, never hard-coded here.

## Branch map

| Branch | What's there |
|---|---|
| `main` | clean slate, scaffolding only |
| `lcp-cortex` | original LCP-anchor (match.cpp) pipeline at HEAD `ae0c50e` |
| `lcp-cortex-fast` | reads-only Cortex caller with fast gzip+LCP, but with the broken chain+union-find decomposition that hurt recall (kept for diff/reference) |
| `lcp-cortex-merged` | **this branch** — fast gzip+LCP plumbing from `lcp-cortex-fast`, simple BFS connected-components bubble decomposition from the 2026-05-26 high-recall version. Target: ~99.6% recall / ~78% precision (the user-described "high recall / lower precision" regime) |
| `vg-rewrite` | prior assembly-input LCP-anchor rewrite (pre-reads) |
| `backup` | previous monolithic LCP-graph impl. Source of `zstr.hpp`, `MurmurHash3.{h,cpp}` |
| `experimental/modular-pipeline` | `svp2/` 6-stage modular pipeline; entry: `svp2/STATUS.md` |
| `kmer`, `no-compression` | dead-end experiments |

Don't port files across branches without explicit instruction.

## Build & run

```
make install-deps    # one-time: builds lib/lcptools (liblcptools.a)
make                 # builds ./svp
make run             # whole-genome read run: reads/normal vs reads/wgs -> outputs/svp/wgs/
make clean           # moves build/, svp, outputs/ -> ~/trash/svp-make-clean/
```

`make run` overrides: `LEVEL` (lcptools depth, 5), `FLANK` (±bp signature window, 10),
`MINMULT` (error-prune threshold, 3), `THREADS` (16), `READS_N`/`READS_T` (FASTQ paths).

## CLI — three resumable, checkpointed stages

The expensive LCP pass is paid once; later stages are seconds and re-runnable.

```
svp lcp     <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5 --flank 10 --margin 10 --threads 16]
svp graph   <dir> [--min-mult 3]
svp bubbles <dir> [--out <dir>/bubbles.tsv]
svp run     <normal.fq.gz> <tumor.fq.gz> <dir> [all opts]   # runs all, skipping done stages
```

`run` skips any stage whose checkpoint already exists, so re-tuning `--min-mult` re-runs only
graph+bubbles (no re-LCP). Checkpoints are written atomically (`.tmp` + `rename`): an
interrupted write never leaves a half-file that resume would mistake for complete.

## Pipeline (`src/`, 546 LOC across 5 cpp files: graph 316 / cores 82 / main 64 / MurmurHash3 56 / io 28)

1. **lcp** (`run_lcp`, `graph.cpp`) — one reader thread per FASTQ shard (so both samples
   inflate in parallel CPUs) + N-2 workers consuming colour-tagged batches. The reader is a
   block-based gzread + manual `\n` scan (line-by-line `gzgets` is 3-5× slower because zlib
   state crosses every line + per-line `std::string` ops dominate). Per read: uppercase,
   call `extract_fast` (lcptools perf-session2 SoA API — inlined into `process_read`, no
   intermediate `Core[]`/`uint64_t[]` vectors), keep interior cores (≥`margin` bp from both
   ends), inline the flank-window `murmur64` for the signature. Workers emit `SigCol{sig,col}`
   records + within-read `(sig_prev, sig_curr)` edges to thread-local vectors; the pool
   merges via LSB radix sort + run-length reduce. Atomic-rename checkpoint:
   `dir/counts.bin` (16 B per record: `sig, n, t`), `dir/edges.bin` (16 B per record).
2. **graph** (`run_graph`, `graph.cpp`) — drop signatures with `max(n,t) < min_mult`
   (sequencing errors are multiplicity-1; `min_mult=3` removes 99.9% of them). Colour each
   survivor: 1 = normal-only, 2 = tumor-only (candidate somatic), 3 = shared. Map signature
   edges → node-id edges, dedupe. → `dir/graph.bin`. Seconds.
3. **bubbles** (`run_bubbles`, `graph.cpp`) — **plain BFS connected-components on non-shared
   nodes**. Backbone (colour=3) is NEVER traversed; each BFS component is one bubble; type
   from the colour bits seen during the walk (`D` = only colour-1, `I` = only colour-2, `C`
   = both). No anchor-merging / union-find — that variant over-merges adjacent loci. →
   `dir/bubbles.tsv` (`type  n_nodes  normal_cov  tumor_cov`). Seconds.

`src/svp.h` holds the one shared POD (`Core`) + the LCP primitives. `extract_cores` /
`compute_sigs` are used by `tools/` only (the caller uses `extract_fast` directly for speed);
`io.cpp::read_fasta_flat` is for the assembly-based GT/analysis tools only — never the caller.

Read `docs/PERFORMANCE.md` for the gzip+lcptools speed details, and
`docs/BUBBLE_CALLER_COMPARISON.md` for the bubble-algorithm history (assembly chr19 → reads WG).

### Measured WG performance & recall (lcp-cortex-merged, AMD 9950X, 8 threads)

| metric | value |
|---|---|
| inputs | `simulation/reads/{normal,wgs}/hifi.fastq.gz` (29 GB each, 30× HiFi) |
| wall clock | **12:02** (full lcp+graph+bubbles, 8 threads) |
| CPU usage | 466% (gzip-bound at 2 streams; one reader per file) |
| peak RAM | 129 GB |
| normal reads / tumor reads | 6,176,526 / 6,174,744 |
| distinct sigs (post-LCP) / edges | 489,021,068 / 628,505,021 |
| nodes after `min_mult=3` / node-edges | 85,309,098 / 86,145,767 |
| bubbles total | **112,977** (57,001 I-somatic + 55,976 D-normal-only + 0 C-mixed) |

GT oracle = `simulation/ref.fasta` (germline) vs `simulation/truth/wgs/out.fa` (insilicoSV tumor):

| metric | value |
|---|---|
| truth somatic sigs \|S\| | 96,903 |
| called somatic sigs \|C\| (n=0, t≥3) | 122,948 |
| \|C ∩ S\| | 96,532 |
| **sig RECALL** | **99.62%** |
| **sig PRECISION** | **78.51%** |
| **bubble PRECISION** (I-bubbles with ≥1 sig ∈ S) | **63.00%** (35,911 / 57,001) |

Tool: `tools/bubble_recall <counts.bin> <graph.bin> <normal.fa> <tumor.fa>` — re-derives
node→sig from `counts.bin`, replays the BFS decomposition, and reports both sig-level (the
union-of-I-nodes metric, identical to `tools/recall`) and bubble-level (per-component) TP/FP.
Runs in ~5 min on WG.

**Why bubble-precision (63%) trails sig-precision (78.5%):** residual-error sigs that survive
`min_mult=3` cluster into small (1-2 node) singleton bubbles. Each such cluster bumps the
bubble denominator by 1 even though the sig count of the cluster is tiny. The fix is upstream
(more aggressive pruning / reverse-complement signature canonicalization to fold strand-pair
sigs into one node, halving the bubble count) — NOT a chain-merge that over-fuses real loci.

**Next:** map nodes back to genome positions for recall vs a truth VCF, and trim bubble allele
sequences to call simple SNP/INS/DEL — both may use truth assemblies as the GT oracle (rule 1).

## Code style — strict

Read `docs/CODE_STYLE.md` before editing. C++17, no OOP, no user templates, no try/catch
(fail-fast `[FATAL]` + `std::abort()`), avoid if-else trees, ≤30 LOC/function, ≤700 LOC total
across `src/*.cpp`. Comments only when the WHY is non-obvious. `.o`/`.d` go to `build/`.

## Simulation (reads — the only caller input)

`simulation/Makefile` is independent (delegated by the top-level Makefile); layout is flat.
Reads are generated whole-genome:

- `scripts/sim_reads_wg_ccs.sh` — realistic 30x HiFi: pbsim3 multi-pass CLR → real `ccs`
  consensus (correlated, context-dependent errors). → `reads/wgs/`, `reads/normal/`.
- `scripts/sim_reads_wg.sh` / `tools/sim_hifi` — fast direct Q30 injector with exact known
  errors; for QC, not realism.
- Truth assemblies/VCF (`simulation/truth/<id>/`, insilicoSV) — GT/analysis only.

**Long runs survive disconnect** in a detached `screen` (the realistic ccs sim is multi-day):
`screen -dmS <name> -L -Logfile <log> <cmd>` (`-L -Logfile` avoids `>` redirects).

## Dependencies

- `lib/lcptools` — LCP core extraction (`make install-deps`). Static lib is stale-prone;
  rebuild if you see segfaults on struct fields after pulling `lib/`.
- `lib/pbsim3` — long-read simulator (autotools; `autoreconf -fi` for `config.guess/sub`).
- `lib/ccs` — PacBio CCS v6.4.0 static binary.
- `lib/insilicoSV` — truth simulator. `zlib` — gz FASTA/FASTQ IO. `truvari` — VCF
  benchmarking (once a VCF emitter exists).

`lib/`, `simulation/`, `outputs/`, `build/`, `tools/`, `svp` are gitignored.
