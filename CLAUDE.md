# CLAUDE.md

Guidance for Claude Code working in this repository.

## Project

`svp` is a somatic structural-variant caller. Active branch: **`graph-inspector`** (this branch).
It builds a **2-colour LCP-core graph directly from reads** (normal vs tumor HiFi FASTQ) and
emits candidate somatic bubbles. No assembly, no alignment, no reference graph — a single
binary. Output is **coordinate-free** today (per-node multiplicity, not genome positions);
mapping back to coordinates / a VCF is the next step.

This branch adds the `inspect/` graph-inspector toolkit and the in-tree Mode A unitig
compaction (`src/modea.cpp` → `svp modea`) on top of the `lcp-cortex-merged` caller. The
caller-only branch (`lcp-cortex-merged`) is kept clean for the merge base.

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
| `lcp-cortex-merged` | merge base for this branch — fast gzip+LCP plumbing from `lcp-cortex-fast`, simple BFS connected-components bubble decomposition from the 2026-05-26 high-recall version. ~99.6% recall / ~78% precision |
| `graph-inspector` | **this branch** — `inspect/` toolkit (seqs_probe / bubble_census / inspect + Python viz), Mode A unitig compaction in `src/modea.cpp`, plus the phantom-bubble and RC-pairing investigation tools in `tools/` |
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

## CLI — four resumable, checkpointed stages

The expensive LCP pass is paid once; later stages are seconds and re-runnable.

```
svp lcp     <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5 --flank 10 --margin 10 --min-mult 3 --threads 16]
svp graph   <dir>
svp modea   <dir>
svp bubbles <dir> [--out <dir>/bubbles.tsv]
svp run     <normal.fq.gz> <tumor.fq.gz> <dir> [all opts]   # runs all, skipping done stages
```

`--min-mult` is an `lcp`-stage flag (the prune runs inside `run_lcp` before `counts.bin` is
written), NOT a `graph` flag — re-tuning it requires re-running `lcp`. `run` skips any stage
whose checkpoint exists (`counts.bin`+`edges.bin`, `graph.bin`, `modea.bin`). Checkpoints are
written atomically (`.tmp` + `rename`): an interrupted write never leaves a half-file that
resume would mistake for complete.

## Pipeline (`src/`, 1001 LOC across 6 cpp files: graph 517 / modea 199 / cores 82 / main 69 / MurmurHash3 56 / io 28)

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
3. **modea** (`run_modea`, `modea.cpp`) — Mode A unitig compaction on the colour-1/2 (non-
   shared) subgraph. Walks each maximal same-colour chain (in_deg/out_deg constraints), splices
   the inline node bytes via exact KMP suffix/prefix overlap (k==0 is FATAL — implies upstream
   data is broken), and emits one super-node per chain with its per-position coverage runs.
   Pure cycles are preserved as singleton super-nodes ([[modea_preserves_cycles]]). → `dir/modea.bin`.
   Mode B (allele-diff refinement on top of Mode A) lives in `inspect/`, not in `src/`.
4. **bubbles** (`run_bubbles`, `graph.cpp`) — **plain BFS connected-components on non-shared
   nodes**. Backbone (colour=3) is NEVER traversed; each BFS component is one bubble; type
   from the colour bits seen during the walk (`D` = only colour-1, `I` = only colour-2, `C`
   = both). No anchor-merging / union-find — that variant over-merges adjacent loci. →
   `dir/bubbles.tsv` (`type  n_nodes  normal_cov  tumor_cov`). Seconds.

`src/svp.h` holds the one shared POD (`Core`) + the LCP primitives. `extract_cores` /
`compute_sigs` are used by `tools/` only (the caller uses `extract_fast` directly for speed);
`io.cpp::read_fasta_flat` is for the assembly-based GT/analysis tools only — never the caller.

Relevant docs (read on demand, not all at once):
- `docs/PERFORMANCE.md` — gzip+lcptools speed details
- `docs/BUBBLE_CALLER_COMPARISON.md` — bubble-algorithm history (assembly chr19 → reads WG)
- `docs/GRAPH_MODELS.md` — Mode A / Mode B compaction design
- `docs/INSPECTOR_FINDINGS.md` — first-run findings from the inspector
- `docs/VARIANT_CALLING_RESEARCH.md` — Severus/Sniffles2 audit, posmap/resolve/vcf plan
- `docs/READ_ERROR_QC.md`, `docs/LCP_TOPOLOGY.md`, `docs/CORE_ANALYSIS.md` — read-mode QC, LCP probes, core overlap
- `docs/VIZ_RESEARCH.md` — visualization background for the inspector

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

## inspect/ — graph-inspector toolkit (this branch only)

Diagnostic + visualization for the 2-colour graph. Separate `inspect/Makefile`; the parent
build supplies `build/io.o`, `build/cores.o`, `build/MurmurHash3.o` (auto-invoked).

```
cd inspect && make            # builds seqs_probe, bubble_census, inspect
./inspect/seqs_probe   <dir>                            # Stage 0: confirm seqs.bin schema
./inspect/bubble_census <dir> <truth.vcf> <assets>      # Stage 1: WG sweep (~10 min) -> census.{bin,tsv}
./inspect/inspect      <dir> {sv|bubble|region} <key>   # Stage 2: per-subgraph raw/modeA/modeB.gfa + meta.json
python3 inspect/src/inspect_viz.py <dir>/inspect/<key>  # Stage 3: PNGs + index.html
```

Mode B (allele-diff refinement) is implemented inside `inspect`, not `src/` — it extends
anchor super-nodes by shared prefixes/suffixes in place ([[modeB_on_modeA]]). Process every
bubble including pathological ones — no `MAX_INTERNAL` truncation ([[no_artificial_caps]]).

## Code style — strict

Read `docs/CODE_STYLE.md` before editing. C++17, no OOP, no user templates, no try/catch
(fail-fast `[FATAL]` + `std::abort()`), avoid if-else trees, ≤30 LOC/function. Comments only
when the WHY is non-obvious. `.o`/`.d` go to `build/`.

The historical "≤700 LOC across `src/*.cpp`" cap from the caller-only branch no longer holds
on `graph-inspector` (modea adds ~200 LOC); the spirit (keep `src/` tight; push diagnostics
to `tools/` or `inspect/`) still applies.

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
