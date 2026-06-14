# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`svp` is a **reference-free somatic structural-variant caller**. It builds a 2-colour
LCP-core graph directly from reads (normal vs tumor HiFi FASTQ) and progressively
simplifies it: raw sig graph → unitig-compacted super-graph (Mode A) → bubble-zipped
super-graph at base-pair resolution (Mode B). Mode A and Mode B are **one `modeb` stage**
— Mode A is built in RAM and consumed in place, never written to disk. No assembly, no
alignment, no reference graph — a single binary. Output is **coordinate-free** (per-node
sequence + per-base coverage, not genome positions).

Active branch: **`graph-inspector`**. The pipeline currently ends at **Mode B** — the
downstream VCF-emission tail (`posmap`/`anchor`/`vcf`), the `inspect/` toolkit, and the
exploratory `tools/` were removed (see git history `b68a980`); the repo is intentionally
trimmed to the read → graph → Mode B core.

### Two hard rules — do not violate

1. **The caller consumes READS ONLY** — two FASTQ inputs (normal, tumor). An
   assembly/`.fa` is **never** a caller input. Assemblies are allowed only as a
   ground-truth / analysis oracle (extract their LCP cores, match
   read-derived cores by signature — LCP is deterministic, so equal signature = exact
   match, no alignment).
2. **Whole-genome only.** Nothing chromosome-specific in the caller. Sim scope is chosen
   by the simulation config, never hard-coded.

## Build & run

```
make install-deps    # one-time: builds lib/lcptools (liblcptools.a) on the biolcp branch
make                 # builds ./svp from src/*.cpp
make run             # whole-genome read run: reads/normal vs reads/wgs -> outputs/svp/wgs/
make clean           # moves build/, svp, outputs/ -> ~/trash/svp-make-clean/
```

`make run` overrides: `LEVEL` (lcptools depth, 5), `FLANK` (±bp signature window, 10),
`MINMULT` (error-prune threshold, 3), `THREADS` (16), `READS_N`/`READS_T` (FASTQ paths).

**Build flag `-DLARGE_GENOME` is load-bearing and must match `lib/lcptools`'s build**
(64-bit LCP core start/end on a whole-genome reference). `make install-deps-lcptools`
archives lcptools in-place with the same flag. Rebuild the static lib if you see
segfaults on struct fields after pulling `lib/`.

## CLI — four resumable, checkpointed stages

The expensive LCP pass is paid once; later stages are seconds-to-minutes and re-runnable.
`run` skips any stage whose checkpoint already exists. Checkpoints are written atomically
(`.tmp` + `rename`).

```
svp lcp     <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5 --flank 10 --margin 10 --min-mult 3 --threads 16]
svp graph   <dir>
svp modeb   <dir>     # Mode A unitig compaction (in RAM) then Mode B bubble-zip -> modeb.bin
svp bubbles <dir> [--out <dir>/bubbles.tsv]
svp run     <normal.fq.gz> <tumor.fq.gz> <dir> [all opts]   # runs all, skipping done stages
```

`--min-mult` is an **`lcp`-stage** flag (the prune runs inside `run_lcp` before
`counts.bin` is written), NOT a `graph` flag — re-tuning it requires re-running `lcp`.
`pre_prune.bin` (the per-sample-filtered merged sig table) is written by `lcp`; if it
exists, re-running `lcp` resumes from it and skips the gzip+pool+merge stages.

**Schema changes invalidate old checkpoints**: a `modeb.bin` from before a struct change
FATALs cleanly on a size mismatch — regenerate the stage.

## Pipeline (`src/`, ~1630 LOC across 11 files; graph 639 / modeb 410 / modea 182)

A FASTQ input is a single file or a directory of shards. Stages and their on-disk schemas:

1. **lcp** (`run_lcp`, `graph.cpp`) — reads → per-sig multiplicity + within-read edges.
   A pool of reader threads (work-steal shards, block gzread + manual `\n` scan) feed a
   bounded queue drained by worker threads. **Per read, two LCP extractions**: forward,
   and reverse-complement (remap+reverse through a 256-byte complement LUT). Each surviving
   core's signature = `murmur64(core ± flank bytes, seed 0)`. Runs **per-sample two-pass**
   (normal then tumor; halves peak pool RAM), each sample filtered at `count >= min_mult`
   **before** the merge. A second pass (`fill_pool`) re-LCPs to fill each kept sig's window
   bytes and emit per-worker-deduped adjacent-sig edges.
   → `dir/counts.bin` (`u64 N`, `N × SigCount{u64 sig, u32 n, u32 t}`, `u64 total`,
   `N × u16 seqlen`, `payload` = the ±flank window bytes per sig) + `dir/edges.bin`
   (`u64 M`, `M × IEdge{u32,u32}` sig pairs) + `dir/pre_prune.bin` (resume checkpoint).
2. **graph** (`run_graph`, `graph.cpp`) — colour each sig: **1 = normal-only, 2 =
   tumor-only (candidate somatic), 3 = shared**. Map sig-edges → node-id edges, dedupe.
   → `dir/graph.bin` (`u64 K`, `color[K]`, `ncov[K]` u32, `tcov[K]` u32, `u64 E`,
   `IEdge[E]`). Seconds.
3. **modeb** (`run_modeb`, `modeb.cpp`) — the single graph-simplification stage. It first
   builds the **Mode A** super-graph **in RAM** (`build_super_graph`, `modea.cpp`) — unitig
   compaction: walk each maximal same-colour chain (in/out-degree 1), splice the inline node
   bytes via exact KMP suffix/prefix overlap (`suffix_prefix`; `k==0` is FATAL — adjacent L5
   cores always overlap ~85 bp), emit one super-node per chain with its per-base coverage
   runs; pure cycles become singleton super-nodes. The node schema is
   `SuperNode{u32 block_len, u32 runs_len, u8 color}[S]` + concatenated `blocks` (bytes) +
   `CovRun{u32 end_pos, n_cov, t_cov}[]` + `super_edges` (`IEdge[]`); per-node slices are
   found by prefix-summing `block_len`/`runs_len`. **This Mode A graph is never written to
   disk** — modeb consumes it directly.
   Then **Mode B**: emit the whole graph, further simplified (it does NOT extract per-bubble
   records — that was the old design). Every Mode A super-node passes through; each
   `std_both` bubble (a `(S1,S2)` colour-3 anchor pair with parallel c1 and c2 arms) is
   **zipped to base-pair resolution**: the prefix/suffix common across **all** of an anchor's
   bubble arms is absorbed into the anchor (grown), and each arm becomes one residual
   super-node holding just its divergent middle. EQ bubbles (identical alleles) collapse to a
   direct edge. Coverage is max-merged where alleles share a region. Output schema = the Mode
   A node schema + edges that carry an **explicit overlap**: `GEdge{u32 from, u32 to, u32
   overlap}` — reconstruct a path by dropping `overlap` head bytes of `to` before
   concatenating. **Bubble seams abut (overlap 0)** — the two arms diverge at exactly one base
   so no overlap can span it; pass-through edges keep the Mode A suffix/prefix overlap. An
   anchor is grown only when every edge on the grown side is a zipped arm (else growing its
   tail would corrupt a surviving backbone edge's overlap). **v1 scope: `std_both` only** —
   orphan buckets (tumor-only INS / normal-only DEL) and tangles (multi-path, cycles) pass
   through unzipped. → `dir/modeb.bin` (node schema above + `GEdge[]`).
   *Implementation notes:* arm coverage is kept RLE (`CovRun`), not flat per-base `uint32`,
   and only the slices actually read are expanded — this is what keeps peak RSS bounded.
   Adjacency is CSR; the bubble BFS/DFS use generation-stamped scratch (no per-call O(S)
   alloc).
4. **bubbles** (`run_bubbles`, `graph.cpp`) — plain BFS connected-components on the
   non-shared nodes of `graph.bin` (backbone colour-3 never traversed). Each component is
   one bubble; type `D` (colour-1 only), `I` (colour-2 only), `C` (both). Operates on
   `graph.bin`, **independent of modeb**. → `dir/bubbles.tsv`.

`src/svp.h` holds the shared POD (`Core`, `SigCount`, `IEdge`) + LCP primitives +
checkpoint IO helpers (`xopen`/`wr`/`rd`/`atomic_close`). `extract_cores`/`compute_sigs`
(`cores.cpp`) and `read_fasta_flat` (`io.cpp`) serve the assembly-based GT/analysis tool
only — never the caller, which uses `extract_fast` (`extract_fast.cpp`, the biolcp adapter)
directly.

## Load-bearing facts (verify in code before changing)

- **The ±flank window is load-bearing.** The signature is `murmur64(core ± 10bp flank)`,
  not the core alone — the flank makes each node locus-specific, so the stored window
  bytes are consistent across reads and the modea/modeb byte-splice works. Measured:
  dropping the flank craters recall and produces zero-overlap splice seams.
- **Dual fw+RC extraction → two nodes per locus.** Every locus appears as both `W` and
  `RC(W)`, each carrying ~full both-strand coverage. This pools strand coverage (lets
  `min_mult` go high to prune ~90% of error sigs) but **doubles** the node/bubble count
  and makes inversions invisible as somatic: the inverted *body* becomes colour-3 (shared
  with normal's RC pass); only the breakpoint *junction* cores are tumor-only, and those
  are novel sequence so they have no reference position. Folding strand pairs is a known
  future step.
- **Coverage representation.** Each super-node carries a per-base `(normal, tumor)` depth
  track, run-length-encoded as `CovRun{end_pos, n_cov, t_cov}` (a run covers
  `[prev_end, end_pos)`). Overlaps between constituent nodes resolve by **max per colour**.
  Both `n` and `t` live in one run because a shared (c3) base genuinely has both depths;
  on c1/c2 nodes one field is always 0 (redundant with colour, kept for the c3/merge case).
- **`margin` guard.** A core must sit ≥`margin` from both read ends. With the default
  `margin == flank` the ±flank window is always full and in-bounds. Footgun: `flank >
  margin` truncates edge-proximal windows → the same locus gets two signatures → split
  node. Keep `margin >= flank`.
- **Per-sample `min_mult` gate is the central precision/recall knob.** Each sample is
  filtered independently *before* the merge, so a shared locus with normal count below
  `min_mult` is recorded tumor-only (somatic). More aggressive at higher `min_mult`.

## Known gaps / next steps

- **modeb v1**: orphan buckets and tangles pass through unzipped — orphan (INS/DEL)
  zipping against the backbone detour is the obvious v2.
- **Mode A build peak**: the merged `modeb` peak is now entirely the in-RAM Mode A build,
  which holds `counts.bin` payload + output blocks simultaneously (≈15 GB + 5.82 GB at WG).
  `mmap`-ing the payload region (read-only) would make it evictable page cache and lift the
  hard allocation ceiling for COLO829 full-coverage — output stays byte-identical.
- **No sequence-reconstruction verifier** for the Mode B output graph yet (only structural
  accounting was checked). A reader that splices the output back and compares to the
  alleles would be the hard proof.
- **`run_lcp` gzread bug** (`graph.cpp:216`, `:345`): `gzread` error (`-1`) is treated as
  EOF, so a corrupt/truncated `.fastq.gz` is silently accepted — add a `got < 0` FATAL.
- **Dead code in `graph.cpp`**: `radix_sort_u64`, `reduce_sigs_u64`, `prune_by_mult` are
  defined but never called (the pool counts in a hashmap and filters in `run_pool`).

### Measured (AMD 9950X)

- `lcp` (WG, dual fw+RC, 8 threads): ~24 min wall, ~96 GB peak RSS; ~887 M raw sigs,
  ~85 M kept after the prune, ~86 M edges.
- Mode A build (`wgs_mm10`): 85,152,754 nodes → **931,910 super-nodes**, 85.9 M edges →
  1.69 M super-edges, 5.82 GB blocks.
- Mode B (`wgs_mm10`): 931,910 → 930,045 nodes, 37,654 bubbles (33,176 `std_both` zipped,
  4,478 passed), 64,491 residual arms, 1,861 EQ collapses.
- **modea+modeb merge** (smoke set, K=12.25 M, S=132,840; warm cache): the old two-stage
  path peaked at modeb **13.27 GB** RSS (flat per-base `uint32` arm coverage = 9.6 GB of it);
  the merged single `modeb` peaks at **4.27 GB** (now entirely the Mode A build) in ~34 s
  vs ~41 s, with no 884 MB `modea.bin` round-trip and byte-identical `modeb.bin`. The
  remaining ceiling is the Mode A build holding `counts.bin` payload + output blocks together
  (≈15 GB + 5.82 GB at WG) — `mmap`-ing the payload is the open next step for COLO829 scale.

## Code style — strict

Read `docs/CODE_STYLE.md` before editing. C++17, no OOP, no user templates, no try/catch
(fail-fast `[FATAL]` + `std::abort()`), avoid if-else trees, ≤30 LOC/function. Comments
only when the WHY is non-obvious. `.o`/`.d` go to `build/`. Keep `src/` tight.

## Dependencies

- `lib/lcptools` — LCP core extraction. **Must be on the `biolcp` branch** (DNA-optimised),
  built with `-DLARGE_GENOME`. `src/extract_fast.cpp` adapts biolcp's
  `init_lps_offset`/`lps_deepen` to svp's SoA `extract_fast`. See `docs/FOUNDATIONS.md`.
- `lib/Badread` (HiFi read simulator, venv at `lib/Badread/venv`),
  `lib/insilicoSV` (truth simulator). `zlib` for gz IO.
- `lib/`, `simulation/`, `outputs/`, `build/`, and `svp` are gitignored.

## Simulation (reads — the only caller input)

`simulation/Makefile` is independent (delegated by the top-level Makefile). Reads are
generated whole-genome. Long runs survive disconnect in a detached `screen`:
`screen -dmS <name> -L -Logfile <log> <cmd>`.

- `scripts/sim_reads_badread.sh` / `gen_badread_reads_wg.sh` — Badread PacBio HiFi over an
  insilicoSV genome (current realistic benchmark). Q30/stdev-3 (`--identity 30,3`) on the
  recent `pacbio2021` model; reads-only, contig-sharded across threads (Badread is
  single-threaded).
- Truth assemblies/VCF (`simulation/truth/<id>/`, insilicoSV) — GT/analysis only (rule 1).

## Branch map

| Branch | What's there |
|---|---|
| `graph-inspector` | **this branch** — read → graph → Mode B (Mode A built in RAM) core (VCF tail, inspect/, tools trimmed) |
| `lcp-cortex-merged` | merge base: fast gzip+LCP + simple BFS bubble decomposition; ~99.6% recall / ~78% precision |
| `lcp-cortex`, `lcp-cortex-fast` | earlier LCP-anchor pipelines (kept for diff/reference) |
| `main` | clean slate, scaffolding only |
| `backup` | previous monolithic LCP-graph impl (source of `zstr.hpp`, `MurmurHash3`) |
| `experimental/modular-pipeline` | `svp2/` 6-stage modular pipeline; entry `svp2/STATUS.md` |

Don't port files across branches without explicit instruction.

Docs (the repo was trimmed to bare minimum): `docs/PIPELINE.md` (current stages +
on-disk schemas + algorithms), `docs/FOUNDATIONS.md` (LCP cores, the signature, why the
flank is load-bearing, structural ceilings), `docs/CODE_STYLE.md`. Older analysis/VCF docs
were removed; their distilled conclusions live in those two summaries.
