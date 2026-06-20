# CLAUDE.md

_Project-specific guidance. Global behavioral rules + cross-project map: `/home/claude/CLAUDE.md`._

---

`svp` — reference-free somatic SV caller. Builds a 2-colour LCP-core graph straight from reads (normal vs tumor HiFi FASTQ), simplifies it, and emits a coordinate-free graph (per-node sequence + per-base coverage, no genome positions). C++17, one binary. No assembly, no alignment, no reference. Flow: **`run` = `graph → simplify`** (two stages, that's all). `min_mult=10`, no prune. State: `main` @ `eb0194e` **plus uncommitted working-tree changes** — `graph` decodes `.gz` with **rapidgzip** (parallel, in-process) and skips pass-2 LCP via a **cross-pass core-boundary cache**; the legacy `lcp`/`counts.bin` path, the `bubbles` diagnostic, and the orphaned assembly helpers have been removed.

## Hard rules
- **Reads only** (two FASTQ). An assembly/`.fa` is a ground-truth oracle, never a caller input.
- **Whole-genome only**; sim scope comes from the config, never hard-coded.

## Foundations
- **LCP core** (lcptools `biolcp`, level 5): a locally-consistent parse — the same locus yields the same core in any context. Adjacent L5 cores overlap ~85 bp, so nodes splice by exact KMP suffix/prefix (**zero overlap = FATAL**).
- **Signature** = `murmur64(core ± 10 bp flank, seed 0)`. The ±flank makes each sig locus-specific so a node's stored window bytes are consistent across reads (the splice depends on it); dropping it craters recall.
- **2 colours**: `(n && t) ? 3 : n ? 1 : 2` → **1 normal-only, 2 tumor-only (somatic), 3 shared**.
- **Dual fwd+RC** extraction per read (the parse isn't RC-symmetric): pools both-strand coverage (so `min_mult` prunes errors) but doubles node count and hides inversions as somatic. The two strands produce **independent** sigs (no per-emission canonicalization — that would break the KMP byte-splice). See `docs/FOUNDATIONS.md`.

## Pipeline
`run` = `graph → simplify`. Both stages are resumable — `run` skips a stage whose output (`graph.bin` / `simplify.bin`) already exists; every checkpoint is `.tmp` + `rename`. **Interactive map: `docs/pipeline.html`** (serve it; it diagrams both stages and the on-disk formats).

### `graph <normal.fq.gz> <tumor.fq.gz> <dir>` → `graph.bin`  (`run_graph_reads`)
Two passes over the same FASTQs, self-contained:
1. **Count** (`run_pool_both`): both samples in **one** pool. One pipe-reader per `.gz` reads plain FASTQ from a **rapidgzip** decode process (`open_decode` → `popen`); workers LCP each read (fwd + RC) → `emit_sig` → count in a sharded **bloom cascade**, one cascade per colour (`min_mult-1` plain bloom filters + exact survivor map). Per-colour `min_mult` gate happens **inside** the cascade, before the n/t merge (a normal-error sig drops, leaving the locus tumor-only). Bloom only over-keeps (FN = 0). Survivors sorted by sig → node id = sorted index (raw sigs are **not** persisted).
2. **Fill** (`fill_pool`): re-decode the same FASTQs; for each survivor fill its ±flank window bytes (atomic single-writer claim) and emit adjacent-signature edges (per-worker deduped, periodic radix-compact, final dedup after resolve).
→ `graph.bin` **v2**: `K, color[K], ncov[K], tcov[K], E, IEdge[E], total_payload, seqlen[K], payload[]` (window bytes inline). Pass-2 abort if any kept sig's slot is empty.

**rapidgzip cross-pass cache** (`CoreCache`, `SVP_NOCACHE` disables): pass 1 records each read's filtered core boundaries (`uint16 start,width`, fwd+rc); pass 2 replays them instead of re-running `extract_fast` (drops the fill pass's LCP cost to ~0). Keyed by `(shard, read-index)` — deterministic because the same `.gz` → same rapidgzip bytes → same parse order. Reads > 65535 bp → flagged `reext`, pass 2 falls back to `extract_fast`.

### `simplify <dir>` → `simplify.bin`  (`run_simplify`; compaction via `build_super_graph`)
The single graph-simplification stage — two phases, the second the namesake:
1. **Compaction — in RAM** (`build_super_graph`, `compact.cpp`): reads `graph.bin`, walks each maximal same-colour in/out-degree-1 chain, splices member windows by KMP overlap, RLE per-base coverage (`max` on overlaps); pure cycles → singletons. Never written to disk; consumed directly. (Log tag `[compact]`.)
2. **Bubble zip** (`simplify.cpp`): enumerate `(S1,S2)` colour-3 anchor bubbles (budgeted DFS, `NODE_BUDGET=500`), bucketed std_both/T_only/N_only. **v1 zips only `std_both`**; direct/cycle/multi-path pass through. For each zip: BFS the c1 + c2 arms, absorb the prefix/suffix shared across all of an anchor's arms into the **grown** anchor (only when "private" — every neighbour on that side is a consumed arm), leave each arm's divergent middle as a residual node, collapse equal alleles to a direct edge. Emits the **whole** simplified graph. (Log tag `[simplify]`.)
→ `simplify.bin`: `S + SuperNode[S] + blocks + runs(CovRun) + GEdge{from,to,overlap}` (drop `overlap` head bytes of `to` to reconstruct; bubble-seam edges have overlap 0, pass-through edges recompute KMP overlap on output bytes). Kept lean by CSR adjacency, **RLE arm coverage** (not flat per-base), generation-stamped BFS/DFS scratch. `SuperNode` has `pad[3]={}` so the file is byte-reproducible.

> "Bubble" here is the graph-topology object the zip operates on (an anchor pair with two allele arms) — the heart of the stage. It is unrelated to the deleted `bubbles.tsv` diagnostic.

## Build & run
```
make install-deps     # lib/lcptools (biolcp) + Badread venv + rapidgzip venv
make                  # ./svp from src/*.cpp + btllib bloom objects
./svp run <n.fq.gz> <t.fq.gz> <dir> [--min-mult 10 --threads 16]
```
Run `svp` from the repo root (the `make run` target does) so the default rapidgzip path resolves. **`-DLARGE_GENOME` is load-bearing — it must match `lib/lcptools`'s build** (64-bit core start/end) or structs corrupt. Defaults: `--level 5 --flank 10 --margin 10 --min-mult 10 --threads 16`. Keep `margin ≥ flank`. No prune — build at the threshold you want. svp no longer links zlib (`-lz` dropped); all decode is rapidgzip.

**rapidgzip** is a runtime dependency installed by `make install-deps` into `lib/rapidgzip/venv` (self-contained venv, repo-local — survives reboots, unlike the old `/tmp` location). The `graph` stage `popen`s it; default `lib/rapidgzip/venv/bin/rapidgzip`, override with `SVP_RAPIDGZIP`. Other tunables: `SVP_DECP` (rapidgzip threads/file), `SVP_FILL_NW` (fill worker count), `SVP_NOCACHE` (disable the cross-pass cache).

## Libraries
- **rapidgzip** (`lib/rapidgzip/venv`, v0.16) — parallel single-stream gzip inflate (~1.5 GB/s @ ~4 cores vs zlib's ~120 MB/s/core). The gz decoder for `graph`; invoked at runtime via `popen` (`open_decode`, `graph.cpp`), **not** linked.
- **lcptools** (`lib/lcptools`, branch `biolcp`, `-DLARGE_GENOME`) — LCP cores; `extract_fast.cpp` adapts `init_lps_offset`/`lps_deepen`.
- **btllib** — only the plain atomic `BloomFilter` headers (svp re-implements a non-atomic `PlainBloom` for the count cascade; the cascade runs under a per-shard lock so atomics were redundant). **MurmurHash3** (vendored) — the signature hash.
- Sim/eval (in `lib/`): **insilicoSV** (truth genome+VCF), **Badread** (HiFi read sim), **Severus** (intended baseline, not wired).
- `tools/` (git-tracked, **not** built by `make`; compiled out-of-band): `lcp_bench` (in-RAM thread-scaling micro-benchmark of the LCP hot path) and `lcp_stream` (LCP-only CPU cost over the real `.gz`, zlib decode). Both **VERBATIM-copy** `graph.cpp`'s `process_read` (it lives in an anon namespace) — they silently drift if the production parser changes. (`.gitignore` still whitelists a long-gone `tools/bubble_recall.cpp`.)

## Source layout (`src/`)
`graph.cpp` (the `graph` stage + readers/pools/bloom/cache/IO helpers), `compact.cpp` (`build_super_graph` — unitig compaction, phase 1 of `simplify`), `simplify.cpp` (`run_simplify` — bubble zip, phase 2), `extract_fast.cpp` (biolcp SoA adapter), `main.cpp` (CLI), `MurmurHash3.*`, `svp.h` (shared PODs + IO + `suffix_prefix`), `graph.h` (stage decls). One binary, one Makefile, no CMake/CI.

## Benchmarks
`ref.fasta → insilicoSV (truth out.fa + out.vcf) → Badread (normal from ref, tumor from out.fa) → svp`. Active sim `hard` (WG ~30×). `scripts/sim_reads_badread.sh`, `scripts/gen_badread_reads_wg.sh`, `simulation/Makefile`. `colo829/` is real tumor/normal HiFi (~107 GB / ~112 GB) + a lifted somatic-SV truthset; its last `run.log` was **OOM-killed** after the normal pool — no completed colo829 output yet (and that runner passed a `--ref` flag the CLI never accepted). See `BENCHMARK.md` for the full `graph`-stage record.

## Performance (measured 2026-06-17, `hard` WG, mm=10, 16 threads pinned `taskset -c 0-15`)
**`graph` = 18.3 min** (count 512 s + fill 502 s + serial tail ~82 s), 73.9% of the 16-core budget, peak RSS ~52 GB — vs 23.4 min for committed `eb0194e` (zlib + atomic bloom + re-LCP), output-equivalent (~86.4 M survivors, 88 M edges, 14.9 GB payload; `graph.bin` ≈ 16.5 GB, `simplify.bin` ≈ 7 GB). The two irreducible costs are LCP `extract_fast` (~3.66k core-s, count pass only now) + gzip decode (~1.67k core-s); the remaining lever is the fill-pass feed coupling and the `idx.find` on the 86 M-entry survivor map. `simplify` ≈ a few min after `graph`. **ALWAYS `--threads 16` + `taskset -c 0-15`** (see `BENCHMARK.md` — `--threads 32` breaks the budget and oversubscribes the readers + rapidgzip children).

## Code style
`docs/CODE_STYLE.md`: C++17, no OOP / no user templates / no try-catch (fail-fast `[FATAL]` + `abort`), avoid if-else trees, ≤30 LOC/function, comments only when the WHY is non-obvious. 0-based half-open coords (`core.end` exclusive). One binary, one Makefile, no CMake/CI.

## State note (working tree, not committed)
On top of `main` @ `eb0194e`: rapidgzip parallel decode + cross-pass `CoreCache` + plain bloom, **and** this cleanup — removed the orphaned `cores.cpp`/`io.cpp`, the legacy `lcp`/`graph <dir>` (`run_lcp`/`run_graph`/`run_pool`/`merge_sorted_nodes`/`pack_payload`, `counts.bin`/`edges.bin`), and the `bubbles` diagnostic (`run_bubbles`/`decompose`/`bubbles.tsv`); renamed `modea.cpp`→`compact.cpp`, `modeb.cpp`→`simplify.cpp`, command `modeb`→`simplify`, `modeb.bin`→`simplify.bin`, "Mode A"→compaction, "Mode B"→bubble-zip; moved rapidgzip out of `/tmp`. Build is output-equivalent to `eb0194e`. **`docs/PIPELINE.md` was rewritten to match; `BENCHMARK.md` `graph`-stage numbers are current.** Stable validation reference: svp `b7e0663` + lcptools `perf-core-extraction` (build the lib manually, no `-DLARGE_GENOME`). Don't port files across branches without instruction.
