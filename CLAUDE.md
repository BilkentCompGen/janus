# CLAUDE.md

_Project-specific guidance. Global behavioral rules + cross-project map: `/home/claude/CLAUDE.md`._

---

`svp` — reference-free somatic SV caller. Builds a 2-colour LCP-core graph straight from reads (normal vs tumor HiFi FASTQ), error-corrects it, and calls somatic variants off the coordinate-free core graph (per-node sequence + per-base coverage, no genome positions). C++17, one binary. No assembly, no alignment, no reference. Flow: **`run` = `graph → refine → call`** (three stages). Canonical `min_mult=3`. `graph` decodes `.gz` with **rapidgzip** (parallel, in-process) and skips pass-2 LCP via a **cross-pass core-boundary cache**.

## Hard rules
- **Reads only** (two FASTQ). An assembly/`.fa` is a ground-truth oracle, never a caller input.
- **Whole-genome only**; sim scope comes from the config, never hard-coded.

## Foundations
- **LCP core** (lcptools `biolcp`, level 5): a locally-consistent parse — the same locus yields the same core in any context. Adjacent L5 cores overlap ~85 bp, so nodes splice by exact KMP suffix/prefix (**zero overlap = FATAL**).
- **Signature** = `core_sig` = `XXH3-64(core bytes)` (no flank); the stored payload IS those core bytes. The byte-splice in `refine`/`call` relies on the native ~85 bp adjacent-core overlap, not a flank. (Earlier murmur+flank design + its measurements: `docs/FOUNDATIONS.md`.)
- **2 colours**: `(n && t) ? 3 : n ? 1 : 2` → **1 normal-only, 2 tumor-only (somatic), 3 shared**.
- **Dual fwd+RC** extraction per read (the parse isn't RC-symmetric): pools both-strand coverage (so `min_mult` prunes errors) but doubles node count and hides inversions as somatic. The two strands produce **independent** sigs (no per-emission canonicalization — that would break the KMP byte-splice). See `docs/FOUNDATIONS.md`.

## Pipeline
`run` = `graph → refine → call`. All stages resumable — `run` skips a stage whose output (`graph.bin` / `refined.bin` / `calls.tsv`) already exists; every checkpoint is `.tmp` + `rename`. (Simplify was removed from the pipeline; `docs/pipeline.html` still diagrams the old 2-stage flow and is stale.)

### `graph <normal.fq.gz> <tumor.fq.gz> <dir>` → `graph.bin`  (`run_graph_reads`)
Two passes over the same FASTQs, self-contained:
1. **Count** (`run_pool_both`): both samples in **one** pool. One pipe-reader per `.gz` reads plain FASTQ from a **rapidgzip** decode process (`open_decode` → `popen`); workers LCP each read (fwd + RC) → `emit_sig` → count in a sharded **bloom cascade**, one cascade per colour (`min_mult-1` plain bloom filters + exact survivor map). Per-colour `min_mult` gate happens **inside** the cascade, before the n/t merge (a normal-error sig drops, leaving the locus tumor-only). Bloom only over-keeps (FN = 0). Survivors sorted by sig → node id = sorted index (raw sigs are **not** persisted).
2. **Fill** (`fill_pool`): re-decode the same FASTQs; for each survivor write its core bytes (atomic single-writer claim) and emit adjacent-signature edges (per-worker deduped, periodic radix-compact, final dedup after resolve).
→ `graph.bin` **v2**: `K, color[K], ncov[K], tcov[K], E, IEdge[E], total_payload, seqlen[K], payload[]` (core bytes inline). Pass-2 abort if any kept sig's slot is empty.

**rapidgzip cross-pass cache** (`CoreCache`, `SVP_NOCACHE` disables): pass 1 records each read's filtered core boundaries (`uint16 start,width`, fwd+rc); pass 2 replays them instead of re-running `extract` (drops the fill pass's LCP cost to ~0). Keyed by `(shard, read-index)` — deterministic because the same `.gz` → same rapidgzip bytes → same parse order. Reads > 65535 bp → flagged `reext`, pass 2 falls back to `extract`.

### `refine <dir>` → `refined.bin`  (`run_refine`, refinement.cpp)
Topology + counts error-correction over the mm-floored core graph; **same v2 format** as `graph.bin` (so `call` reads it transparently). Each structure judged by its read count RELATIVE to local coverage (never a global gate). Passes: (1) **tips** — short dead-end chains far below their junction + low-count isolated islands; (2) **low-cov side-branches** — leaf sprouts far below their strongest neighbour; (3) **bubble pop** — the low, sub-het, non-somatic minority arm (germline-het + low-VAF somatic arms kept, flank-relative); then a tip re-pass. `abs_keep` floor (default 10) protects well-covered cores. Env: `REF_TIP_MAX/REF_TIP_FRAC/REF_ISLAND_CNT/REF_LOWCOV_FRAC/REF_BUB_MAX/REF_HET_FRAC/REF_SOM_FRAC/REF_ABS_KEEP/REF_ROUNDS`; `REF_TIP_NUKE=1` → aggressive tip removal to a distinct `refined_nuketip.bin` (never clobbers `refined.bin`). (Log `[refine]`.)

### `call <dir>` → `calls.tsv` + `calls.altsigs.bin`  (`run_call`, caller.cpp)
Reference-free somatic calls straight off `refined.bin` (colour-2 = tumor-only somatic). ONLY graph-intrinsic signals (colour, ncov/tcov, topology, payload). **BUBBLE**: a two-arm bubble with ≥1 tumor arm (`T/N`, `T/S`, `T/T`) → SNP/INS/DEL typed+sized by the **signed spliced arm-length Δ** (`|Δ|≤2`→SNP, `Δ>0`→INS, `Δ<0`→DEL). **BND**: a tumor-only run **bridging two shared (real) anchor cores at ~full depth** = a novel adjacency with no reference counterpart → DUP/INV/large-indel/transloc breakend (traced via the highest-depth tumor path to the nearest shared core each side; works inside the repeat tangle). **TIP**: tumor-only dead-end → insertion-end/breakend candidate. fwd+RC dups collapsed by a canonical (min seq/revcomp) haplotype hash. Env: `CALL_MINSUP` (default 3), `CALL_JMAX`, `CALL_DEPTH_FRAC`. (Log `[caller]`.)

## Build & run
```
make install-deps     # lib/lcptools (biolcp) + Badread venv + rapidgzip venv
make                  # ./svp from src/*.cpp + btllib bloom objects
./svp run <n.fq.gz> <t.fq.gz> <dir> [--min-mult 3 --threads 16]
```
Run `svp` from the repo root (the `make run` target does) so the default rapidgzip path resolves. **`-DLARGE_GENOME` is load-bearing — it must match `lib/lcptools`'s build** (64-bit core start/end) or structs corrupt. Defaults: `--level 5 --margin 10 --min-mult 3 --threads 16`. svp no longer links zlib (`-lz` dropped); all decode is rapidgzip.

**rapidgzip** is a runtime dependency installed by `make install-deps` into `lib/rapidgzip/venv` (self-contained venv, repo-local — survives reboots, unlike the old `/tmp` location). The `graph` stage `popen`s it; default `lib/rapidgzip/venv/bin/rapidgzip`, override with `SVP_RAPIDGZIP`. Other tunables: `SVP_DECP` (rapidgzip threads/file), `SVP_FILL_NW` (fill worker count), `SVP_NOCACHE` (disable the cross-pass cache).

## Libraries
- **rapidgzip** (`lib/rapidgzip/venv`, v0.16) — parallel single-stream gzip inflate (~1.5 GB/s @ ~4 cores vs zlib's ~120 MB/s/core). The gz decoder for `graph`; invoked at runtime via `popen` (`open_decode`, `graph.cpp`), **not** linked.
- **lcptools** (`lib/lcptools`, branch `biolcp`, `-DLARGE_GENOME`) — LCP cores; `extract.cpp` adapts `init_lps_offset`/`lps_deepen`.
- **btllib** — only the plain atomic `BloomFilter` headers (svp re-implements a non-atomic `PlainBloom` for the count cascade; the cascade runs under a per-shard lock so atomics were redundant). **xxHash** (`/usr/include/xxhash.h`, `XXH_INLINE_ALL`) — the `core_sig` identity hash (in `extract.cpp` + `caller.cpp`).
- Sim/eval (in `lib/`): **insilicoSV** (truth genome+VCF), **Badread** (HiFi read sim), **Severus** (intended baseline, not wired).
- `tools/` (git-tracked, **not** built by `make`; compiled out-of-band): `lcp_bench` (in-RAM thread-scaling micro-benchmark of the LCP hot path) and `lcp_stream` (LCP-only CPU cost over the real `.gz`, zlib decode). Both **VERBATIM-copy** `graph.cpp`'s `process_read` (it lives in an anon namespace) — they silently drift if the production parser changes. (`.gitignore` still whitelists a long-gone `tools/bubble_recall.cpp`.)

## Source layout (`src/`)
`graph.cpp` (the `graph` stage + readers/pools/bloom/cache/IO helpers), `extract.cpp` (biolcp SoA adapter + the XXH3 `core_sig`), `refinement.cpp` (`run_refine`), `caller.cpp` (`run_call`), `main.cpp` (CLI), `svp.h` (shared PODs + IO helpers), `graph.h` (stage decls). One binary, one Makefile, no CMake/CI.

## Benchmarks
`ref.fasta → insilicoSV (truth out.fa + out.vcf) → Badread (normal from ref, tumor from out.fa) → svp`. Active sim `hard` (WG ~30×). `scripts/sim_reads_badread.sh`, `scripts/gen_badread_reads_wg.sh`, `simulation/Makefile`. `colo829/` is real tumor/normal HiFi (~107 GB / ~112 GB) + a lifted somatic-SV truthset; its last `run.log` was **OOM-killed** after the normal pool — no completed colo829 output yet (and that runner passed a `--ref` flag the CLI never accepted). See `BENCHMARK.md` for the full `graph`-stage record.

## Performance (measured 2026-06-17, `hard` WG, mm=10, 16 threads pinned `taskset -c 0-15`)
**`graph` = 18.3 min** (count 512 s + fill 502 s + serial tail ~82 s), 73.9% of the 16-core budget, peak RSS ~52 GB — vs 23.4 min for committed `eb0194e` (zlib + atomic bloom + re-LCP), output-equivalent (~86.4 M survivors, 88 M edges, 14.9 GB payload; `graph.bin` ≈ 16.5 GB). The two irreducible costs are LCP extraction (~3.66k core-s, count pass only now) + gzip decode (~1.67k core-s); the remaining lever is the fill-pass feed coupling and the `idx.find` on the 86 M-entry survivor map. `refine` + `call` ≈ a few min each after `graph`. **ALWAYS `--threads 16` + `taskset -c 0-15`** (see `BENCHMARK.md` — `--threads 32` breaks the budget and oversubscribes the readers + rapidgzip children).

## Code style
`docs/CODE_STYLE.md`: C++17, no OOP / no user templates / no try-catch (fail-fast `[FATAL]` + `abort`), avoid if-else trees, ≤30 LOC/function, comments only when the WHY is non-obvious. 0-based half-open coords (`core.end` exclusive). One binary, one Makefile, no CMake/CI.

## State note (branch `changes`)
Pipeline = `graph → refine → call` (committed `ac21697`): caller added as `src/caller.cpp` (the `call` stage); the now-skipped simplify stage (`simplify.cpp`, `resolve.cpp` + their hooks) and dead `MurmurHash3` were removed. `graph.cpp`/`refinement.cpp` algorithms are unchanged by that cleanup; `refined.bin` + `calls` are byte-reproducible given the same `graph.bin`. **Caveat:** `graph.bin` itself is NOT byte-reproducible run-to-run — the parallel bloom cascade's false-positive over-keep is order-dependent (~250 junk error cores / 0.0003% jitter); the **real-core set + caller metrics are deterministic** (`hard` WG sim, mm=3: recall 95.7%, precision 84.2%; ground-truth tooling lives in scratch — see the `svp-error-correction` memory). `docs/pipeline.html` still diagrams the old 2-stage flow and is stale (not yet regenerated). The opt-in `SVP_SPLIT_REPEATS` repeat-capture in `graph.cpp` is left dormant. Don't port files across branches without instruction.
