# CLAUDE.md

_Project-specific guidance. Global behavioral rules + cross-project map: `/home/claude/CLAUDE.md`._

---

`svp` — reference-free somatic SV caller. Builds a 2-colour LCP-core graph straight from reads (normal vs tumor HiFi FASTQ), error-corrects it, and calls somatic variants off the coordinate-free core graph (no genome positions, and since **v3** no stored sequence either — a node IS its 64-bit core signature + its per-colour read counts). C++17, one binary. No assembly, no alignment, no reference in the caller. Flow: **`svp run` = `graph → refine`** (two stages); the call + anchor stages are out-of-band tools (`tools/call_v3`, `tools/anchor_v3`). Canonical `min_mult=3`. `graph` decodes `.gz` with **rapidgzip** (parallel, in-process).

> **DIRECTION CHANGE, 2026-07-28.** The reference is leaving the caller for good. `refidx` (strand vote +
> repeat blacklist/merge) is **deleted** — both strands go in, repeats are ordinary nodes. Next: an output
> vocabulary of **DEL / INS / DUP / INV / CPX / SMV** (complex events are NOT decomposed — CPX is the answer,
> not a fallback), calls in **normal-graph coordinates** instead of reference positions, our own scorer (which
> also absorbs the mirror-twin duplicate-call problem), and a way back from a call to its reads so resolution
> below the ~152 bp L5 core can be bought where it is actually needed. The reference-anchored baseline this
> replaces — how Truvari scores and what the old pipeline measured — is recorded in `docs/TRUVARI.md`.
> Everything below this line still describes the graph stage as built; the caller sections describe the OLD path.

## Hard rules
- **Reads only** (two FASTQ / `.svb` sample dirs). An assembly/`.fa` is a ground-truth oracle or an anchoring reference — never a caller input.
- **Whole-genome only**; sim scope comes from the config, never hard-coded.
- **The v2 payload format is retired. Do NOT restore `src/caller.cpp` / `src/anchor.cpp`** (deleted in the working tree, preserved in backup commit `2a59122`). They reached ~85% accuracy (per user) and remain the reference *spec*, but they require the v2 on-disk payload (`color` + `seqlen` + raw core bytes) that no longer exists. `tools/call_v3.cpp` + `tools/anchor_v3.cpp` are their faithful sig-only ports — extend those, don't re-add the originals to `src/`.
- Don't port files across branches without instruction.

## Foundations
- **LCP core** (lcptools `biolcp`, level 5): a locally-consistent parse — the same locus yields the same core in any context. Adjacent L5 cores overlap heavily (~85–105 bp measured), so a graph edge is a genuine sequence overlap, and the *core step* (bp per core) is far below the core length — this is what lets `anchor_v3` convert core counts → bp.
- **Signature** = `core_sig` = **MurmurHash64A over the core's own bytes** (seed 0, no flank), computed once in `extract()` (`src/extract.cpp`). The sig **is** the node identity and the only per-node payload stored on disk. (Earlier murmur+flank / XXH3-of-payload designs and their measurements: `docs/FOUNDATIONS.md`.)
- **2 colours, DERIVED from counts** (never stored): `(ncov && tcov) ? 3 : ncov ? 1 : 2` → **1 normal-only, 2 tumor-only (somatic), 3 shared**.
- **NO STRAND RESOLUTION — both strands go in** (2026-07-28: the reference index was deleted). The parse isn't RC-symmetric, so a read's fwd and RC bytes give two *independent* sig families with different core boundaries; nothing reference-free can say which one a read was sequenced in, so **both are inserted**. Every locus therefore has a **mirror twin** in the other family: two nodes, different sigs, **no core in common**. Pairing the twins is the scorer's job — a somatic event called once per family is ONE event, not two. (The retired reference vote halved the graph and is preserved in git history; the project's best measured result, hard F1 0.737, came from a dual-strand graph anyway.)
- **A large inversion's INTERIOR is shared, not somatic** — the inverted sequence's reads land on the normal reads' mirror-family nodes, so the interior looks like normal material. What is callable is the two BREAKPOINT junctions.

## Graph format (v3, sig-only)
`graph.bin` and `refined.bin` share ONE layout — `refine` and the v3 tools read either transparently:
```
u64 K | u64 sig[K] | u32 ncov[K] | u32 tcov[K] | u64 E | CEdge{u32 src,dst,ncov,tcov}[E]
```
- **No `color`, no `seqlen`, no `payload`** (v2 carried all three; v3 carries none). No sequence is stored anywhere, so all sizes downstream are **CORE COUNTS**, not bp — only `anchor_v3` turns them into bp.
- **Edges carry per-colour counts** (new in v3): `ncov`/`tcov` = how many normal/tumor reads cross that adjacency. Used by refinement's edge gate + peel loop and by `call_v3`'s depth logic.
- Node id = index into the sig-sorted survivor array; `sig[]` is a bijection with the id, so `sig2node` reconstructs it anywhere.
- 16 B/node + 16 B/edge: e.g. `outputs/svp/graph_mm1` (K=95.6 M, E=98.9 M) = 3.1 GB, vs ~16.5 GB for the v2 payload-era graph of comparable size.

## Pipeline
`svp run` = `graph → refine`, both resumable — `run` skips a stage whose output (`graph.bin` / `refined.bin`) already exists; every checkpoint is `.tmp` + `rename`. Calling is out-of-band: `tools/call_v3` → `tools/anchor_v3` → `calls.vcf` → Truvari (`scripts/truvari.sh`).

### `graph <normal> <tumor> <dir>` → `graph.bin`  (`run_graph_reads`, graph.cpp)
Each sample arg is a single file **or a directory of shards** (`.fastq` / `.fq` / `.svb`, read in parallel). Two passes over the same inputs, self-contained:
1. **Count** (`run_pool_both`): both samples in **one** pool. One pipe-reader per `.gz` reads plain FASTQ from a **rapidgzip** decode process (`open_decode` → `popen`); `.svb` shards are `fopen`ed and 2-bit-decoded directly (`svb_next`). Workers LCP each read (fwd + RC), drop cores within `margin` of a read end, and count each sig in a sharded **bloom cascade**, one cascade per colour (`min_mult-1` plain bloom filters + exact survivor map). Per-colour `min_mult` gate happens **inside** the cascade, before the n/t merge (a normal-error sig drops, leaving the locus tumor-only). Bloom only over-keeps (FN = 0). Survivors merged into one `(sig,n,t)` record, sorted by sig → node id = sorted index.
2. **Edge count** (`edge_count_pool`): re-decode + re-LCP the same reads; for every consecutive pair of *surviving* cores in a read, increment the `(prev,cur)` adjacency in the reader's colour (per-worker map → 64 sig-sharded global maps). A core that didn't survive breaks the chain. Sig-pair edges are then resolved to node ids → `CEdge[E]` (self-loops dropped, sorted by `(src,dst)`).

**No reference index, no blacklist, no repeat merge.** A read's kept cores enter as `extract()` produced them, both strands. Repeat cores are ordinary nodes again (high-degree hubs are back); the only repeat handling left is coverage-based (`SVP_SPLIT_REPEATS`, below) — a run collapsed from C loci carries ~C× the single-copy depth. Measured cost of the removal on `complex` (same reads, refidx on → off): nodes 48,311,744 → **99,287,043** (2.055×), edges 50,320,569 → **108,268,381**, `graph.bin` 1.58 → 3.32 GB, wall 752.9 s → **849.8 s (only 1.13×** — the bottleneck is LCP + gz decode, both already paid twice).
Index: 86.4M cores, 363 MB, 4.19 B/entry (22-bit bucket + 29-bit fingerprint + 3 flag bits). FW-only 42,926,357 · RC-only 42,926,620 · BOTH 589,004 (0.68%) · REPEAT 1,552,839 (1.80%). Verified against the two hash sets it replaced: **0 strand disagreements / 87M sigs, 0 missed + 0 invented repeats**, FP rate 2e-8 (≈9 spurious votes genome-wide). Memory 1.56 GB → 0.36 GB, load 3.9 s → 0.5 ms; lookup 127 → 135 ns (a wash).

**`SVP_COUNT_ONLY=1`** → dump every observed survivor as `counts.bin` (`u64 K` + `SigCount{sig,n,t}[K]`) and stop before the edge pass (analysis mode: what would any mm-floor drop?).

**rapidgzip cross-pass cache** (`CoreCache`): pass 1 records each read's filtered core boundaries (`uint16 start,width`, fwd+rc, keyed by `(shard, read-index)`), now ON by default (`SVP_NOCACHE=1` disables). **Still vestigial** — it is populated but *never consumed*: pass 2 (`edge_count_pool`) re-extracts unconditionally. Rewire or drop it before trusting the old "fill-pass LCP ≈ 0" claim.

### `refine <dir>` → `refined.bin`  (`run_refine`, refinement.cpp)
Topology + counts error-correction over the mm-floored graph; **same v3 format** in and out. Every structure is judged by its read count RELATIVE to local coverage (never a global gate). Default path:
1. **tips** — short dead-end chains far below their junction + low-count isolated islands (iterated, `REF_ROUNDS`).
2. **low-cov side-branches** — leaf sprouts far below their strongest neighbour.
3. **bubble pop** — the low, sub-het, non-somatic minority arm (germline-het + low-VAF somatic arms kept, judged **flank-relative** per colour against the bubble's anchor coverage); then a tip re-pass.
4. **edge gate** (`pass_edges`, v3): edges carry raw traversal counts (min 1) while node counts are mm-floored, so a "tumor-only" edge can rest on a single tumor read — the dominant somatic-junction FP. HARD floor `max(ncov,tcov) >= REF_EDGE_MM` (default 3) + RELATIVE gate (drop an edge carrying < `REF_EDGE_FRAC`, default 0.10, of its endpoints' coverage in *every* colour it has).

`abs_keep` floor (default 10) protects well-covered cores in tips/low-cov. Env: `REF_TIP_MAX/REF_TIP_FRAC/REF_ISLAND_CNT/REF_LOWCOV_FRAC/REF_BUB_MAX/REF_HET_FRAC/REF_SOM_FRAC/REF_ABS_KEEP/REF_ROUNDS/REF_EDGE_MM/REF_EDGE_FRAC`. (Log `[refine]`.)

**Opt-in variants, each to a DISTINCT file (never clobbers `refined.bin`):**
- `REF_ITER=<N>` → **iterative node↔edge peel** to `refined_iter.bin`: one-shot flow-conservation filter (`REF_FLOW=1`: drop nodes whose summed incident edge counts exceed the node's own count — ~95% error cores), then up to N rounds of *gate edges → subtract each dropped edge's counts from its endpoint nodes → remove zero-coverage nodes → re-run the node passes*, stopping at convergence. `REF_EDGE_CLAMP=1` (clamp edges to `min(endpoint counts)`) is OFF by default — measured to recolour shared edges to tumor-only (+~55k false somatic edges).
- `REF_TIP_NUKE=1` → aggressive tip removal to `refined_nuketip.bin`.

### `tools/call_v3 <graph_v3.bin> <out_prefix> [minsup=3]` → `<prefix>.tsv` + `<prefix>.altsigs.bin`
Reference-free somatic calls straight off `refined.bin` (colour-2 = tumor-only somatic). ONLY graph-intrinsic signals (derived colour, ncov/tcov, topology). **BUBBLE**: a two-arm bubble with ≥1 tumor arm (`T/N`, `T/S`, `T/T`) → SNP/INS/DEL typed+sized by the **signed arm CORE-COUNT Δ** (`|Δ| ≤ SNP_EPS=0` → SNP, `Δ>0` → INS, `Δ<0` → DEL). **BND**: a tumor-only run **bridging two shared (real) anchor cores at ~full depth** (`CALL_DEPTH_FRAC`, default 0.30) = a novel adjacency with no reference counterpart → DUP/INV/large-indel/transloc breakend, traced via the highest-tcov tumor path to the nearest shared core each side (works inside the repeat tangle); flagged `BND_CNgain` if an anchor's T/N depth ratio ≥ 1.5. **TIP**: tumor-only dead-end (≤16 cores) → insertion-end/breakend candidate. Dedup = hash of the alt core-sig SET. `altsigs.bin` = per call `{u64 id, u32 nsig, u64 sig[]}`. Env: `CALL_MINSUP` (3), `CALL_JMAX` (12), `CALL_DEPTH_FRAC`. (Log `[call]`.)

### `tools/anchor_v3 <graph_v3.bin> <ref.fasta> <calls_prefix> <out.vcf> [threads=16]` → VCF
The calls are coordinate-free; anchoring maps them onto a reference so Truvari can score them. A somatic alt is tumor-only (absent from the reference), but its shared/normal-only graph NEIGHBOURS are in the reference — so: collect each call's anchor nodes, **contig-parallel re-extract the reference** (same level/margin/repeat merge; needs `ref.fasta.fai`), match by sig, drop sigs that land at ≥2 distinct positions (repeat → ambiguous), and bracket the call between its located anchors on the consensus contig. Core counts → bp via the **average core step measured on that scan** (`tot_bp / tot_cores`). Re-types/re-sizes BND↔DEL/INS from the anchored reference span, dedups same-class calls at the same locus (keeps highest support), emits `SVTYPE/SVLEN/END/NANCHOR/CALLSRC/SUBTYPE` + `LOCALSUB`/`TIP` filters. (Log `[anchor]`.)

### Scoring
`scripts/truvari.sh <truth.vcf> <calls.vcf> <ref.fa> [outdir]` (needs `.venv`; runs `truvari bench --pctseq 0`, prints TP/FP/FN + a per-type FN/FP breakdown), or `cd simulation && make bench/hard CALL=<vcf>`.

## Build & run
```
make install-deps     # lib/lcptools (biolcp) + Badread venv + rapidgzip venv
make                  # ./svp from src/*.cpp + btllib bloom objects
./svp run <n.fq.gz> <t.fq.gz> <dir> [--min-mult 3 --threads 16]
```
Run `svp` from the repo root (the `make run` target does) so the default rapidgzip path resolves. **`-DLARGE_GENOME` is load-bearing — it must match `lib/lcptools`'s build** (64-bit core start/end) or structs corrupt. CLI defaults: `--level 5 --margin 10 --min-mult 3 --threads 16` (the `make run` target passes `MINMULT=1`). svp does not link zlib; all decode is rapidgzip. PGO is wired in the Makefile (`make PGO=generate` → run → `make PGO=use`).

**rapidgzip** is a runtime dependency installed by `make install-deps` into `lib/rapidgzip/venv` (self-contained venv, repo-local — survives reboots). The `graph` stage `popen`s it; override with `SVP_RAPIDGZIP`. Other tunables: `SVP_DECP` (rapidgzip threads/file), `SVP_FILL_NW` (edge-pass worker count), `SVP_NOCACHE`.

**The `tools/*_v3` binaries are NOT built by `make`** and no build script for them is tracked — compile out-of-band. `call_v3.cpp` is self-contained (no lcptools). `anchor_v3.cpp` includes `svp.h`/`unordered_dense.h` and needs `build/extract.o` + `liblcptools` (`-DLARGE_GENOME`) so its parse and sigs are byte-identical to the graph's.

## Libraries
- **rapidgzip** (`lib/rapidgzip/venv`, v0.16) — parallel single-stream gzip inflate (~1.5 GB/s @ ~4 cores vs zlib's ~120 MB/s/core). The gz decoder for `graph`; invoked at runtime via `popen` (`open_decode`, `graph.cpp`), **not** linked.
- **lcptools** (`lib/lcptools`, branch `biolcp`, `-DLARGE_GENOME`) — LCP cores; `extract.cpp` adapts `init_lps_offset`/`lps_deepen`.
- **btllib** — only the plain atomic `BloomFilter` headers (svp re-implements a non-atomic `PlainBloom` for the count cascade; the cascade runs under a per-shard lock so atomics were redundant). **ankerl::unordered_dense** (vendored `src/unordered_dense.h`) is the flat hash map everywhere. No xxHash any more — node identity is MurmurHash64A in `extract.cpp`.
- Sim/eval (in `lib/`): **insilicoSV** (truth genome+VCF), **Badread** (HiFi read sim), **Severus** (intended baseline, not wired).
- `tools/` sources are **all tracked now** (force-added 2026-07-28 in commit `70ec568`, before the direction change); binaries stay ignored. None of it is built by `make` — compile out-of-band.
  - v3-aware: `call_v3`, `anchor_v3`, `census_v3` (topology census of a v3 bin: degrees, SCCs, components, unitig N50), `fasta_cores` (a FASTA's effective-core sig set — the node-identity ground truth), `svb_pack`/`svb_dump` (FASTQ ↔ 2-bit `.svb`).
  - v2-only, **cannot run on a v3 bin**: `trace` (needs `color`/`seqlen`/payload), `lcp_bench`, `lcp_stream` (both VERBATIM-copy an old `graph.cpp` `process_read` and silently drift).
  - Active LCP probes (`lcp_symmetry`, `lcp_union`, `lcp_consistency`, `lcp_warmup`, `uniquify_levels`) — the RC-symmetry / union-cut / boundary-warm-up investigation behind the dual-strand question. See `docs/FOUNDATIONS.md`.

## Source layout (`src/`)
`graph.cpp` (the `graph` stage + readers/pools/bloom/`.svb`/`CoreCache`/IO helpers), `extract.cpp` (biolcp SoA adapter + the MurmurHash64A `core_sig`), `refinement.cpp` (`run_refine`), `posmap.{h,cpp}` + `placement.{h,cpp}` (the opt-in `SVP_POSAWARE` (sig,locus) node model — the ONLY reference-dependent code left, and it has never been run), `main.cpp` (CLI: `graph` / `refine` / `run` — **no `call`, no `anchor`**), `svp.h` (shared PODs + `Span` + IO helpers), `graph.h` (stage decls), `unordered_dense.h` (vendored). One binary, one Makefile, no CMake/CI.

## Benchmarks
`ref.fasta → insilicoSV (truth out.fa + out.vcf) → Badread (normal from ref, tumor from out.fa) → svp`. Baseline sim `hard` (WG ~30×, hypermutated: ~118k SNVs, ~6.6k small indels, ~600 SVs, `min_intersv_dist=1` so small variants pepper the junctions). Plus **5 complex sims layered on the same hard field**: `tra_intra` / `tra_inter` / `tra_inv` (intra/inter/inverted translocations), `compound` (2–3-SV structures, every interval ≥1 kb), `chrom_cnv` (chromosome-scale arm→whole-chrom DEL/DUP, 10–80 Mb). Each appends its class to the verbatim `hard` variant field and emits truth **novel-adjacencies (`adjacencies.bed`) + `sim.paf`** for a junction-level evaluator (no re-sim needed). Complex-sim tumor reads are compact 2-bit `.svb` shards (Badread → `tools/svb_pack`); the FASTQ **normal** reads (`reads/normal/`) are shared across all sims.
```
cd simulation && make reads-cplx INSILICOSV=/home/yuio/svp/.venv/bin/insilicosv THREADS=16   # resumable
```
`scripts/sim_reads_badread.sh`, `scripts/gen_badread_reads_wg.sh`, `simulation/Makefile`. See `BENCHMARK.md` for the full `graph`-stage record.
*(Carried over, UNVERIFIED: a `colo829/` real tumor/normal HiFi set (~107 GB / ~112 GB) + lifted somatic-SV truthset whose last `run.log` was OOM-killed after the normal pool. There is **no `colo829/` in the working tree** and nothing about it in `BENCHMARK.md` or `scripts/` — treat as historical until it reappears.)*

## Performance (measured 2026-06-17 — on the **v2** build; NOT re-measured for v3)
`hard` WG, mm=10, 16 threads pinned `taskset -c 0-15`: **`graph` = 18.3 min** (count 512 s + pass-2 fill 502 s + serial tail ~82 s), 73.9% of the 16-core budget, peak RSS ~52 GB — vs 23.4 min for committed `eb0194e`. The two irreducible costs are LCP extraction (~3.66k core-s) + gzip decode (~1.67k core-s). **ALWAYS `--threads 16` + `taskset -c 0-15`** (see `BENCHMARK.md` — `--threads 32` breaks the budget and oversubscribes the readers + rapidgzip children).
**What v3 changes:** pass 2 is now edge counting, not payload fill — it still re-decodes and re-LCPs both samples (the boundary cache is not wired into it), so its LCP cost is *back*, while the payload write and the 86 M-entry `idx.find`-per-core fill are gone; `graph.bin` shrinks ~5×. Treat the wall-clock above as the v2 record until a v3 run is timed.

## Code style
`docs/CODE_STYLE.md`: C++17, no OOP / no user templates / no try-catch (fail-fast `[FATAL]` + `abort`), avoid if-else trees, ≤30 LOC/function, comments only when the WHY is non-obvious. 0-based half-open coords (`core.end` exclusive). One binary, one Makefile, no CMake/CI.

## Caveats
- `graph.bin` is NOT byte-reproducible run-to-run (order-dependent bloom over-keep, ~0.0003% junk-core jitter); the real-core set is deterministic.
- Stale comments in code (v2 leftovers, harmless but misleading): `refinement.cpp`'s header still says "IDENTICAL v2 format"; `svp.h` still declares `RepTriple` and the `SVP_SPLIT_REPEATS` repeat-capture, but **nothing writes `threads.bin` any more** (only `tools/trace` still reads it).
- `.gitignore` whitelists a long-gone `tools/bubble_recall.cpp`.
- Uncommitted working tree on branch `changes` (backup commit `2a59122`). Do NOT commit without instruction.
