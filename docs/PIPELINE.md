# Pipeline — current stages and algorithms

`svp` is a reference-free somatic SV caller. It builds a 2-colour LCP-core graph directly from
reads (normal vs tumor HiFi FASTQ), error-corrects it, then calls somatic variants off it. No
assembly, no alignment, no reference. The graph is **coordinate-free**: per-node core sequence +
per-base coverage, not genome positions.

End-to-end benchmark flow:
`reference → insilicoSV (truth) → Badread (reads) → FASTQ → [svp: graph → refine → call]`.

## CLI — three resumable, checkpointed stages
The expensive LCP pass is paid once in `graph`; `run` skips any stage whose checkpoint already
exists. Checkpoints are atomic (`.tmp` + rename).

```
svp graph  <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5 --margin 10 --min-mult 3 --threads 16]
svp refine <dir>
svp call   <dir>
svp run    <normal.fq.gz> <tumor.fq.gz> <dir>   # graph → refine → call, skipping completed stages
```

## 1. graph  (`run_graph_reads`, graph.cpp)  — reads → self-contained `graph.bin`
Two passes over the same FASTQs.

**Decode.** One pipe-reader per `.gz` reads plain FASTQ from a **rapidgzip** decode process
(`open_decode` → `popen`; parallel single-stream inflate, ~1.5 GB/s @ ~4 cores vs zlib's
~120 MB/s/core). Tunables: `SVP_RAPIDGZIP` (binary path, default `lib/rapidgzip/venv/bin/rapidgzip`),
`SVP_DECP` (decode threads/file), `SVP_FILL_NW` (fill workers), `SVP_NOCACHE` (disable the cache).

**Pass 1 — count (`run_pool_both`).** Both samples in one pool. Workers run **two LCP extractions
per read, forward and reverse-complement** (the parse is not RC-symmetric). Each core's identity is
`core_sig` = **XXH3-64 over the core's own bytes** (no flank). Each sig is counted in a sharded
**bloom cascade** — one cascade per colour (`min_mult-1` plain bloom filters + an exact survivor
map). **Survivors** = sigs seen **≥ `min_mult`** times; they become the graph's nodes. The
per-colour `min_mult` gate happens *inside* the cascade, **before** the normal/tumor merge (a
normal-error sig drops, leaving its locus tumor-only). Bloom collisions can only over-keep (false
survivor); they never drop a real one (FN = 0). Survivors are sorted by sig, **node id = sorted
index** — the raw 64-bit sigs are not persisted. Cores within `margin` of a read end are dropped.

**Pass 2 — fill (`fill_pool`).** Re-decode the same FASTQs; for each survivor write its **core
bytes** once (atomic single-writer claim) and emit adjacent-signature edges (per-worker deduped,
periodic radix-compaction, final dedup after resolving sig→node-id). A **cross-pass core-boundary
cache** (`CoreCache`) recorded each read's filtered core `(start,width)` in pass 1 so pass 2 replays
them instead of re-extracting (deterministic because the same `.gz` → same rapidgzip bytes → same
parse order; reads > 65535 bp fall back to re-extraction).

**Colour.** `1 = normal-only, 2 = tumor-only (candidate somatic), 3 = shared.`
→ `graph.bin` **v2**: `u64 K · u8 color[K] · u32 ncov[K] · u32 tcov[K] · u64 E · IEdge{u32,u32}[E]
· u64 total_payload · u16 seqlen[K] · u8 payload[]` (the core bytes inline — same region as the sig,
no flank).

## 2. refine  (`run_refine`, refinement.cpp)  — `graph.bin` → `refined.bin`
Topology + counts error-correction over the mm-floored core graph; **same v2 format**. Each structure
is judged by its read count **relative to local coverage** (never a global gate — that is what an
over-aggressive `min_mult` does, distorting real bubbles). Passes:
1. **tips** — short dead-end core-chains whose count ≪ the junction they hang off, plus tiny
   low-count isolated contigs (pure error islands). Iterated (removal exposes new tips).
2. **low-cov side-branches** — leaf-ish sprouts whose count is a small fraction of their strongest
   neighbour path (recurrent-error cores clinging to a real path).
3. **bubble popping** — the low, error-like minority arm of a short bubble. Conservative: a balanced
   bubble is a germline het (keep both); a low tumor-only arm above a VAF floor is a candidate
   somatic (keep); only a clearly imbalanced, non-somatic minority arm is popped. Then a tip re-pass.

An `abs_keep` floor (default 10) protects well-covered cores in tips/low-cov. Env-tunable:
`REF_TIP_MAX REF_TIP_FRAC REF_ISLAND_CNT REF_LOWCOV_FRAC REF_BUB_MAX REF_HET_FRAC REF_SOM_FRAC
REF_ABS_KEEP REF_ROUNDS`. `REF_TIP_NUKE=1` = aggressive tip removal, written to a distinct
`dir/refined_nuketip.bin` (never clobbers `refined.bin`). (Log `[refine]`.)

## 3. call  (`run_call`, caller.cpp)  — `refined.bin` → `calls.tsv` + `calls.altsigs.bin`
Reference-free somatic variant calling straight off the refined core graph (no simplify). Uses
**only graph-intrinsic signals** — colour, per-node ncov/tcov, topology, payload. A somatic allele
always involves a tumor-only (colour-2) core. Call sources:
- **BUBBLE** — a two-arm bubble carrying ≥1 tumor arm → `T/N`, `T/S`, `T/T`. The tumor arm is the alt
  allele; the SNP/INS/DEL type and size come from the signed **spliced arm-length Δ**
  (`|Δ|≤2` → SNP, `Δ>0` → INS of size Δ, `Δ<0` → DEL of size −Δ).
- **BND** — a tumor-only run that **bridges two shared (real) anchor cores at ~full tumor depth** =
  a novel adjacency with no reference counterpart (so it never forms a bubble) → DUP / INV /
  large-indel / translocation breakend. Recovered by tracing the highest-depth tumor path to the
  nearest shared core on each side (works inside the repeat tangle, no degree-1 requirement).
- **TIP** — a tumor-only dead-end chain → insertion-end / breakend candidate.

fwd+RC strand duplicates are collapsed by a canonical (min seq/revcomp) haplotype hash. Env:
`CALL_MINSUP` (default 3), `CALL_JMAX`, `CALL_DEPTH_FRAC`. Outputs `calls.tsv` (per-call features +
sequences) and `calls.altsigs.bin` (per-call alt-allele core sigs, for recall/precision scoring).
(Log `[caller]`.)

## Source layout (`src/`)
`graph.cpp` (graph stage + readers/pools/bloom cascade/`CoreCache`/checkpoint IO), `extract.cpp`
(biolcp SoA adapter + the XXH3 `core_sig`), `refinement.cpp` (`run_refine`), `caller.cpp`
(`run_call`), `main.cpp` (CLI), `svp.h` (shared PODs + IO helpers), `graph.h` (stage decls).
One binary, one Makefile, no CMake/CI.

## Build flag
`-DLARGE_GENOME` is load-bearing and must match `lib/lcptools`'s build (64-bit core start/end on a
whole-genome reference). Rebuild the static lib if you see segfaults on struct fields. svp links no
zlib — all decode is rapidgzip (installed by `make install-deps` into `lib/rapidgzip/venv`).

## Known gaps / next steps
- **Short insertions (<50 bp)** are below LCP-core resolution (median INS ~25 bp) → low recall;
  ≥50 bp DEL/DUP/INV recover ~95–98%.
- **DUP vs INV sub-typing** and the **somatic-vs-error FP filter** (the T/S-SNP and tip false
  positives) are downstream of the current high-recall call set.
- Inversions/DUPs are limited by dual fwd+RC extraction (see FOUNDATIONS.md) — strand folding is the
  precondition for cleaner recovery.
