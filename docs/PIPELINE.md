# Pipeline — current stages and algorithms

`svp` is a reference-free somatic SV caller. It builds a 2-colour LCP-core graph directly from
reads (normal vs tumor HiFi FASTQ) and simplifies it. No assembly, no alignment. Output is
**coordinate-free**: per-node sequence + per-base coverage, not genome positions.

End-to-end benchmark flow:
`reference → insilicoSV (truth) → Badread (reads) → FASTQ → [svp: graph → simplify]`.

## CLI — two resumable, checkpointed stages
The expensive LCP pass is paid once; `run` skips any stage whose checkpoint already exists.
Checkpoints are atomic (`.tmp` + rename).

```
svp graph    <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5 --flank 10 --margin 10 --min-mult 10 --threads 16]
svp simplify <dir>
svp run      <normal.fq.gz> <tumor.fq.gz> <dir>   # runs both, skipping done stages
```

## 1. graph  (`run_graph_reads`, graph.cpp)  — reads → self-contained `graph.bin`
Two passes over the same FASTQs.

**Decode.** One pipe-reader per `.gz` reads plain FASTQ from a **rapidgzip** decode process
(`open_decode` → `popen`; parallel single-stream inflate, ~1.5 GB/s @ ~4 cores vs zlib's
~120 MB/s/core). Tunables: `SVP_RAPIDGZIP` (binary path, default `lib/rapidgzip/venv/bin/rapidgzip`),
`SVP_DECP` (decode threads/file), `SVP_FILL_NW` (fill workers), `SVP_NOCACHE` (disable the cache).

**Pass 1 — count (`run_pool_both`).** Both samples in one pool. Workers run **two LCP
extractions per read, forward and reverse-complement** (remap+reverse through a 256-byte
complement LUT). Each surviving core's signature = `murmur64(core ±flank, seed 0)`. Each sig is
counted in a sharded **bloom cascade** — one cascade per colour (`min_mult-1` plain bloom
filters + an exact survivor map). **Survivors** = signatures seen **≥ `min_mult`** times; they
become the graph's nodes. Counting fewer-than-`min_mult` sigs are dropped as likely errors. The
per-colour `min_mult` gate happens *inside* the cascade, **before** the normal/tumor merge (a
normal-error sig drops, leaving its locus tumor-only). Bloom collisions can only over-keep
(false survivor); they never drop a real one (FN = 0). Survivors are sorted by sig, and **node
id = sorted index** — the raw 64-bit sigs are not persisted.

**Pass 2 — fill (`fill_pool`).** Re-decode the same FASTQs; for each survivor fill its ±flank
window bytes once (atomic single-writer claim) and emit adjacent-signature edges (per-worker
deduped, periodic radix-compaction, final dedup after resolving sig→node-id). A **cross-pass
core-boundary cache** (`CoreCache`) recorded each read's filtered core `(start,width)` in pass 1
so pass 2 replays them instead of re-running `extract_fast` (deterministic because the same
`.gz` → same rapidgzip bytes → same parse order; reads > 65535 bp fall back to `extract_fast`).

**Colour.** `1 = normal-only, 2 = tumor-only (candidate somatic), 3 = shared.`
→ `graph.bin` **v2**: `u64 K · u8 color[K] · u32 ncov[K] · u32 tcov[K] · u64 E · IEdge{u32,u32}[E]
· u64 total_payload · u16 seqlen[K] · u8 payload[]` (the ±flank window bytes inline).

## 2. simplify  (`run_simplify`, simplify.cpp)  — compaction (in RAM) then bubble zip
**Compaction** (`build_super_graph`, compact.cpp), built straight into RAM — **never written to
disk**: walk each maximal same-colour chain (in/out-degree 1), splice the inline node bytes via
exact KMP suffix/prefix overlap (`suffix_prefix`; `k==0` is FATAL — adjacent L5 cores always
overlap ~85 bp), emit one super-node per chain carrying per-base coverage runs. Pure cycles are
preserved as singleton super-nodes. Node schema: `SuperNode{u32 block_len, u32 runs_len,
u8 color}[S]` + concatenated `blocks` (bytes) + `CovRun{u32 end_pos, n_cov, t_cov}[]` +
`super_edges` (IEdge[]); per-node slices found by prefix-summing block_len/runs_len. (Log `[compact]`.)

**Bubble zip** — emit the whole graph, further simplified. Every super-node passes through; each
`std_both` bubble (a `(S1,S2)` colour-3 anchor pair with parallel c1 and c2 arms) is **zipped to
base-pair resolution**: the prefix/suffix common across **all** of an anchor's bubble arms is
absorbed into the anchor (grown, only when "private"), and each arm becomes one residual
super-node holding just its divergent middle. Equal alleles collapse to a direct edge. Coverage
is max-merged where alleles share a region. Output schema = the compacted node schema + edges
with an **explicit overlap**: `GEdge{u32 from, to, overlap}` — reconstruct a path by dropping
`overlap` head bytes of `to` before concatenating. Bubble seams abut (overlap 0); pass-through
edges recompute the KMP overlap on output bytes. **v1 scope: `std_both` only** — orphan buckets
(tumor-only INS, normal-only DEL) and tangles (multi-path, cycles) pass through unzipped.
→ `simplify.bin`. (Log `[simplify]`.) Memory: arm coverage is kept RLE (not flat per-base) and
only read-slices are expanded; adjacency is CSR; the bubble BFS/DFS use generation-stamped scratch.

> "Bubble" is the graph-topology object the zip operates on (anchor pair + two allele arms) — the
> core of this stage, unrelated to any TSV.

## Source layout (`src/`)
`graph.cpp` (graph stage + readers/pools/bloom cascade/`CoreCache`/checkpoint IO), `compact.cpp`
(`build_super_graph`), `simplify.cpp` (`run_simplify`), `extract_fast.cpp` (biolcp SoA adapter,
the LCP hot path), `main.cpp` (CLI), `MurmurHash3.*`, `svp.h` (shared PODs + IO + `suffix_prefix`),
`graph.h` (stage decls).

## Build flag
`-DLARGE_GENOME` is load-bearing and must match `lib/lcptools`'s build (64-bit core start/end on a
whole-genome reference). Rebuild the static lib if you see segfaults on struct fields. svp links no
zlib — all decode is rapidgzip (installed by `make install-deps` into `lib/rapidgzip/venv`).

## Known gaps / next steps
- **zip v1:** orphan buckets and tangles pass through unzipped; orphan INS/DEL zipping against the
  backbone detour is the obvious v2.
- No sequence-reconstruction verifier for the `simplify` output graph yet (only structural
  accounting checked).
- Inversions/DUPs limited by dual fw+RC extraction (see FOUNDATIONS.md) — strand folding is the
  precondition for recovering them.
