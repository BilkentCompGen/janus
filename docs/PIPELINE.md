# Pipeline — current stages and algorithms

`svp` is a reference-free somatic SV caller. It builds a 2-colour LCP-core graph directly from
reads (normal vs tumor HiFi FASTQ) and progressively simplifies it. No assembly, no alignment.
Output is **coordinate-free**: per-node sequence + per-base coverage, not genome positions.

End-to-end flow:
`reference -> insilicoSV (truth) -> Badread (reads) -> FASTQ -> [svp: lcp -> graph -> modeb]
-> Truvari (evaluation)`.

## CLI — four resumable, checkpointed stages
The expensive LCP pass is paid once; later stages are seconds-to-minutes and re-runnable.
`run` skips any stage whose checkpoint already exists. Checkpoints are atomic (`.tmp` + rename).
A schema change makes an old `modeb.bin` FATAL cleanly on a size mismatch -> regenerate.

```
svp lcp     <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5 --flank 10 --margin 10 --min-mult 3 --threads 16]
svp graph   <dir>
svp modeb   <dir>     # Mode A unitig compaction (in RAM) then Mode B bubble-zip -> modeb.bin
svp bubbles <dir> [--out <dir>/bubbles.tsv]
svp run     <normal.fq.gz> <tumor.fq.gz> <dir>   # runs all, skipping done stages
```

## 1. lcp  (`run_lcp`, graph.cpp)  -- reads -> per-sig multiplicity + within-read edges
Reader-thread pool work-steals FASTQ shards (block gzread + manual `\n` scan) into a bounded
queue drained by workers. **Per read, two LCP extractions: forward and reverse-complement**
(remap+reverse through a 256-byte complement LUT). Each surviving core's sig = `murmur64(core
+/- flank, seed 0)`. Runs **per-sample two-pass** (normal then tumor; halves peak pool RAM),
each sample filtered at `count >= min_mult` **before** the merge. A second pass (`fill_pool`)
re-LCPs to fill each kept sig's window bytes and emit per-worker-deduped adjacent-sig edges.
`--min-mult` is an **lcp-stage** flag (prune runs before `counts.bin` is written), not a graph
flag. `pre_prune.bin` is a resume checkpoint (skips gzip+pool+merge if present).
-> `counts.bin` (u64 N, N x SigCount{u64 sig,u32 n,u32 t}, u64 total, N x u16 seqlen, payload
   = +/-flank window bytes) + `edges.bin` (u64 M, M x IEdge{u32,u32}) + `pre_prune.bin`.

## 2. graph  (`run_graph`, graph.cpp)  -- colour each sig
**1 = normal-only, 2 = tumor-only (candidate somatic), 3 = shared.** Map sig-edges -> node-id
edges, dedupe. Seconds. -> `graph.bin` (u64 K, color[K], ncov[K] u32, tcov[K] u32, u64 E, IEdge[E]).

## 3. modeb  (`run_modeb`, modeb.cpp)  -- Mode A unitig compaction (in RAM) then Mode B zip
**Mode A** (`build_super_graph`, modea.cpp), built straight into RAM — **never written to
disk**: walk each maximal same-colour chain (in/out-degree 1), splice the inline node bytes via
exact KMP suffix/prefix overlap (`suffix_prefix`; `k==0` is FATAL -- adjacent L5 cores always
overlap ~85 bp), emit one super-node per chain carrying per-base coverage runs. Pure cycles are
preserved as singleton super-nodes. Node schema: `SuperNode{u32 block_len, u32 runs_len,
u8 color}[S]` + concatenated `blocks` (bytes) + `CovRun{u32 end_pos, n_cov, t_cov}[]` +
`super_edges` (IEdge[]); per-node slices found by prefix-summing block_len/runs_len.

**Mode B** -- emit the whole graph, further simplified. Every Mode A super-node passes through;
each `std_both` bubble (a `(S1,S2)` colour-3 anchor pair with parallel c1 and c2 arms) is
**zipped to base-pair resolution**: the prefix/suffix common across **all** of an anchor's bubble
arms is absorbed into the anchor (grown), and each arm becomes one residual super-node holding
just its divergent middle. EQ bubbles (identical alleles) collapse to a direct edge. Coverage is
max-merged where alleles share a region. Output schema = the Mode A node schema + edges with an
**explicit overlap**: `GEdge{u32 from, to, overlap}` -- reconstruct a path by dropping `overlap`
head bytes of `to` before concatenating. **Bubble seams abut (overlap 0)**; pass-through edges
keep the Mode A overlap. An anchor is grown only when every edge on the grown side is a zipped
arm (else growing corrupts a surviving backbone edge's overlap). **v1 scope: `std_both` only** --
orphan buckets (tumor-only INS, normal-only DEL) and tangles (multi-path, cycles) pass through
unzipped. -> `modeb.bin`.
Memory: arm coverage is kept RLE (not flat per-base `uint32`) and only read-slices are expanded;
adjacency is CSR; the bubble BFS/DFS use generation-stamped scratch (no per-call O(S) alloc).

## 4. bubbles  (`run_bubbles`, graph.cpp)  -- plain BFS connected-components
On the non-shared nodes of `graph.bin` (backbone colour-3 never traversed). Each component is
one bubble; type `D` (colour-1 only), `I` (colour-2 only), `C` (both). Operates on `graph.bin`,
independent of modeb. **Must stay plain BFS** -- chain+union-find-by-anchor over-merges and
drops recall. -> `bubbles.tsv`.

## Source layout (src/, ~1630 LOC)
`graph 639 / modeb 410 / modea 182`. `svp.h` holds the shared POD (`Core`, `SigCount`, `IEdge`)
+ LCP primitives + checkpoint IO (`xopen`/`wr`/`rd`/`atomic_close`). `extract_fast.cpp` adapts
biolcp's `init_lps_offset`/`lps_deepen` to svp's SoA extractor (the caller's hot path).

## Build flag
`-DLARGE_GENOME` is load-bearing and must match `lib/lcptools`'s build (64-bit core start/end
on a whole-genome reference). Rebuild the static lib if you see segfaults on struct fields.

## Known gaps / next steps
- **modeb v1:** orphan buckets and tangles pass through unzipped; orphan INS/DEL zipping against
  the backbone detour is the obvious v2.
- No sequence-reconstruction verifier for the Mode B output graph yet (only structural
  accounting checked).
- Inversions/DUPs limited by dual fw+RC extraction (see FOUNDATIONS.md) -- strand folding is the
  precondition for recovering them.
