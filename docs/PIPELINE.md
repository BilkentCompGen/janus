# Pipeline — current stages and algorithms

`svp` is a reference-free somatic SV caller. It builds a 2-colour LCP-core graph directly from
reads (normal vs tumor HiFi), error-corrects it, then calls somatic variants off it. No assembly,
no alignment, no reference in the caller. The graph is **coordinate-free** *and*, since the **v3
sig-only** format, **sequence-free**: a node is a 64-bit core signature plus its per-colour read
counts. Sequence and genome coordinates only re-enter at the very end, in the anchoring step, and
only to write a VCF a benchmark can score.

End-to-end benchmark flow:
`reference → insilicoSV (truth) → Badread (reads) → [svp graph → svp refine] → call_v3 → anchor_v3 → Truvari`.

## The v3 on-disk graph (one format, two files)
`graph.bin` (build output) and `refined.bin` (refine output) are byte-identical in layout, so every
consumer reads either one:

```
u64 K                         # nodes
u64 sig[K]                    # MurmurHash64A of the core's bytes = the node's identity
u32 ncov[K], u32 tcov[K]      # normal / tumor read counts on the core
u64 E
CEdge{u32 src, u32 dst, u32 ncov, u32 tcov}[E]    # directed adjacency + per-colour traversal counts
```

- **No `color[]`, no `seqlen[]`, no payload.** Colour is *derived* on demand:
  `ncov && tcov → 3 shared`, `ncov → 1 normal-only`, `tcov → 2 tumor-only (candidate somatic)`.
- **Edges carry counts** (how many normal/tumor reads crossed that adjacency). Node counts are
  mm-floored at build time; edge counts are raw (min 1) — the gap between the two is what `refine`'s
  edge gate exists to close.
- **No sequence anywhere** ⇒ every "size" the caller produces is a **core count**, not bp. Only
  `anchor_v3` converts to bp, using the average core step it measures on the reference.

## CLI — two resumable, checkpointed stages
The expensive LCP pass is paid in `graph`; `run` skips any stage whose checkpoint already exists.
Checkpoints are atomic (`.tmp` + rename). **`call` and `anchor` are no longer subcommands** — they
live as out-of-band tools (below).

```
svp graph  <normal> <tumor> <dir> [--level 5 --margin 10 --min-mult 3 --threads 16]
svp refine <dir>
svp run    <normal> <tumor> <dir>   # graph → refine, skipping completed stages
```

Each sample argument is a single FASTQ/`.gz` **or a directory of shards** (`.fastq` / `.fq` / `.svb`),
read in parallel.

## 1. graph  (`run_graph_reads`, graph.cpp)  — reads → self-contained `graph.bin`
Two passes over the same reads.

**Decode.** One pipe-reader per `.gz` reads plain FASTQ from a **rapidgzip** decode process
(`open_decode` → `popen`; parallel single-stream inflate, ~1.5 GB/s @ ~4 cores vs zlib's
~120 MB/s/core). `.svb` shards (2-bit packed, `tools/svb_pack`) are `fopen`ed and decoded inline —
no pipe, no gzip. Tunables: `SVP_RAPIDGZIP` (binary path, default `lib/rapidgzip/venv/bin/rapidgzip`),
`SVP_DECP` (decode threads/file), `SVP_FILL_NW` (pass-2 workers), `SVP_NOCACHE`.

**Pass 1 — count (`run_pool_both`).** Both samples in one pool. Workers run **two LCP extractions
per read, forward and reverse-complement** (the parse is not RC-symmetric), drop any core within
`margin` (10) of a read end, and take each core's identity `core_sig` = **MurmurHash64A over the
core's own bytes** (no flank). Each sig is counted in a sharded **bloom cascade** — one cascade per
colour (`min_mult-1` plain bloom filters + an exact survivor map). **Survivors** = sigs seen
**≥ `min_mult`** times; they become the graph's nodes. The per-colour `min_mult` gate happens *inside*
the cascade, **before** the normal/tumor merge (a normal-error sig drops, leaving its locus
tumor-only). Bloom collisions can only over-keep (false survivor); they never drop a real one (FN = 0).
Survivors are merged into one `(sig, n, t)` record and sorted by sig — **node id = sorted index**.

`SVP_COUNT_ONLY=1` stops here and dumps every survivor as `counts.bin` (`u64 K` + `{sig,n,t}[K]`),
so any candidate mm-floor can be studied offline.

**Pass 2 — edge count (`edge_count_pool`).** Re-decode and re-LCP the same reads (same dual-strand
extraction, same margin, same blacklist merge, so the adjacencies match the node set). For each read,
walk its cores in order and, for every consecutive pair where **both** cores survived the mm gate,
increment that `(prev, cur)` adjacency in the read's colour (a non-surviving core breaks the chain).
Per-worker maps are periodically merged into 64 sig-sharded global maps. Afterwards the sig-pair edges
are resolved to node ids (self-loops dropped) and sorted → `CEdge[E]`.

**Repeat-core blacklist (ON by default, `src/blacklist.cpp`).** Cores whose sig occurs more than
*minfreq* times in the reference (built by `tools/build_blacklist`) are repeats. A maximal run of them
is absorbed into its non-blacklisted flanks as two **overlapping** merged cores — `(A+X)` and `(X+B)` —
so `A→X→B` collapses to `A→B` and the bare repeat core never becomes a node (the two merged cores
still share the run's bytes, so the edge is a real overlap). `SVP_BLACKLIST=<path>` overrides the file;
`SVP_NOBLACKLIST=1` disables it. **The identical merge must be applied by anything that re-parses the
reference** (`anchor_v3`, `tools/fasta_cores`) or merged reference cores will not match merged graph
nodes and anchoring silently fails.
*Gotcha:* the default path `outputs/svp/blacklist_l5_mm10.bin` is not in the tree (only
`blacklist_l5_mm2.bin` is), and a missing blacklist is a fail-fast abort — pass `SVP_BLACKLIST=` or
`SVP_NOBLACKLIST=1`.

**Cross-pass core-boundary cache (`CoreCache`) — currently vestigial.** Pass 1 can record each read's
filtered core `(start, width)` pairs, keyed by `(shard, read-index)`, so a second pass could replay
them instead of re-extracting. Nothing consumes it: `edge_count_pool` re-extracts unconditionally, and
it is disabled outright whenever the blacklist is on (which is the default), because the merge rewrites
core boundaries. Don't count on it for a cost estimate.

## 2. refine  (`run_refine`, refinement.cpp)  — `graph.bin` → `refined.bin`
Topology + counts error-correction; **same v3 format in and out**. Each structure is judged by its read
count **relative to local coverage** (never a global gate — that is what an over-aggressive `min_mult`
does, distorting real bubbles). Default passes:

1. **tips** — short dead-end core-chains whose count ≪ the junction they hang off, plus tiny
   low-count isolated contigs (pure error islands). Iterated (removal exposes new tips).
2. **low-cov side-branches** — leaf-ish sprouts whose count is a small fraction of their strongest
   neighbour (recurrent-error cores clinging to a real path).
3. **bubble popping** — the low, error-like minority arm of a short bubble, judged **per colour and
   flank-relative**: the arm is compared against the bubble's anchor (flank) coverage in *each* sample,
   not an absolute count. A balanced arm is a germline het (keep both); an arm present in tumor and
   absent in normal is a candidate somatic (keep); only a clearly imbalanced, non-somatic minority arm
   is popped. Then a tip re-pass.
4. **edge gate** (`pass_edges`) — node counts are mm-floored but edge counts are raw, so a
   "tumor-only" edge can rest on a *single* tumor read: the dominant somatic-junction false positive.
   Two local, per-colour gates: a **hard floor** (`max(ncov,tcov) ≥ REF_EDGE_MM`, default 3 — the edge
   analogue of the node mm floor) and a **relative gate** (drop an edge carrying less than
   `REF_EDGE_FRAC`, default 0.10, of its endpoints' coverage in *every* colour it has). Edges with a
   killed endpoint are dropped regardless.

An `abs_keep` floor (default 10) protects well-covered cores in the tip/low-cov passes. Env-tunable:
`REF_TIP_MAX REF_TIP_FRAC REF_ISLAND_CNT REF_LOWCOV_FRAC REF_BUB_MAX REF_HET_FRAC REF_SOM_FRAC
REF_ABS_KEEP REF_ROUNDS REF_EDGE_MM REF_EDGE_FRAC`. (Log `[refine]`.)

**Opt-in modes — each writes a DISTINCT file, never clobbering `refined.bin`:**
- **`REF_ITER=<N>` → `refined_iter.bin`: the iterative node↔edge convergence loop.** Nodes and edges
  are not independent — an error edge's reads also inflated its endpoints' node counts. So: first a
  one-shot **flow-conservation filter** (`REF_FLOW=1`, default on inside this mode) removing nodes
  whose summed incident edge counts *exceed* their own read count (an impossibility: every read
  through a node is partitioned among its edges — empirically ~95% error cores, real cores violate
  ~0.01%); it must run on still-consistent counts, i.e. before any deflation. Then up to N rounds of:
  gate the edges → **subtract each dropped edge's per-colour counts from its endpoint nodes** → remove
  nodes left at zero coverage in both colours → rebuild adjacency → re-run the node passes. The peel is
  monotonic (dead nodes only accumulate, edges only shrink), so it converges; N is just the bound, and
  the loop exits early when a round drops nothing. Known approximation: a read crosses a node once but
  up to two of its edges, so a node with one kept and one dropped edge is over-deflated.
  `REF_EDGE_CLAMP=1` (clamp each surviving edge to `min(endpoint counts)`) is **off** — measured to
  recolour shared edges to tumor-only once deflation zeroes a node's `ncov` (~55k false-somatic edges);
  the flow filter is the safe alternative because it *drops* over-flow structures instead of
  recolouring them.
- **`REF_TIP_NUKE=1` → `refined_nuketip.bin`**: aggressive tip removal (ignores `abs_keep`, removes
  every tip ≤ `REF_TIP_MAX`). Measured: clears ~116k extra tip errors while eroding ~0.03% of somatic
  cores (short SV-end trims, not lost calls).

## 3. tools/call_v3  — `refined.bin` → `<prefix>.tsv` + `<prefix>.altsigs.bin`
Standalone (not a `svp` subcommand, not built by `make`). Reference-free somatic calling straight off
the refined graph, using **only graph-intrinsic signals** — derived colour, per-node `ncov`/`tcov`,
topology. A somatic allele always involves a tumor-only (colour-2) core. Call sources:

- **BUBBLE** — a two-arm bubble carrying ≥1 tumor arm → `T/N`, `T/S`, `T/T`. The tumor arm is the alt
  allele; type and size come from the signed **arm core-count Δ** = `|alt| − |ref|` cores
  (`|Δ| ≤ SNP_EPS` (0) → SNP, `Δ > 0` → INS of Δ cores, `Δ < 0` → DEL of −Δ cores). Support = the alt
  arm's minimum `tcov`; VAF = mean alt `tcov` / anchor `tcov`.
- **BND** — a tumor-only run that **bridges two shared (real) anchor cores at ~full tumor depth**
  (`CALL_DEPTH_FRAC`, default 0.30) = a novel adjacency with no reference counterpart (so it never
  forms a bubble) → DUP / INV / large-indel / translocation breakend. Recovered by tracing the
  highest-`tcov` tumor path to the nearest shared core on each side (works inside the repeat tangle;
  no degree-1 requirement), run length capped at `CALL_JMAX` (12). Tagged `BND_CNgain` when an anchor's
  tumor/normal depth ratio ≥ 1.5.
- **TIP** — a tumor-only dead-end chain (≤16 cores) → insertion-end / breakend candidate.

Duplicates (including fwd/RC re-emissions of the same allele) are collapsed by a hash of the **alt
core-sig set**; consumed nodes are not reused by a later pass. Env: `CALL_MINSUP` (default 3),
`CALL_JMAX`, `CALL_DEPTH_FRAC`. Outputs `<prefix>.tsv` (per-call features; sizes in **cores**) and
`<prefix>.altsigs.bin` (`{u64 id, u32 nsig, u64 sig[]}` per call — the alt cores' stored sigs).
(Log `[call]`.)

## 4. tools/anchor_v3  — calls + reference → `calls.vcf`
`anchor_v3 <graph_v3.bin> <ref.fasta> <calls_prefix> <out.vcf> [threads=16]` (needs `ref.fasta.fai`).
The calls are coordinate-free; this stage gives them coordinates so a benchmark can score them.

The alt cores are tumor-only and by definition **absent** from the reference — but their shared /
normal-only graph **neighbours** are present. So: for each call, collect its anchor nodes from the
graph; **re-extract the reference contig-parallel** with the *same* level/margin/blacklist merge and
the *same* `core_sig`, both strands; match anchors by sig (a sig seen at ≥2 distinct positions is a
repeat → dropped as ambiguous); bracket the call between its located anchors on the contig where most
of them agree. Core counts become bp via the **average core step measured on that same scan**
(`total_bp / total_cores`) — nothing is assumed about core length. BND records are re-typed and
re-sized against the anchored reference span (a bridged span that re-parses as mostly normal-only
cores is really a DEL; a sub-50 bp span with a long junction is really an INS), a `LOCALSUB` filter
flags junctions that merely reproduce the span they bridge (local SNP/MNP misfires), and same-class
calls at the same locus are deduped keeping the highest-support representative. Emits VCF 4.2 with
`SVTYPE / SVLEN / END / NANCHOR / CALLSRC / SUBTYPE`. (Log `[anchor]`.)

## 5. Scoring
```
scripts/truvari.sh <truth.vcf> <calls.vcf> <ref.fa> [outdir]     # truvari bench --pctseq 0 + FN/FP breakdown
cd simulation && make bench/hard CALL=<calls.vcf>                # same, against truth/hard/out.vcf
```

## Source layout (`src/`)
`graph.cpp` (graph stage + readers / pools / bloom cascade / `.svb` reader / `CoreCache` / checkpoint
IO), `extract.cpp` (biolcp SoA adapter + the MurmurHash64A `core_sig`), `blacklist.cpp` (the shared
repeat `(A+X)→(X+B)` merge), `refinement.cpp` (`run_refine`), `main.cpp` (CLI: `graph`/`refine`/`run`),
`svp.h` (shared PODs + IO helpers), `graph.h` (stage decls), `unordered_dense.h` (vendored flat hash
map). One binary, one Makefile, no CMake/CI. The v3 call/anchor tools live in `tools/` and are compiled
out of band; `anchor_v3` must link the same `extract.o` + `blacklist.o` + `lcptools` so its parse and
sigs are byte-identical to the graph's.

## Build flag
`-DLARGE_GENOME` is load-bearing and must match `lib/lcptools`'s build (64-bit core start/end on a
whole-genome reference). Rebuild the static lib if you see segfaults on struct fields. svp links no
zlib — all decode is rapidgzip (installed by `make install-deps` into `lib/rapidgzip/venv`).

## Known gaps / next steps
- **Sizes are core counts, not bp, until anchoring.** Everything upstream of `anchor_v3` reasons in
  cores; a size in bp only exists after the average-core-step conversion, which is a *mean* — expect
  size error on short events.
- **Short insertions (<50 bp)** are below LCP-core resolution (median INS ~25 bp) → low recall.
  **[historical, v2 payload caller]** ≥50 bp DEL/DUP/INV recovered ~95–98%; the v2 caller reached
  ~85% overall. The sig-only v3 caller is a port of it, **not yet re-benchmarked** — do not quote the
  old numbers for it.
- **DUP vs INV sub-typing** is unresolved: both land as `BND`. The BND/TIP false-positive filter is
  the other open front.
- Inversions/DUPs are limited by dual fwd+RC extraction (see `FOUNDATIONS.md`) — strand folding is the
  precondition for cleaner recovery, and cannot be done by per-emission canonicalization.
