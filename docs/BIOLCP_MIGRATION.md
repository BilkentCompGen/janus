# lcptools → biolcp migration + signature experiments (2026-06-11/12)

Switched the LCP backend to lcptools' **`biolcp`** branch and tested two requested signature
changes: (1) drop the ±10bp flank, (2) use biolcp's native core **label** instead of the
flank-window murmur. Both were run end-to-end on the full Sim-it data (mm=10) and benchmarked.

**Result:** biolcp is a real extraction win and ships (= baseline recall). **Dropping the
flank fails** and **the label fails** — both measured below. The flank is load-bearing; it is
not redundant. The earlier "content-addressing collapses repeats" explanation in this file was
wrong — see the corrected mechanism below.

## What shipped (committed)

- `lib/lcptools` on **`biolcp`**, built in-place by `make install-deps-lcptools`
  (`-O3 -march=native -DLARGE_GENOME`; its Makefile installs to /usr/local, needs root).
- `src/extract_fast.cpp` adapts biolcp's `init_lps_offset`/`lps_deepen` to svp's SoA
  `extract_fast` interface (inclusive→exclusive end, one-time `LCP_INIT`, thread-local SoA).
- **Signature unchanged: `murmur64(core ± 10bp flank window)`.** biolcp speeds extraction; the
  signature and the whole downstream pipeline are byte-for-byte the proven version.

biolcp is genuinely faster at extraction: at 8 threads the old `extract_fast` made the pool
**LCP-bound** (94.5% util, readers 32% idle on queue-full); with biolcp the workers finish
faster and the pool becomes **reader(inflate)-bound** (pool_n: 83.3% util, reader-wall 318s,
queue-full 5.3s). The bottleneck moved off LCP onto the gzip reader — that is biolcp working.

## The three signatures, measured (mm=10, Sim-it, truvari)

| signature | ≥50bp recall / prec | indel recall / prec | SNP recall | modea zero-overlap | posmap collisions |
|---|---|---|---|---:|---:|
| `murmur64(core+flank)` — baseline | 0.495 / 1.000 | 0.894 / 0.966 | 0.913 | **0** | 5.74% |
| `murmur64(core)` — no flank, 64-bit full content | 0.385 / 0.027 | **0.416** / 0.586 | **0.40** | 58,011 | 5.96% |
| biolcp `label` — 32-bit, boundary-only | 0.045 / 0.010 | 0.20 / 0.93 | — | 193,218 | 7.78% |

Both flank-free variants crater recall. So the flank is necessary, and using the label gains
nothing (the murmur is not a bottleneck — §speed).

## Corrected mechanism (why the flank is load-bearing)

The signature is `murmur64` over the core **plus ±10bp of flanking context**. That context
makes the signature **position-discriminating**: two genomic loci with the same LCP core but
different neighbours get **different** sigs, so each graph node ≈ **one locus**. Consequently
a node's stored ±flank window bytes are **consistent** (they always come from that one locus).

modea/modeb reconstruct alleles by KMP suffix/prefix splicing of those stored bytes. After
`min_mult` pruning, a chain edge can connect two surviving cores with a pruned core **between**
them — a coordinate gap; the splice only overlaps because the ±10 flank extends across the gap.
With the baseline sig the node's stored flank IS that locus's flank, so the bridge holds → **0
zero-overlap**.

Drop the flank (either `murmur64(core)` or the label) and same-core/different-context loci
collapse into one node whose stored window is from an **arbitrary** occurrence. A skip-gap
chain edge from a *different* read then splices against the wrong flank → no overlap → a seam.
Measured: 58,011 (core-only) / 193,218 (label) zero-overlap modea edges vs 0 baseline. Those
seams corrupt the spliced allele sequences → recall craters. (posmap collisions barely move —
5.96% vs 5.74% — so positioning ambiguity is NOT the cause; the byte-splice corruption is.)

The label additionally is a poor identity: biolcp's `label` is `MurmurHash3_**32**` over only
the *first / second-last / last child + length* (`core.c init_core3`) — boundary-only and
32-bit, so genuinely different cores collide (measured 2.1% on a 300K-read probe: 9,448,997
distinct true contents → 9,248,176 distinct labels). `bit_rep` is worse — truncated to ≤63
bits at L5, only 5,689 distinct values. Neither is the "global unique id" it's described as.

## Speed: it is NOT gzip-bound, and sharding is a documented dead end

Corrected against `docs/PERFORMANCE.md §4–5` (measured 2026-06-02, the authoritative record):

- At the 8-thread sweet spot (old design) the pool is **94.5% utilized and LCP-bound**, not
  gzip-bound — readers sit 32% of wall on queue-full, i.e. they produce faster than workers
  consume. The "low CPU%" comes **entirely from the 287s single-threaded merge tail**
  (radix-sort sigs 64s + edges 128s + graph resolve 27s + adjacency 14s), memory-bandwidth-bound.
- **Sharding plateaus:** 62 shards + 32 threads → 895% CPU = **9 effective cores of 32**;
  "brute-force parallelism plateaus around 8 cores; beyond that, algorithmic changes are
  required." Net benefit too small for the 346 GB / 8:59 decompress cost.
- The current **per-sample two-pass** prune (commit 9fb6503) inflates normal→tumor
  **sequentially** (`pool_n` 339s → `pool_t` 337s, one reader each, reader-bound) plus a heavy
  `fill_pool` — three single-reader passes, vs the old design's one pass with two parallel
  sample readers.

The real speed levers (none is the signature, none is "just shard"): **parallel/bypassed
merge** (256 thread-local `sig&0xFF` partitions, no global sort — documented, deferred),
**un-serializing the per-sample passes**, and a faster single-stream decompressor (libdeflate)
— though gzip inflation of a *single* stream is inherently serial.

## Reproduce
- `outputs/svp/biolcp_murmur/` — baseline (ships): recall 0.495 / 0.894 / 0.913.
- `outputs/svp/coresig/` — `murmur64(core)` no flank: recall 0.385 / 0.416 / 0.40.
- `outputs/svp/biolcp_label/` — 32-bit label: recall 0.045 / 0.20.
- `tools/biolcp_probe.cpp` — label vs content collision probe (gitignored scratch).
