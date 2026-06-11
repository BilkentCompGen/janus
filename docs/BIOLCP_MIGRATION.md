# lcptools → biolcp migration + the label-signature investigation (2026-06-11)

Switched the LCP backend from the in-tree `perf-core-extraction` fork to lcptools'
**`biolcp`** branch (author's DNA-optimised rewrite), and investigated replacing the
flank-window murmur signature with biolcp's native core **label**, as requested.

**Outcome:** biolcp is in (drop-in, verified identical recall). The bare-label signature
is **not viable** — it breaks the byte-splice reconstruction and craters recall (measured).
Separately, the lcp stage is **gzip-inflation-bound**, so neither biolcp nor the label can
make the run "much faster" with the current single-file input — the real lever is elsewhere.

## What shipped

- `lib/lcptools` on branch **`biolcp`**, built in-place by `make install-deps-lcptools`
  with `-O3 -march=native -DLARGE_GENOME` (its own Makefile installs to `/usr/local`, needs
  root; we archive in-place instead). `-DLARGE_GENOME` = 64-bit core `start`/`end` so posmap
  over the 3.2 Gbp reference cannot truncate.
- `src/extract_fast.cpp` — a thin SoA adapter exposing the existing `extract_fast(...)`
  interface on top of biolcp's `init_lps_offset` + `lps_deepen`. It self-initialises the
  alphabet (`LCP_INIT` once) and converts biolcp's **inclusive** core end → the **exclusive**
  end the rest of svp assumes. graph.cpp / cores.cpp are otherwise unchanged; the murmur
  signature stays.
- `svp` CXXFLAGS gain `-DLARGE_GENOME` (ABI match with the lib) and `-lm` (biolcp uses libm).

**Verification (mm=10, full Sim-it):** biolcp+murmur reproduces the baseline exactly —
geq50 TP=99 / **FP=0** / recall 0.495 / prec 1.000; indel TP=840 / recall 0.894; SNP
TP=15489 / recall 0.913; 16,762 calls (baseline 16,768). modea/modeb reported
`zero-overlap-breaks=0` / `zero_ovl=0`, i.e. the FATAL invariant holds. So biolcp produces
LCP cores equivalent to the old fork — a safe swap.

## Why the bare label fails (measured, not theorised)

biolcp's `struct core` exposes a 32-bit `label` and a 64-bit `bit_rep`. Probing 300K real
reads (`tools/biolcp_probe.cpp`) at L5:

| candidate sig | distinct (of 9.45M true cores) | merge loss |
|---|---|---|
| `bit_rep` | 5,689 | **99.94%** (truncated to ≤63 bits at L5 — degenerate) |
| `label` | 9.25M | 2.1% (mostly structural — it hashes only boundary child-labels) |

So `bit_rep` is unusable and `label` looks good *as a content id*. But the label is
**content-only**, and that is exactly the problem. The old signature was
`murmur64(core ± 10bp flank window)` — hashing the local **context** makes each node ≈ one
genomic locus, so graph-adjacent nodes always share byte overlap. modea's unitig splice and
modeb's arm splice both rely on that overlap (both `abort()` on zero overlap). The label
collapses identical core content from anywhere in the genome into one node, so a degree-1
graph edge can now join byte windows from unrelated loci with **zero overlap**. No flank size
fixes this — it is the content-addressing itself.

Full mm=10 label run (with modea/modeb taught to skip zero-overlap instead of crashing):

| metric | bare label | murmur (baseline) |
|---|---|---|
| modea zero-overlap chain edges | 193,218 / 86.3M (0.22%) | 0 |
| modeb bubbles skipped (non-contiguous arm) | **25,824 / 45,250 = 57%** | 0 |
| posmap sig collisions | 7.78% | 5.74% |
| **≥50bp SV recall / precision** | **0.045 / 0.0098** | 0.495 / 1.000 |
| **indel recall / precision** | **0.20 / 0.93** | 0.894 / 0.966 |

A 0.22% per-edge zero-overlap rate compounds across multi-node arm paths into 57% per-bubble
loss (bubbles sit at the complex/repetitive loci where content-addressing collapses hardest),
and the surviving calls are mostly mispositioned (precision 0.0098). The label is dead.

**If native labels are still wanted** (the author's intent — skip the byte re-hash), the way
to keep locality is a **label-context** signature: hash the core's label together with its
neighbours' labels (a few u32s, cheaper than the byte murmur) instead of the bare label.
That preserves position-discrimination, but it re-introduces a context window — i.e. it does
not "get rid of the ±flank" in spirit. Not done here; flagged for discussion.

## The real speed bottleneck: gzip inflation, not LCP

The mm=10 run log shows lcp **worker utilisation 38.5%**, "reader(inflate)-bound": the single
zlib reader spends ~315 s/pass inflating the 22 GB input, and the file is re-inflated 3× (two
pool passes + one fill pass) ≈ **91% of the 1036 s lcp stage**. The 15 LCP workers sit idle
62% of the time. So optimising LCP extraction (biolcp) or dropping the per-core murmur (label)
**cannot** move wall-time much — the workers are already starved on gzip.

Levers that would actually speed the run (all orthogonal to the signature):
1. **Shard the input** — the caller already reads a *directory* of FASTQ shards with one
   reader thread each (`collect_shards`); a single `.fastq.gz` gets one reader. Splitting each
   sample into N shards gives N parallel inflators (gzip can't parallelise inflation of a
   *single* stream, so this is the only way). Biggest, lowest-risk win, and **already measured
   here** (2026-06-01, whole genome): single file + 8 threads (2 readers) = **19:44**, vs
   62 shards + 32 threads (8 readers + 24 workers) = **7:47** (`outputs/svp/wg_sharded_*.log`,
   `wg_shardmap_*.log`) — ~2.5×. This is where "much faster" comes from, not the signature.
2. **Faster decompressor** — swap zlib for libdeflate/igzip (~2–3× inflate throughput).
3. **Fewer inflation passes** — the count-then-fill design re-reads the input 3×; caching is
   RAM-prohibitive (~120 GB inflated), but a one-pass restructure would help (touches run_lcp).

## Reproduce

- `tools/biolcp_probe.cpp` — collision + speed + overlap probe (gitignored scratch).
- `outputs/svp/biolcp_label/` — bare-label run + benchmark (recall 0.045).
- `outputs/svp/biolcp_murmur/` — verified biolcp+murmur run + benchmark (= baseline).
