# Notes — measured numbers and deferred ideas

Every number here says where it came from and when. Nothing in this file is an estimate.

## The cost budget (hcc1954, 16 threads, machine idle)

Source: `docs/RESULTS.md`, snape arm re-measured end to end 2026-08-13 with the shipped binary.

| stage | wall | peak RSS | share of 39:01 |
|---|---|---|---|
| `graph` (both read sets, 391 Gbp) | 10:21 | **32.8 GB** | 26.5% |
| `bubbles` | 04:02 | 2.9 GB | 10.3% |
| `junctions` | 02:06 | 2.9 GB | 5.4% |
| minimap2, 222,387 contigs / 2.04 Gbp | 04:00 | 19.1 GB | 10.3% |
| `call` (candidates) | 00:02 | 0.04 GB | — |
| `pick` | 02:01 | 0.12 GB | 5.2% |
| minimap2, 212,754 normal reads / 3.79 Gbp | **16:25** | 20.9 GB | **42.1%** |
| `call` | 00:05 | 0.04 GB | — |
| **total** | **39:01** | **32.8 GB** | |

Two things follow. The peak RSS is `graph`'s and nothing else comes near it. The wall is the
aligner's: the two minimap2 calls are 20:25, 52.3%, and the normal alignment alone is 42.1% — it
spends 16:25 on 3.79 Gbp while the contigs' 2.04 Gbp cost 4:00, i.e. **4x per base**. That ratio has
not been explained and is the largest single unexamined item in the pipeline.

Throughput per base, from the same table:

| stage | Gbp | wall | Mbp/s | strands parsed |
|---|---|---|---|---|
| `graph` | 391 | 621 s | 630 | both (`extract2`) |
| `junctions` | 196 | 126 s | 1556 | one (`extract`) |
| `pick` | 195 | 121 s | 1612 | one (`extract`) |

`graph` is ~2.5x slower per base because it reads twice the data, parses both strands of it, and
inserts where the other two only look up. Against the parse rate the code records (one core of
`extract2` = 65 Mbp/s, so single-strand `extract` ~130), `pick` runs at 83% of what its 15 workers
could parse and `graph` at 65% — the gap is `graph`'s insert and merge work.

## `graph` internals (exp1 core sweep, t=16, 2026-08-28)

Source: `results/exp1/sweep_t16.log`. Input 600,000 reads per sample, ~10.2 Gbp each — a subsample,
not the full pair.

```
sample 1/2: 71.385 M distinct sigs -> 35.035 M >= 3   (49.1%)
sample 2/2: 66.723 M distinct sigs -> 32.688 M >= 3   (49.0%)
merged nodes 48.773 M | adjacencies 48.407 M
core emissions 412.585 M
count phase: wall 57.1 s, CPU 835.3 s, peak RSS 20.48 GB
  of that CPU: extract 345.2 s (41.3%), flush 73.4 s (8.8%), UNACCOUNTED 416.7 s (49.9%)
sort_edges 0.5 s, write 4.5 s
```

Two findings worth acting on:

1. **Half the count-phase CPU is neither the parse nor the flush.** It is the per-core worker-local
   map insert, `seq_put`, and the edge-set insert. No counter exists for it.
2. **Half of every sample's distinct signatures die at the `min_mult` gate** (49% pass). They are
   read errors, and the exact per-sample table holds all of them until the sample ends. This is the
   most likely single driver of the 32.8 GB peak, and it is not measured — only inferred.

`eff` in that table divides CPU by the box's core count (32 from affinity), not by `--threads`. At
`--threads 16` the reported 45.7% is 14.6 busy threads out of 16, i.e. ~91% thread efficiency.
`graph` is not thread-starved; it is doing more work.

## `graph.bin` v4 anatomy (measured from `tmp/night/exp1_out_t16/graph.bin`, 2026-08-28)

S=2, K=48,772,867, E=48,407,350, file 3,768,707,185 B (3.51 GiB). Computed length matches the file
byte for byte.

| section | bytes | share |
|---|---|---|
| header (40 B: magic, version, S, reserved, ov_s, ov_c) | 40 | ~0 |
| `sig[K]` | 390 MB | 10.4% |
| `cov[S][K]` | 97.5 MB | 2.6% |
| `slen[K]` | 195 MB | 5.2% |
| `seq` | 2.51 GiB | **71.6%** |
| `esrc`+`edst` | 387 MB | 10.3% |

- Core raw length: mean **219.8 bp**, median 218, 5–95% 179–267, max 2,474. None zero, none stored
  verbatim in this graph.
- Measured consecutive-core overlap: 50,259,604,774 bp over 410,457,397 pairs = **122.45 bp**.
- Somatic cores (`cov[0]==0`): 13,737,394 = 28.17%; after `drop.sigs` 13,621,249.
- Only `cov[S][K]` grows with the sample count. At S=21 it is ~1 GB; everything else is fixed.
- `seq` exists for one consumer, `bubbles`, which asks about 580,542 + 585,988 = **1,166,530 cores,
  2.4% of K**. The other 97.6% is stored and never read.

## `junctions` output shape (same run)

```
chunks: 1,391,087 total
        333,007 carried by some read
        137,916 with a contig
        195,091 carried but no read had 300 bp of anchor on both sides
        355,026 windows rejected for short flanks
185,183 contigs, 2.39 Gbp from 10.22 Gbp of reads (4.28x less sequence)
179,231 distinct reads produced a contig = 29.9% of the 600,000
```

Contigs per chunk, counted from `ctg.t.fa`: 1 contig 74.3%, 2 19.4%, 3 4.8%, 4 1.0%, **5 (the cap)
0.6%**. The `per_chunk = 5` cap therefore almost never binds; it is a ceiling on alignment volume,
not the mechanism by which a call reaches its support floor.

Support is aggregated across chunks at the alignment stage, not within a chunk. Checked on three
calls from the same run — contigs landing at the call's locus, and their `chunk=` fields:

```
chr17:39956249  SUPP=25   30 contigs -> 14 distinct chunks
chr17:39980058  SUPP=10   15 contigs -> 10 distinct chunks
chr5:181242308  SUPP=10   15 contigs -> 11 distinct chunks
```

## Deferred: a reverse index instead of `seq`

Proposed 2026-08-28. **Not adopted. Current system stands.**

Store, per core, up to N occurrences of `(sample, read_no)` — 6 + 28 bits, 34 rounded to 36 — instead
of the core's 2-bit packed bytes. At 5 occurrences that is 180 bits against the 442 bits/core `seq`
costs today (55.3 B measured above); even 11 occurrences stays under it.

The goal was never file size. It was (1) lower peak RSS and (2) less time, by fetching and parsing
only the reads a stage actually needs instead of streaming every read set again.

What it would cost and buy, against the budget above:

- **RSS.** The peak is `graph`'s 32.8 GB and the index is built in `graph`. Swapping the 55.3 B/core
  `seq` arena for a 22.5 B/core index is ~2.3 GB on 70.4 M nodes, ~7%. `junctions` (2.9 GB) and
  `pick` (0.12 GB) do not touch the peak at all.
- **Time.** `junctions` + `pick` are 4:07 of 39:01 = 10.6%; that is the ceiling. Realistically
  `junctions` 2:06 -> ~1:00 (only 29.9% of reads produce a contig) and `pick` 2:01 -> ~0:20 (it
  scans 11.34 M reads to select 176,912). ~2:40, ~7%.
- **`bubbles` pays it back.** Today `seq` serves it with no IO and no parsing. With an index it must
  fetch a read per core for 1,166,530 cores and, without a stored offset, re-parse each one to find
  the core — and it needs the core's END as well as its start, so the offset alone would not do.
  `derive` is single-threaded today. `bubbles` is 4:02, twice `junctions`.

Blockers that would have to be solved either way:

- `.svb` has no index and its records are variable length, so `(sample, read_no)` is not addressable.
  A per-read offset table is needed: ~8 B x reads, ~100 MB per sample. Cheap, but not optional.
- A 16-bit offset field would cap read length at 65,535. Sampled max was 39,991 bp over 5,817 reads
  of `exp1_fastq/t.fastq`; the true max is unmeasured and overflow would corrupt silently.
- Which occurrences get recorded depends on thread scheduling. Today's equivalent race is harmless
  because the same signature implies the same bytes; offsets and read numbers are not invariant, so
  `graph.bin` would stop being content-reproducible.

Verdict: the arithmetic is sound and the bits are cheaper, but it addresses ~7% of RSS and ~7% of
wall while the aligner is 52.3% of the wall. Revisit only if the alignment cost is dealt with first,
or if the profiling below shows the peak is somewhere the index would actually move.

## Open: profile `graph` from the inside

The 32.8 GB peak and the 49.9% of count-phase CPU that is unattributed are both inferences, not
measurements. Nothing in the code reports per-structure bytes. See the `[eff]` table for what exists
today: wall, CPU and RSS per phase, and nothing below that.

## h2009 alignment, measured 2026-08-28

`results/measure/h2009/mm2_t.time`. Input plain FASTQ (`hifi/h2009/t.fastq`, 377 GB, 12,683,263
reads), 16 threads, machine otherwise idle.

| | |
|---|---|
| minimap2 + `samtools sort` + `samtools index` | **2:21:51** wall, 114,818 CPU-s (1349%) |
| minimap2 alone (its own `[M::main]` line) | 1:59:19 wall, 106,589 CPU-s |
| peak RSS | 19.48 GB (minimap2's own; `sort` stayed below it) |
| output | 75.9 GiB BAM + 22.7 MiB index, exit 0 |

**Comparability trap.** `results/castle_keep/*.time` and `results/measure/*.time` do not measure the
same thing:

- `castle_keep/t.mm2.time` (hcc1954, 4:38:39) timed **minimap2 only**, fed from a **pipe**
  (`snape dump` off `.svb`), with `sort` and `index` outside the measurement — and it ran at 1163%
  CPU, i.e. the box was shared. Its `File system inputs` is 11.3 M blocks; h2009's is 1,285 M.
- `results/measure/h2009/mm2_t.time` wraps `minimap2 | samtools sort && samtools index` and reads a
  plain file.

The like-for-like number for hcc1954 is `RESULTS.md`'s **2:02:57** (minimap2 + sort + index,
196 Gbp, 2026-08-11), against h2009's 2:21:51 on a 3% larger input. Do not put 4:38:39 next to
2:21:51 — that gap is contention and a different command, not a property of the data.

## What `graph` actually reproduces (measured 2026-08-28)

The same binary, the same input (`tmp/night/exp1_samples.tsv`), two runs, compared section by
section:

| section | run to run |
|---|---|
| header, `sig[K]`, `cov[S][K]`, `E`, `esrc`, `edst` | **byte-identical** |
| `slen[K]` | 2,242,044 of 195,091,468 bytes differ (1.15%) |
| `seq`, `sbytes` | differ, including total length (2,698,628,214 vs 2,698,626,852) |

So the node set, the counts and the topology are fully deterministic; what varies is **which read
donated a core's bytes**. A signature is the hash of the *homopolymer-compressed* bytes but `seq`
stores the *raw* ones, so two reads that spell the same core with different run lengths give the
same node different bytes. That is the whole of `CLAUDE.md`'s "not byte-reproducible in its junk-core
tail".

Any change to `graph` is therefore verified as: `sig` / `cov` / `esrc` / `edst` byte-identical to a
baseline run, and `slen` / `seq` differing no more than two baseline runs differ from each other.

## `graph` profiled from the inside (2026-08-28)

Instrumentation added to `graph.cpp`: a second per-read CPU counter for the bookkeeping the loop does
outside `extract2`, a sampled timer on the byte-packing block, exact per-structure byte accounting at
each sample's end, a 200 ms RSS trace, and a count histogram. All additive. Verified: `sig`, `cov`,
`E`, `esrc`, `edst` byte-identical to a baseline run; `slen`/`seq` within the run-to-run band.

### Where the count phase's CPU went (subsample, 600k reads/sample, 16 threads, 832 CPU-s)

| | CPU-s | share |
|---|---|---|
| `extract2` | 347 | 41.7% |
| **byte packing** | **~306** | **36.8%** |
| shard flush | 78 | 9.4% |
| reader thread | 17.3 | 2.1% |
| unattributed remainder | ~84 | 10.1% |

The packing block ran on **400.8 M of 412.6 M core emissions — 97.2%**. The worker-local map is
cleared at every flush, so the same core is "first sight" again and again and pays the full pack each
time. `clock_gettime(CLOCK_THREAD_CPUTIME_ID)` is a syscall here, not a vDSO read, so the sampled
timer must have its own cost calibrated and subtracted (80–90 ns per pair) or the estimate exceeds
the loop containing it.

### The RSS peak is occupancy, not a growth transient

The 200 ms trace is a **ramp**, not a spike: RSS climbs through each sample's pass and drops when
`merge_sample` frees that sample's tables (13.18 -> 9.64 GB, then 19.32 -> 12.22 GB). The earlier
"flat map doubling its dense array" hypothesis is **wrong**; `reserve()` will not buy the 8–13 GB it
predicted.

At the peak (second sample), exact bytes: sample-map 3.00, sample-seq 2.45, sample-edges 3.00,
graph-map 1.50, graph-cov 0.12, graph-seqref 1.00, graph-seq 2.72, graph-edges 1.50, plus
worker-local ≤6.16 GB (sum of per-worker high water). That accounts for the 19.5 GB peak with nothing
left over.

Signature count histogram, per sample: **30.1% seen once, 20.8% seen twice** — 50.9% can never reach
`min_mult = 3`, yet the exact table carries them to the end of the sample. That is the memory lever,
and it needs an exact filter, not a probabilistic one.

### `SNAPE_FLUSH` is a free dial (measured, no code change)

| FLUSH | peak RSS | worker-local | count wall | flush CPU |
|---|---|---|---|---|
| 2^20 (default) | 19.54 GB | 6.16 GB | 57.3 s | 78.0 |
| 2^18 | 17.91 GB | 3.00 GB | 57.2 s | 78.0 |
| 2^16 | 17.26 GB | 2.53 GB | 58.8 s | 86.5 |

2^18 costs nothing measurable and returns 1.6 GB. Not made the default yet: the dedup a bigger table
buys may matter more at 391 Gbp than at 20, and that has not been measured.

## BUG, found and fixed 2026-08-28: reverse-strand cores stored the wrong bytes

`extract2` returns a REVERSE-strand core's `[st, en)` **in the reverse complement's own frame**.
`graph.cpp` sliced the FORWARD read at those offsets, so every core donated by a reverse-strand
sighting — about half the graph — carried a real but unrelated 220 bp span from the mirror position
in the read.

Direct test (`tmp/prof/t_extract.cpp`, synthetic 40 kb read, level 5):

```
forward cores 437, signature reproduced 436     (the 1 miss is the last core's known-bad `en`)
reverse cores 430:
  rc(forward_read[st:en])       <- what graph.cpp did :   0/430
  rc_read[st:en]                                      : 430/430
  rc(forward_read[L-en:L-st])                         : 430/430
```

On the real subsample graph, unpacking each stored core and hashing its homopolymer-compressed bytes
reproduced the node's own signature for **50.5%** of cores. After the fix: **100.0%**.

Lengths were never affected (`w = en - st` is the same in either frame), so `bubbles`' base-span
clause was right all along. What was wrong is its **minimizers**, for half the cores — hence
`drop.sigs`, hence which bubbles are dropped, hence recall. Every accuracy number measured before
2026-08-28 was measured with this defect.

The fix is one line: slice at `len - en[i]` for the reverse strand.

## `graph` optimisation, 2026-08-28

Two changes to `graph.cpp`, both verified against a baseline run of the unmodified binary
(`sig`, `cov`, `E`, `esrc`, `edst` byte-identical; `slen`/`seq` within the run-to-run band):

1. **LUT-driven packing.** `b2()` was a four-way ternary chain called twice per base — once to scan
   for a non-ACGT byte, once to pack — and the reverse strand paid two further per-byte loops
   (reverse in place, then complement). Now one 256-entry LUT, one pass, four bases per store, and
   the escape test falls out of the same loads; the reverse complement is built in a single pass.
2. **The reverse-strand slice fix** (see the BUG section above).

Measured on the subsample (600k reads/sample, 16 threads), count phase:

| | wall | CPU-s | per-core bookkeeping |
|---|---|---|---|
| before | 57.3 s | 832.3 | 306.2 |
| after | **44.8 s** | **643.8** | **121.6** |

`graph` end to end 62.8 -> 49.8 s, **-20.7%**. Packing went from ~306 to ~96 CPU-s; at 88.2 Gbp
packed that is 1.09 ns/base.

Tried and worth recording: removing the 220-byte `cb.assign` copy in front of the packer — 88 GB of
memcpy on this subsample — is worth **nothing** (121.6 against 123.6 CPU-s). The cost is the packing
loop, not the copy feeding it. The one-pass form was kept because it is shorter, not because it is
faster.

What is left in the count phase: `extract2` ~347 CPU-s (54%), bookkeeping 122 (19%), flush ~80 (12%),
reader 17. The packing block still runs on **97.2% of all core emissions**, because the worker-local
map is cleared at every flush and the same core is "first sight" again and again. Making it *rare*
rather than *fast* is the next idea, and it needs a way to know a core will survive before its bytes
are gone.

### Full scale changes the memory picture (h2009 normal, 13.4 M reads)

```
180.069 M distinct -> 69.107 M >= 3 (38.4%)
  seen once  109.371 M (60.7%)
  seen twice   1.591 M ( 0.9%)
before merge: sample-map 10.00 | sample-seq 5.27 | sample-edges 8.00 GB, RSS 21.79 GB
```

On the subsample the split was 30.1% / 20.8%; at full coverage it is **60.7% / 0.9%**. The subsample
is a poor proxy for this question — errors are ~unique, so their share of the distinct set grows with
depth while the true cores saturate.

61.6% of every sample's exact table can never reach `min_mult = 3`. Dropping a signature the first
time it is seen and admitting it on the second — with the entry created at count 2 — keeps every
survivor's count exactly right. The obstacle is knowing "seen before" without a table: a Bloom filter
does it, but the worker map is cleared per flush, so a second filter is needed to remember which
signatures have already been admitted, and a false positive THERE undercounts by one and can drop a
true count-3 core. That risk has not been quantified and the change has not been made.

### Full-scale validation, h2009 (26.1 M reads, 771 GB FASTQ, 16 threads)

`results/measure/h2009/snape_graph.{time,log}`, optimised binary, 2026-08-28.

```
graph          12:26 wall, 1518% CPU, peak RSS 34.56 GB
nodes 70,518,415  edges 75,492,862  core sequence 3.673 GiB
count phase    738.8 s wall, 11,321.2 CPU-s
  extract2     7,087.7 (62.6%)   bookkeeping 2,484.8 (21.9%, of which packing ~1,975)
  flush        1,182.8 (10.4%)   reader 361.7 (3.2%)   workers blocked 0.3 s
```

Against the last full-scale run of the old binary (hcc1937, `results/castle_keep/hcc1937.snape.run.log`),
normalised per G core emissions — 7.23 G there, 8.32 G here:

| | old | new | per emission |
|---|---|---|---|
| count CPU-s | 13,077.7 | 11,321.2 | **−24.8%** |
| `extract2` | 6,198.3 | 7,087.7 | −0.6% (unchanged, as it should be) |
| bookkeeping | 5,930.5 | 2,484.8 | **−63.5%** |
| count wall | 872.6 s | 738.8 s | **−26.4%** |
| peak RSS | 33.53 GB | 34.56 GB | ~unchanged |

The two pairs are not the same data, so only the per-emission column is a comparison. `extract2`
landing within 0.6% is the control: it was not touched, and it did not move.

### Per-stage, h2009, optimised binary (each stage run alone under `/usr/bin/time -v`)

| stage | wall | peak RSS | CPU |
|---|---|---|---|
| `graph` | 12:26 | 34.56 GB | 1518% |
| `bubbles` | **0:09.06** | 7.57 GB | 151% |

`bubbles` at 9 seconds against `RESULTS.md`'s 04:02 is not a 27x speed-up from anything measured
here. The likely explanation is that the 2026-08-13 cost table predates `graph.bin` carrying core
sequence, when `bubbles` still re-read every FASTQ to recover the arms' bytes. If so that row of the
cost table is stale, and these per-stage numbers replace it. Worth confirming on a second pair.
| `junctions` | 2:26.63 | 1.97 GB | 1317% |

`junctions` on h2009: 12,683,263 reads (201.81 Gbp) -> **756,552 contigs, 7.07 Gbp, 28.5x less
sequence**. Chunks 272,369 total, 235,641 carried by some read, 228,249 with a contig, 7,392 carried
with no read having 300 bp of anchor on both sides, 164,953 windows rejected for short flanks.

`bubbles` now drops far more: 609,330 of 686,429 clean bubbles, 867,326 cores. Against the same pair
before the reverse-strand fix the contig count rose 666,659 -> 756,552 (+13.5%) — removing cores
splits chunks, and the per-(chunk, line) cap of 5 then admits more windows.
| contig alignment (756,552 contigs / 7.07 Gbp) | 16:05.44 | 17.80 GB | 1361% |
| `call` (candidates, floor 4) | 0:06.31 | 0.15 GB | 99% |

Candidates: 25,669 junctions + 29,961 insertions clustered, 1,484 dropped for a repetitive anchor,
**3,628 candidates** (1,844 junction + 1,784 insertion).
| `pick` | 2:11.76 | 0.11 GB | 1530% |

`pick` on h2009: 355,426 of 13,430,557 normal reads, **2.65%**, 11.57 GB FASTQ.
| normal alignment (355,426 reads / 11.57 GB FASTQ) | 12:11.10 | 18.64 GB | 1353% |
| `call` | 0:10.11 | 0.15 GB | 99% |
| **total** | **45:46** | **34.56 GB** | |

Final call set: 383 junctions + 566 insertions = **949**. 2,050 dropped by the germline check,
1,167 by the repetitive-anchor filter.

### Accuracy after the reverse-strand fix (h2009, CASTLE `minda_4_2`, 925 records)

| | precision | recall | F1 |
|---|---|---|---|
| `RESULTS.md`, buggy binary | 0.854 | 0.876 | **0.865** |
| today, fixed binary | 0.850 | 0.872 | **0.861** |

TP 807, FN 118, FP 142 over 949 calls. The fix **costs 0.004 F1** on this pair: 3 fewer true calls,
4 more false ones.

The likely cause is that `bubbles`' `max_unmatched = 20` was calibrated against the broken
minimizers. With correct ones the LCS of two genuinely similar arms is much higher, so more bubbles
clear the colinearity clause and are dropped as local edits — including real junctions. The threshold
needs re-sweeping now that its input is correct; `SNAPE_BUBBLE_TSV` exists for exactly that, and the
rule can be scored offline instead of re-running the stage per candidate value.

### h1437 from `.svb` (12.5 M + 12.1 M reads, 105 GB packed)

| stage | wall | peak RSS | CPU |
|---|---|---|---|
| `graph` | 12:22.31 | 34.85 GB | 1485% |

nodes 70,291,527, edges 75,199,228, 8,479.7 M core emissions.

Worth noting against h2009's `graph` (12:26.36 from 771 GB of plain FASTQ): the same wall time from
**105 GB of `.svb`**. The input format buys disk, not time — the parse dominates and the reader was
never the bottleneck (workers blocked 0.3 s over the whole h2009 run).
| `bubbles` | 0:07.01 | 7.17 GB | |
| `junctions` | 2:05.25 | 1.84 GB | |

h1437 `junctions`: 12,080,830 reads (199.82 Gbp) -> 392,414 contigs, 3.66 Gbp, **54.7x less
sequence**. Chunks 159,696 total, 132,884 with a contig. `bubbles` dropped 370,406 of 413,858
bubbles, 471,630 cores.

`bubbles` at 7 s here and 9 s on h2009 settles it: the 04:02 in `RESULTS.md` is from before
`graph.bin` carried core sequence, and that row of the cost table must be replaced.
| contig alignment (392,414 contigs / 3.66 Gbp) | 8:57.29 | 17.51 GB | |
| `call` (candidates, floor 4) | 0:03.29 | 0.06 GB | |

h1437 candidates: 15,057 junctions + 14,683 insertions -> 1,702 (1,036 + 666).
| `pick` | 2:16.86 | 0.11 GB | |
| normal alignment (217,462 reads) | 9:44.51 | 18.64 GB | |
| `call` | 0:05.83 | 0.06 GB | |

### Accuracy after the fix, both pairs (CASTLE `minda_4_2`, min_size 50, tolerance 500)

| pair | before (RESULTS.md) | after | ΔF1 |
|---|---|---|---|
| h2009 | 0.854 / 0.876 / **0.865** | 0.850 / 0.872 / **0.861** | −0.004 |
| h1437 | 0.727 / 0.783 / **0.754** | 0.665 / 0.783 / **0.719** | −0.035 |

h1437 keeps every true call (TP 173 either way) and gains 22 false ones; h2009 barely moves. The fix
is not in question — the stored bytes now reproduce their own signature 100% of the time against
50.5% before — but the call set got worse at the shipped thresholds.

**`max_unmatched` is not the lever.** Swept with `bubbles` alone (7–9 s per point):

| unm | h2009 bubbles dropped | h1437 |
|---|---|---|
| 0 | 282,409 (41.1%) | 194,696 (47.0%) |
| 2 | 360,361 (52.5%) | 240,892 (58.2%) |
| 5 | 470,934 (68.6%) | 297,924 (72.0%) |
| 10 | 572,416 (83.4%) | 349,532 (84.5%) |
| **20 (shipped)** | 609,330 (88.8%) | 370,406 (89.5%) |
| 60 | 611,893 (89.1%) | 372,385 (90.0%) |

The curve saturates by 20, so there is nothing above it, and going below it drops FEWER cores, which
leaves more somatic material and produces MORE contigs, not fewer — measured on h1437 through
`junctions`: unm 0 → 759,670 contigs / 7.11 Gbp, 2 → 685,109, 5 → 551,763, 10 → 431,563,
20 → 392,414 / 3.66 Gbp. Lower `unm` therefore costs alignment time and, on the FP evidence, would
add candidates rather than remove them.

So the recalibration this needs is not a slide along `max_unmatched`. Left for a decision.

Type breakdown, for what it is worth: the CASTLE truth sets carry **no INV records at all**, yet
snape reports 68 INV on h2009 (51 counted TP by position, 17 FP) and 64 on h1437 (40/24). That is
not new and not caused by this fix, but it means a sixth of h1437's calls are a type the truth set
does not contain.

## Experiment 2 with the fixed binary, conditions A and C (2026-08-28)

`results/exp2_fixed/`, `tmp/night/exp2_fixed.sh`. Only the two pairs whose reads are still on disk.
A = CASTLE's own `minda_4_2` truth, snape scored from outside, no membership premium. C = the same
rule with snape added to the HiFi group as a 12th caller.

F1 under A:

| caller | h2009 (925) | h1437 (221) |
|---|---|---|
| Severus v1.2 HiFi | 0.903 | 0.810 |
| Severus v1.2 ONT | 0.894 | 0.791 |
| nanomonsv HiFi | 0.877 | 0.759 |
| **snape** | **0.861** | **0.719** |
| nanomonsv ONT | 0.848 | 0.709 |
| SAVANA HiFi | 0.748 | 0.791 |
| SAVANA ONT | 0.759 | 0.769 |
| Sniffles2 2.2 | 0.746 | 0.441 |
| Sniffles2 2.0.7 | 0.622 | 0.331 |
| manta | 0.335 | 0.577 |
| gripss | 0.305 | 0.560 |
| svaba | 0.170 | 0.218 |

snape is 4th of 12 on h2009 and 6th of 12 on h1437.

Under C, snape reaches 0.869 (h2009, 936-record consensus) and 0.738 (h1437, 228), i.e. a membership
premium of +0.008 and +0.019 — small, which is the point of reporting A as the headline.

## h1437 raw reads: the download completed and does NOT match ENA (2026-08-29)

`tmp/night/fetch_h1437.sh` finished 00:43:32 and reported a size warning on both files. The warning
itself was spurious — the sizes hard-coded in the script (98,341,262,931 and 108,932,596,051) are
stale, and the live ENA record now gives exactly what was downloaded:

| accession | on disk | ENA `fastq_bytes` | ENA `fastq_md5` | our md5 |
|---|---|---|---|---|
| SRR28305179 (t) | 91,552,855,387 | 91,552,855,387 | a45dfe1f88d2aa96082d86555f38f2de | 29ea4a1e2b26527dd76791ab7b181c97 |
| SRR28305175 (n) | 101,455,462,629 | 101,455,462,629 | 01da31b87cea3ce7d8a76d608a5ebee0 | 7aec0ddd1c18dde25888af6fea6fdb60 |

Right length, wrong bytes, on both files. Both were resumed from `.part` files left by an earlier
session (42 GB and 383 MB), which is the obvious suspect: aria2c resumes on byte offset and verifies
nothing unless it is given a checksum.

`fetch_h1437.sh` never asked for one. That is the process defect: it checks size only, and the size
it checks against was written down by hand months earlier instead of read from the API. Any
re-download must pass `--checksum=md5=<the value the API returns now>` so aria2c fails loudly instead
of leaving a file that looks finished.

**This does not touch the h1437 accuracy results.** Those came from `hifi/h1437/{n,t}.svb`, packed
2026-08-26 from an earlier, separate download. The `.gz` files are needed only for the timing arm,
which must be measured from plain FASTQ.

The `.svb` files are independently corroborated and the h1437 results stand on them: `svb_pack`
recorded 12,080,830 reads / 199,821,630,627 bases for the tumour and 12,505,827 / 220,503,020,172 for
the normal, and ENA's `read_count` / `base_count` for the two accessions are those numbers to the
digit. `graph` then read exactly those read counts.

## The homopolymer pass: measured, and left alone (2026-08-29)

`tmp/prof/bench_hpc.cpp`, `tmp/prof/bench_hpc2.cpp`, on 20,000 real h1437 reads (0.33 Gbp):

```
hpc      2.78 ns/base      extract2 15.92 ns/base      hpc is 17.4% of extract2
```

That is three times what it was assumed to be, and extract2 is 62.6% of the count phase, so the
homopolymer pass is ~11% of it — about 1,233 CPU-s of h2009's 11,321.

Three scalar shapes, same reads:

| | ns/base |
|---|---|
| as shipped, `push_back` per run | 2.41 |
| buffers held at high water, indexed | 2.16 |
| the same, branchless (map store to a sink slot on a continuation) | 2.05 |

15% of a 17% component, i.e. 2.6% of extract2. **Not taken**, and not only because it is small:
`extract.cpp`'s own header warns that `emit()` reads one byte past the last core of a strand and
lands on whatever the buffers left there, so the shipped `clear`/`reserve`/`push_back` sequence is
load-bearing. Changing the buffer discipline can move `en` for that core, hence whether it clears the
margin, hence the graph. Not worth 2.6%.

A SIMD left-pack (compare against the byte-shifted vector, then `pshufb` with a shuffle-mask table
for the bytes and a permute for the 32-bit offsets) would plausibly reach ~0.8 ns/base, ~4% of the
graph stage. It needs intrinsics, which the house rules do not allow, and it has the same buffer
problem.

minimap2 does not help here: `sketch.c` uses a 256-byte `seq_nt4_table` and an inner loop that skips
each run, and it never materialises the compressed string at all -- it consumes the HPC k-mer on the
fly. We cannot, because lcptools needs a contiguous string and `map[]` to get back to raw
coordinates.

## `call` parameters can be swept without re-aligning (2026-08-29)

The contig BAM and the normal BAM are the expensive part of a variant; `snape call` reads them and
nothing else. So every parameter of `call` -- MAPQ floor, support floor, the two repeat guards, the
germline dirtiness fraction, the clustering tolerance -- can be swept on a finished pair in about
30 seconds a point, against the 24 minutes a bubble-rule variant costs. `tmp/night/call_sweep.sh`.

On h1437 over the align@25 alignments, one parameter at a time from the shipped defaults (F1 0.733):

| | P | R | F1 |
|---|---|---|---|
| `min_mapq` 10 -> 30 | 0.715 | 0.783 | **0.747** |
| `anchor_rep` 0.50 -> 0.45 | 0.712 | 0.783 | **0.746** |
| `min_mapq` 20 | 0.698 | 0.783 | 0.738 |
| `shred_rep` 0.40 -> 0.30 | 0.690 | 0.787 | 0.736 |
| `min_supp` 6 | 0.733 | 0.733 | 0.733 |

and the cliffs, which matter as much: `anchor_rep` 0.60 gives 0.548 and 1.00 gives 0.327,
`max_mm` 0.02 gives 0.583 and 0.05 gives 0.314, `min_supp` 3 gives 0.405, `min_mapq` 0 gives 0.494.
`normal_dirty` and the clustering tolerance barely move anything.

Combined:

| | P | R | F1 | TP | FP |
|---|---|---|---|---|---|
| **mapq 40, anchor_rep 0.45, shred_rep 0.30** | 0.738 | 0.778 | **0.758** | 172 | 61 |
| mapq 30, anchor_rep 0.45, shred_rep 0.30 | 0.735 | 0.778 | 0.756 | 172 | 62 |
| mapq 30, anchor_rep 0.45, min_supp 6 | **0.784** | 0.724 | 0.753 | 160 | 44 |
| shipped defaults | 0.685 | 0.787 | 0.733 | 174 | 80 |

0.758 is above the 0.754 recorded for this pair with the old, defective binary -- and reached with
the stored core bytes correct. **Tuned on one pair only.** Nothing is changed in the code until it
holds on h2009 as well; `docs/RESULTS.md`'s own rule is that no default is tuned per pair.

## Is homopolymer compression justified on CASTLE? (measured 2026-08-30)

Source: `tmp/errprof/run.sh` -> `tmp/errprof/out.log`. hcc1954 NORMAL, first 2.41 M reads of
`hifi/hcc1954/n.svb` (~39 Gbp), minimap2 -ax map-hifi to hg38, MAPQ>=20, chr1 only; then
`scripts/hifi_error_profile.py` on chr1:2,000,000-7,000,000 (3,471 reads over the window).
Wall 21:00:27 -> 21:22:59 = 22.5 min, 14 threads.

The script does NOT compare reads to the reference -- at 0.018% error most reference differences are
real variants. It compares reads to each other: a difference carried by <=20% of the reads covering a
locus of depth >=10 is called an error.

```
substitution  25,867  (23.4%)
insertion     36,114  (32.7%)
deletion      48,382  (43.8%)
indel total   84,496  (76.6%)          indel : substitution = 3.27 : 1
```

Homopolymer enrichment, by reference run length, indels against substitutions as the control:

| run | indel | substitution |
|---|---|---|
| 1 | 0.4x | 1.1x |
| 2 | 1.1x | 0.9x |
| 3 | 1.9x | 0.9x |
| 4 | 2.7x | 0.8x |
| 5 | 3.4x | 2.7x |
| 6 | 3.8x | 0.8x |
| 7 | 4.3x | 1.3x |
| **>=3** | **2.3x** | **1.1x** |

Conclusion: the argument holds on real CASTLE reads. Errors are overwhelmingly indels, indel errors
concentrate in homopolymer runs and do so monotonically in run length, and substitutions -- which
have no reason to prefer a run -- show no enrichment, so what is being measured is run-length slip
and not a general "hard region" effect.

**The 38:1 figure in `src/extract.cpp`'s comment does NOT hold here.** It is 3.27:1 on this data.
38:1 came from `scripts/homopolymer_errors.py` on simulated (Badread) reads; that script no longer
exists on this box and the simulation evidently overstates indel dominance against real Revio HiFi.
The related claim "91.3% of insertion errors merely lengthen a homopolymer run" is NOT measured by
this script either -- it measures run CONTEXT, not whether the inserted base extends the run. Neither
number may be used for CASTLE without a new measurement.

Caveats: one sample, one 5 Mb window. The substitution row at run length 5 reads 2.7x against ~1x
everywhere else and rests on 6.1% of the substitution errors -- most likely small numbers.
