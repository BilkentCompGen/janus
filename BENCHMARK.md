# BENCHMARK — `graph` build

Performance record for the `graph` stage (`run_graph_reads`: count → fill → resolve → write).

## Ground rules (read before measuring)

- **16-core budget. ALWAYS `--threads 16` + `taskset -c 0-15`.** The box exposes 32 logical CPUs but only 16 are ours. Running `--threads 32` breaks the budget, makes wall-times incomparable to the baseline, AND oversubscribes (workers + 2 readers + rapidgzip decode threads > 16) which inflates even the contention-immune CPU-seconds ~1.45× via cache/scheduling contention. `taskset` is required so the rapidgzip child processes inherit the cap.
- **Data:** `hard` WG sim, `min_mult=10`. normal + tumor HiFi, each 77 GB `.gz` → 193.5 GB FASTQ (~96 GB sequence/file). 12,915,155 reads total; ~5.4 G core emissions.
- **How CPU-seconds are measured:** per-thread `CLOCK_THREAD_CPUTIME_ID` (advances only on-CPU → contention-immune, excludes idle/blocked) for readers/workers; `getrusage(RUSAGE_CHILDREN)` for the rapidgzip decode processes; `RUSAGE_SELF` per-phase deltas cross-check the total. `efficiency = CPU / (16 · wall)` = fraction of the 16-core budget kept busy.
- All numbers below: idle box, pinned to `0-15`, measured **2026-06-17**.

## Head-to-head (same box, 16 cores, pinned, threads=16)

| phase | HEAD `eb0194e` (zlib + atomic bloom + re-LCP) | current (rapidgzip + plain bloom + boundary cache) | Δ |
|---|--:|--:|--:|
| count (pass 1) | 669.0 s | 512.6 s | −23.4% |
| fill (pass 2)  | 668.3 s | 501.9 s | −24.9% |
| colour+resolve+write | 64.5 s | 82.3 s | +18 s |
| **TOTAL** | **1401.8 s = 23.4 min** | **1096.8 s = 18.3 min** | **−305 s (−21.8%)** |

Output equivalent (bloom is nondeterministic within a small band): HEAD 86,442,101 survivors / 88,023,070 edges vs current 86,442,239 / 88,023,125; both 14.9 GB payload.

The committed baseline (CLAUDE.md, 2026-06-15) records 24 min @ 16 threads (pass1 685 s + pass2 683 s = 1368 s) — the same-box HEAD run above (1402 s) confirms it within 2.5%, i.e. **no regression in committed code**.

### What changed this session (all uncommitted working-tree)
- **rapidgzip in-process parallel decode** (was zlib `gzread`, reader-bound) → pass 1 −23%.
- **plain (non-atomic) bloom** in the count cascade (btllib's per-bit atomic `fetch_or` was redundant under the per-shard mutex).
- **core-boundary cache**: pass 1 records each read's LCP core `(st,width)` (uint16, ~21.8 GB) so pass 2 replays instead of re-running `extract_fast` → fill `extract_fast` 3657→**0** core-s.

## CPU-seconds breakdown — current, 16 cores (total 12,971 core-s, wall 1096.8 s, eff 73.9%)

| phase | wall | CPU-s | in-proc | decode | eff |
|---|--:|--:|--:|--:|--:|
| count (`pool_both`) | 512.6 s | 7313.8 | 6475.1 | 838.7 | 89.2% |
| colour + idx (serial) | 15.0 s | 14.9 | 14.9 | 0 | — |
| fill (`fill_pool`) | 501.9 s | 5575.3 | 4739.4 | 836.0 | 69.4% |
| resolve_edges (serial) | 33.1 s | 32.9 | 32.9 | 0 | — |
| seqlen + write (serial) | 34.2 s | 34.0 | 34.0 | 0 | — |
| **TOTAL** | **1096.8 s** | **12971.1** | 11296.4 | 1674.7 | **73.9%** |

### By component (sub-phase thread-CPU split, both passes rolled up)

| component | core-s | % of build | notes |
|---|--:|--:|---|
| **LCP `extract_fast`** | **3656.7** | 28.2% | count pass only (fill = 0, cached) |
| **per-read fill**: `idx.find` + slot-write + `seen`-dedup + murmur | **4225.7** | 32.6% | fill pass; `idx.find` on the 86.4 M-entry map dominates |
| **counting**: bloom cascade + survivor map | 2359.4 | 18.2% | count pass |
| **gzip decode** (rapidgzip) | 1674.7 | 12.9% | 838.7 count + 836.0 fill (both passes re-decode) |
| edge insert + radix-compact | 297.5 | 2.3% | fill pass |
| reader FASTQ-parse | 383.7 | 3.0% | 172.1 count + 211.6 fill |
| sig/murmur (count) | 244.3 | 1.9% | |
| **serial tail** (merge + colour + resolve + write) | 110.2 | 0.8% | single-threaded |

Count pass is well-parallelized (89% of 16 cores busy). Fill is feed-bound (69%): once `extract_fast` is cached the workers outrun the 2-reader + 4-decode feed, which can't deliver fast enough on a 16-core budget.

## Theoretical optimum — holding LCP + gzip decode fixed

The two costs we've decided **not** to optimize further are **LCP `extract_fast` (3656.7 core-s)** and **gzip decode (1674.7 core-s)** = **5331.4 core-s irreducible**. Everything else (counting, per-read fill, edges, parse, serial tail) is in principle optimizable toward zero. The hard floor is therefore that irreducible CPU spread over the core budget:

**floor(N) = (LCP + decode) / N = 5331.4 / N** (assumes perfect parallelism and all other work → 0):

| N cores | floor |
|--:|--:|
| 1  | 5331 s (88.9 min) |
| 2  | 2666 s (44.4 min) |
| 4  | 1333 s (22.2 min) |
| 8  |  666 s (11.1 min) |
| 12 |  444 s (7.4 min) |
| **16** | **333 s (5.6 min)** |

### Layered floors at 16 cores (the realistic ladder)

| boundary | wall @ 16c | what it assumes |
|---|--:|---|
| **absolute floor** (LCP + decode only, rest → 0) | **5.6 min** | the only irreducible work is LCP + decode |
| Amdahl floor (all *current* CPU, perfect parallel, serial tail = 110 core-s stays serial) | **15.2 min** | `110 + (12971−110)/16` |
| **current actual** | **18.3 min** | measured |
| old HEAD (zlib + re-LCP) | 23.4 min | measured |

Reading the ladder:
- **18.3 → 15.2 min** (3.1 min): lost to imperfect parallelism — the fill-pass feed coupling (decode↔reader↔queue) and the small serial tail not being hidden. A reader/feed restructure is the lever (not core-count tuning — that was tried and gave <3%).
- **15.2 → 5.6 min** (9.6 min): the "other work" that LCP+decode-fixed allows optimizing — dominated by **per-read fill 4226 core-s** (the `idx.find` on the 86.4 M map: ~5.4 B cache-missing lookups → flat open-addressing table + software prefetch, or a minimal perfect hash) and **counting 2359 core-s**.

### Feed isolation (informs the fill floor)

rapidgzip decode of one 77 GB file (→193.5 GB): **-P4 = 128 s (1514 MB/s)**, -P16 = 99 s (doesn't scale past ~4 threads — memory/IO-bound). Single-thread FASTQ parse ≈ 132 s/file. So the feed's parallel floor is ~132 s; the fill pass running 502 s means it's coupling/IO-bound, not compute-bound, on the current pipeline.
