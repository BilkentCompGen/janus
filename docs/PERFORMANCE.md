# Performance: gzip pipeline + lcptools `extract_fast`

The two speed wins on `lcp-cortex-merged` vs the original `match.cpp` assembly pipeline. Both
are baked into `run_lcp` in `src/graph.cpp`; the second also into `src/cores.cpp` for
analysis tools that still go through `extract_cores`.

## 1. Gzip pipeline (FASTQ inflate + parse + dispatch)

### What changed

The original FASTA-input pipeline (`io.cpp::read_fasta_flat` on `match.cpp`-era code) used
line-by-line `gzgets` on a single `gzFile`. That is **CPU-bound on zlib state crossing and
per-line `std::string` ops** — measured ~3-5× slower per reader thread than block reads.

The read pipeline replaces it with:

- **One reader thread per FASTQ shard.** A "shard" is either a file or a `.fastq`/`.fq`
  member of a directory; `collect_shards()` resolves the input. For the standard
  `simulation/reads/{normal,wgs}/hifi.fastq.gz` setup that's two readers, one per sample,
  inflating **in parallel CPUs** (`gzread` is serial within one stream but parallel between
  streams). Reader count = `max(2, threads/4)` so the inflater pool always has ≥2.
- **Block reads, not `gzgets`.** Each reader pulls 1 MB blocks via `gzread`, scans for `\n`
  manually, and emits 4-line FASTQ records into a bounded SPSC queue (`Q_CAP=256`,
  `BATCH=64`). The condition-variable wait avoids busy-spin when workers fall behind.
- **`threads - readers` worker threads** pop batches and run LCP + signature work. Workers
  emit `SigCol{sig, col}` records and `(sig_prev, sig_curr)` edges into thread-local vectors;
  there is **no global hashmap lock on the hot path** (the per-sig count happens later via
  radix-sort + run-length reduce).

### Why parallel reads matter

Single-stream zlib decompresses at ~31 MB/s in / ~195 MB/s out on this hardware. With two
samples each as one `.fastq.gz`, the prior single-reader code was capped at one stream =
31 MB/s in. The new layout runs both readers concurrently → ~62 MB/s combined input,
~390 MB/s combined output. Workers (LCP + murmur) run on the remaining cores and saturate
on the decompressed data.

`pigz`/`libdeflate` do **not** speed up a single plain-gzip stream — gzip streams are
inherently serial. What helps is **more parallel streams**: a 62-shard plain-FASTQ setup
hit 242 MB/s combined input at 16 inflaters (measured prior).

### Measured

| dataset | threads | reader layout | wall (lcp+graph+bubbles) |
|---|---|---|---|
| WG 30× HiFi, 29 GB × 2 (`hifi.fastq.gz`) | 8 | 2 readers + 6 workers | **12:02** |
| WG 30× HiFi, 62 shards × 2 (plain) — prior reference | 32 | 16 readers + 16 workers | 14:12 |
| WG 30× HiFi, original line-by-line `gzgets` — prior baseline | 16 | 1 reader serial | 56:00 |

Note the 12:02 / 8-threads number **beats** the 14:12 / 32-threads number because the
shard-pool overhead and worker idle cycles dominate when there are more workers than
inflate-saturable streams. For a 2-file setup, 6-8 threads is the sweet spot; more is wasted.

## 2. lcptools `extract_fast` SoA API

### What changed

`cores.cpp::extract_cores` (the analysis-tool path) and `graph.cpp::process_read` (the caller
hot path) used to call `init_lps_offset` + `lps_deepen` per record, write into a malloc'd
`struct lps_core *cores`, then copy into our `std::vector<Core>`. Three issues:

1. **Malloc per record.** `init_lps_offset` + `lps_deepen` + `free_lps` allocates fresh
   `bit_rep` and `cores` arrays every call. WG ≈ 12M reads → ~2.4B mallocs/frees → glibc
   arena contention kills scaling past ~11 threads on a 32-core box.
2. **AoS → SoA copy.** `struct lps_core { uint64_t start, end; uint32_t label; }` is 20 B
   AoS; our `Core` is 24 B (padded). The copy is one cache miss per core.
3. **Intermediate vectors.** `extract_cores` → `compute_sigs` writes `Core[]` + `uint64_t[]`
   into vectors that the caller immediately re-iterates. For the read path, both vectors are
   ephemeral.

`extract_fast` (`lib/lcptools` perf-session2 branch) is the SoA fast-path API:

```c
void extract_fast(const char *str, int len, uint64_t offset, int target_level,
                  int *n_out, const uint64_t **starts_out,
                  const uint64_t **ends_out, const uint32_t **labels_out);
```

Returns **three parallel arrays** (`starts[]`, `ends[]`, `labels[]`) from internal
thread-local scratch — **no allocation in the steady state**, no `bit_rep` rebuild per call
(scratch buffers reused across reads via `thread_local struct lps_scratch`).

### Hot-path integration in `process_read`

`process_read` (`src/graph.cpp:43`) inlines the per-read work end-to-end:

```cpp
extract_fast(s.data(), len, 0, level, &n, &st, &en, &lb);
for (int i = 0; i < n; i++) {
    if (st[i] < margin || en[i] + margin > len) continue;       // interior only
    const uint64_t L = st[i] > flank ? st[i] - flank : 0;
    const uint64_t R = en[i] + flank < len ? en[i] + flank : len;
    const uint64_t g = murmur64(s.data() + L, R - L, lb[i]);    // flank-window sig
    lsigs.push_back({g, col});
    if (have && prev != g) ledges.push_back(...);
    prev = g; have = true;
}
```

No `std::vector<Core>` materialized. The murmur seed is the core label, the window is
`[start-flank, end+flank)`. Adjacent within-read edges are collected in the same loop.

### Measured

Per-read LCP time (single-thread, lcptools `lcp_dump_fast` vs `lcp_dump`):

| API | WG 30× per-read mean |
|---|---|
| `init_lps_offset` + `lps_deepen` + `free_lps` (baseline) | 21.1 s / chunk |
| `init_lps_scratch` (scratch reuse only) | ~13 s / chunk |
| `extract_fast` (full SoA fast path) | **8.50 s / chunk** |

Cumulative 60% reduction vs baseline. WG-level the bottleneck shifts back to **gzip
decompression** (the input is gz-bound, not LCP-bound, at 8+ threads).

## 3. Other plumbing wins (smaller individually, add up)

- **LSB radix sort** for the post-pool reduce: `(sig)` 8 passes × 256 buckets, `(a,b)` 8+8
  passes. Faster than `std::sort` on uniform 64-bit keys; ~2× on the WG sig vector.
- **Run-length reduce.** Sorted `SigCol[]` collapses into `SigCount{sig,n,t}` in one linear
  pass; no hashmap.
- **Atomic-rename checkpoints.** Each stage writes `*.tmp` then `rename()`. An interrupted
  write never leaves a half-written file that resume would treat as complete (a previous
  session lost 12 min of LCP this way before it was wired up).
- **`thread_local` reuse** for `extract_fast`'s scratch, `process_read`'s `std::string`
  batch buffer, the worker's edge vector — every persistent allocation is hoisted out of
  the per-read loop.

## 4. Negative results (don't redo these)

### Sharded-counter map (the rust-mdbg DashMap pattern) — 64% slower

Tried 2026-06-02. Replaced the SigCol-emit + 8-pass radix-sort + run-length-reduce design
with a 256-shard `unordered_map<sig, Counts>` (per-shard `std::mutex`). Workers directly
incremented the shared counter per emit (the same pattern rust-mdbg uses with
`Arc<DashMap<Kmer, DbgEntry>>`).

| design | WG 8 threads, 2 `.fq.gz` | CPU % | bubble output |
|---|---|---|---|
| SigCol + radix-sort (current) | **12:02** | 466% | 112,977 |
| 256-shard `unordered_map` + `mutex` | 19:44 | 412% | 112,977 (identical) |

8.5 M voluntary context switches vs 154 K for SigCol = **55× more** worker stalls. User-time
+45%. Reasons it doesn't transfer from Rust:

1. **Emit density.** LCP cores at L5 fire ~1 per 17 bp; rust-mdbg minimizers fire ~1 per
   50–100 bp. We do **~5× more per-emit work**, so per-emit overhead dominates.
2. **`std::unordered_map` ≠ DashMap shard.** std-unordered-map is closed-chaining with
   per-node `malloc`; DashMap is open-addressing flat-table inside its RwLock. 1.1 B
   `m[g]` calls cost ~25 ns extra each vs `vec.push_back`; the radix-sort end-of-stage
   merge cost is **less** than that aggregate.
3. **Mutex is heavier than RwLock-write.** `std::mutex` lock acquire under contention is a
   futex syscall on Linux → context switch → reschedule. Rust's parking_lot RwLock
   (DashMap default) is a spin-then-park hybrid.

SigCol's "delayed shard-and-count" (radix sort partitions by `sig & 0xFF` per pass = the
hash-shard semantics of DashMap, just batched once at the end) is empirically the right
trade for our dense emit rate. Don't revisit without ALSO swapping `std::unordered_map`
for a flat open-addressing map (e.g. `ankerl::unordered_dense` or `absl::flat_hash_map`)
AND switching to a spin-mutex.

### Pre-shard input → 32 threads / 62 plain shards — only 21% faster than 16t baseline

Tried 2026-06-02. Pre-decompressed both `hifi.fastq.gz` into 31 record-aligned plain
FASTQ shards each (`pigz -dc | split -l 800000`), then ran svp at 32 threads.

| inputs | threads | wall | CPU % |
|---|---|---|---|
| 1 `.fq.gz` × 2 | 8 | 12:02 | 466% |
| 1 `.fq.gz` × 2 | 16 | 9:52 | 622% |
| 62 plain shards | 32 | 7:47 | **895%** |

895% CPU = ~9 effective cores out of 32 → workers still stall (~14 of 24 idle on avg).
With inflate uncapped the new bottleneck is somewhere else (edge-merge phase, allocator
contention on `lsigs.push_back`'s growth, or LCP saturation per malloc tree). Sharding
also costs 8:59 one-time decompress + 346 GB disk. Net benefit too small for the cost.

Conclusion: **on this workload, the brute-force parallelism path plateaus around 8 cores**;
beyond that, algorithmic changes are required (e.g. minimizer-space sketching like
rust-mdbg, or reverse-complement signature canonicalization to halve the work upstream).

## 5. Per-phase timing (2026-06-02, WG 8 threads, lcp-cortex-merged)

Run was instrumented with `Tic` wall-clock marks at every stage boundary inside `run_lcp` /
`run_graph` / `run_bubbles`, plus aggregate per-reader / per-worker accounting in `run_pool`
(busy ns, idle ns, queue-full-wait ns). Same inputs as §1: two `hifi.fastq.gz`, 29 GB each.

| phase | wall | % of total | threads doing work |
|---|---:|---:|---|
| **lcp**: reader+worker pool | **444.86 s** | 61.4% | **7.56 / 8 (94.5%)** |
| **lcp**: radix-sort sigs (2.59 B records, 41.5 GB) | 63.58 s | 8.8% | 1 |
| **lcp**: reduce sigs → SigCount | 9.25 s | 1.3% | 1 |
| **lcp**: radix-sort edges (2.58 B records, 41.3 GB) | 127.77 s | 17.6% | 1 |
| **lcp**: dedup edges | 2.79 s | 0.4% | 1 |
| **lcp**: write counts.bin + edges.bin | 2.35 s | 0.3% | 1 |
| **graph**: read counts.bin | 1.83 s | 0.3% | 1 |
| **graph**: prune+colour | 12.34 s | 1.7% | 1 |
| **graph**: read edges.bin | 2.49 s | 0.3% | 1 |
| **graph**: resolve edges (sig→id hashmap lookup) | 27.02 s | 3.7% | 1 |
| **graph**: sort+dedup id-edges | 1.10 s | 0.2% | 1 |
| **graph**: write graph.bin | 0.19 s | 0.0% | 1 |
| **bubbles**: read graph.bin | 0.18 s | 0.0% | 1 |
| **bubbles**: build adjacency list | 13.76 s | 1.9% | 1 |
| **bubbles**: BFS connected-components | 0.19 s | 0.0% | 1 |
| **bubbles**: write bubbles.tsv | 0.01 s | 0.0% | 1 |
| **total wall** | **725 s = 12:05** | 100% | 467% avg = 4.67 effective cores |

### What the 467% average actually hides

Inside the **444 s pool phase** the picture is good:

- Workers (`process_read`: LCP + sig + edge collect) sum to 2516.6 s across 6 workers — that's
  **5.67 effective cores**, ~100% per-worker utilization (worker-idle = 0.1 s total).
- Readers wall-time sum is 838.7 s across 2 readers = **1.89 effective cores** — 95% per
  reader. **266 s of that was `queue-full` wait** (32% of reader wall): readers can produce
  faster than workers consume. So inflate is *not* the bottleneck; LCP work is.
- Pool effective cores = 5.67 + 1.89 = **7.56 / 8 = 94.5% utilization**.

Everything after the pool (~287 s) is **single-threaded**:

- 192 s of radix-sort on the 41 GB SigCol and 41 GB SEdge vectors (memory-bandwidth-bound,
  8 sequential passes × 41 GB read+write each ≈ 64 s/pass × 8 = 512 GB of memory traffic).
- 14 s reduce / dedup / write.
- 45 s graph stage (resolve_edges sig→id hashmap is the bulk).
- 14 s bubbles stage (build_adj is the bulk).

Averaged: `(444 s × 7.56 cores + 287 s × 1 core) / 725 s wall = 5.02 effective cores`,
matching the observed 467%. **The "low CPU%" is entirely the 287 s of serial tail.**

### What this implies for next round

The merge phase is the clearest target — 287 s × 7 idle cores = ~2000 core-seconds left on
the floor. Two ways to attack it (deferred until after E2E call accuracy is sorted):

1. **Parallel radix sort.** Each pass partitions chunks across threads (count-then-scatter
   with per-thread offsets). Memory bandwidth is shared so it won't be 8×, but ~3-4×
   reduction on the 192 s sort phase is realistic → cuts 287 s → ~150 s, total 12:05 → ~9:35.
2. **Bypass the sort.** Workers emit to **256 thread-local partition buffers** (`sig & 0xFF`)
   instead of one flat vector. At end-of-pool, 8 reducer threads each handle 32 shards,
   producing the final `SigCount`s in parallel. No global sort, no 41 GB monolithic
   vector → also cuts peak RAM. The "real" algorithmic fix and a partial reprise of the
   shard-map idea, but with **no per-emit locks** (each worker has its own 256 vectors).

Same shape for the edge merge.

## How to reproduce

```
make                       # builds ./svp on lcp-cortex-merged
/usr/bin/time -v ./svp run simulation/reads/normal/hifi.fastq.gz \
                          simulation/reads/wgs/hifi.fastq.gz \
                          outputs/svp/wgs --threads 8

make tools/bubble_recall
/usr/bin/time -v ./tools/bubble_recall \
    outputs/svp/wgs/counts.bin outputs/svp/wgs/graph.bin \
    simulation/ref.fasta simulation/truth/wgs/out.fa \
    --level 5 --flank 10 --min-mult 3
```

WG run measured 2026-06-01 on AMD Ryzen 9 9950X, 157 GB RAM, NVMe; reproducible without
external network (all inputs in `simulation/`, all deps under `lib/` after
`make install-deps`).
