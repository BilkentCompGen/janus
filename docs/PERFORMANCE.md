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
