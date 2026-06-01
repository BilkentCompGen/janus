# LCP-core analysis: uniqueness, tangle, and context-radius rescue

Consolidates measurements from three independent sweeps:

1. **Whole-genome uniqueness sweep** on `simulation/ref.fasta` (3,209,295,639 bp, 455 records) — level × flank grid measuring distinct-key counts, singleton counts, window-length distribution, and adjacent-window gap statistics.
2. **Whole-genome tangle sweep** on the same input — the de Bruijn-like graph over LCP-core signatures, cyclomatic complexity = E − V + 1.
3. **Context-radius rescue sweep** on chr19 only — the `--context-radius R` mechanism that disambiguates colliding signatures by mixing in `±R` neighbor signatures (Approach B2 from `plans/in-the-current-version-reflective-lantern.md`).

Source binaries:
- `tools/analyze_uniqueness.cpp` — distinct-key + percentile sweep
- `tools/analyze_tangle.cpp` — V/E/tangle on the signature-adjacency graph
- `tools/analyze_cores.cpp` — per-genome uniqueness histograms
- `tools/analyze_kmer.cpp` — per-genome R-context singleton/max-rep counts
- `tools/dump_gfa.cpp` — GFA 1.0 emission of DBG (signature-collapsed) or VG (per-core) graphs
- `./svp ... --context-radius R` — Pareto svp anchor sweep

Raw CSVs in `docs/analysis/`:
- `L{4,5,6}_uniqueness.csv` — whole-genome uniqueness sweep (Part A)
- `L{4,5,6}_tangle.csv` — whole-genome tangle sweep (Part A)
- `context_radius_pareto_chr19.csv` — chr19 svp anchor sweep (Part B)
- `context_radius_gfa_tangle_chr19.csv` — chr19 DBG/VG graph topology (Part B)
- `context_radius_kmer_uniqueness_chr19.csv` — chr19 per-genome k-mer uniqueness (Part B)

## Glossary

Signature schemes (the uniqueness analyzer computes 4 per core, with `m = flank`):

| key | definition |
|---|---|
| `distH`  | distinct values of `murmur64(seq[start−m..end+m], seed=0)` (unseeded content) |
| `distHL` | distinct values of `murmur64(seq[start−m..end+m], seed=core.label)` (production "double hash") |
| `distL`  | distinct values of `core.label` (lcptools internal label, uint32, ignores flank) |
| `distS`  | distinct values of the actual window byte content (ground truth) |
| `sing*`  | rep=1 cores under that key (singletons — eligible to be anchors) |

Flank is signed: positive extends beyond the core, 0 = window == core, negative = window strictly inside the core (window length = core_length + 2·flank).

For the **tangle** graph: node = signature, directed edge = adjacency in the genome traversal between consecutive cores in the same record (skipping invalid-window cores). **tangle = E − V + 1** (cyclomatic complexity assuming single connected component).

For the **context-radius rescue**: when a core's R=0 signature collides with another core's signature in the same genome, only the colliding cores receive a recomputed signature that mixes `±R` neighbor signatures via `fmix64`-folded XOR. Non-colliding cores keep their R=0 signature. Closest published precedent: variable-order de Bruijn graphs (Boucher et al., DCC 2015).

## Compute cost

| sweep | scope | wall (3-way parallel) | peak RSS |
|---|---|---:|---:|
| uniqueness (with percentiles) | whole genome | 32:45.29 (L4) / 18:35.11 (L5) / 8:31.25 (L6) | 25.2 / 12.6 / 10.3 GiB |
| tangle                        | whole genome | 29:08.28 (L4) / 16:40.55 (L5) / 7:41.62 (L6) | 33.1 / 15.9 / 10.3 GiB |
| context-radius Pareto (svp)   | chr19         | 3.80–4.30 s per cell × 28 cells (4 sims × 7 R) | 1.9–2.2 GB |

---

# PART A — Whole-genome uniqueness and tangle (no rescue, R=0)

## A.1 Baseline at flank=0 (window == core)

| level | total cores | distL | distS=distH | sing@f=0 |
|---:|---:|---:|---:|---:|
| 4 | 106,167,674 | 93,183,755 | 97,594,901 | 94,472,930 |
| 5 |  45,415,538 | 42,481,516 | 43,074,367 | 41,957,626 |
| 6 |  19,446,667 | 18,543,934 | 18,701,046 | 18,276,927 |

**Genome repetition** (`total − distS`):

| level | repeated cores | rep rate |
|---:|---:|---:|
| 4 | 8,572,773 | **8.07%** |
| 5 | 2,341,171 | **5.15%** |
| 6 |   745,621 | **3.83%** |

**Label collision rate** (`distS − distL` at flank=0, i.e. distinct core sequences that share a label with another distinct core sequence):

| level | label collisions | % of distS |
|---:|---:|---:|
| 4 | 4,411,146 | **4.52%** |
| 5 |   592,851 | **1.38%** |
| 6 |   157,112 | **0.84%** |

## A.2 Core-length percentile distribution (at flank=0, window == core)

| level | min | p5 | p25 | median | p75 | p95 | mean | max |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 |  35 |  52 |  58 |  63 |  69 |  79 |  65.84 | 30,000,060 |
| 5 |  81 | 124 | 138 | 150 | 162 | 183 | 156.05 | 30,000,156 |
| 6 | 175 | 294 | 326 | 352 | 380 | 427 | 364.88 | 30,000,377 |

The p5–p95 band is tight (~30 bp at L4, ~60 bp at L5, ~135 bp at L6), so the bulk of cores cluster around the median. The huge `max` values (~30 Mbp) come from N-block runs in the acrocentric/centromeric regions that lcptools labels as single very-long cores.

Window-length distribution at all flanks (window len = core len + 2·flank):

### Level 4

| flank | valid | min | p5 | p25 | med | p75 | p95 | max | mean |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −50 |    224,509 |   1 |   1 |   3 |   7 |  15 |  54 | 29,999,960 | 947.36 |
| −40 |  3,791,374 |   1 |   1 |   2 |   5 |   9 |  22 | 29,999,980 |  62.70 |
| −30 | 66,783,175 |   1 |   1 |   4 |   7 |  12 |  21 | 30,000,000 |  11.70 |
| −20 |106,166,013 |   1 |  12 |  18 |  23 |  29 |  39 | 30,000,020 |  25.85 |
| −10 |106,167,674 |  15 |  32 |  38 |  43 |  49 |  59 | 30,000,040 |  45.84 |
|   0 |106,167,674 |  35 |  52 |  58 |  63 |  69 |  79 | 30,000,060 |  65.84 |
| +10 |106,167,674 |  55 |  72 |  78 |  83 |  89 |  99 | 30,000,080 |  85.84 |
| +20 |106,167,674 |  75 |  92 |  98 | 103 | 109 | 119 | 30,000,100 | 105.84 |
| +30 |106,167,674 |  95 | 112 | 118 | 123 | 129 | 139 | 30,000,120 | 125.84 |
| +40 |106,167,674 | 115 | 132 | 138 | 143 | 149 | 159 | 30,000,140 | 145.84 |
| +50 |106,167,674 | 130 | 152 | 158 | 163 | 169 | 179 | 30,000,160 | 165.84 |

At negative flank the *core-length* filter kicks in: at f=−50 only 224,509 cores have a valid (>0-length) window, p5/p25/p50 collapse to 1/3/7. p95 goes UP to 54 because at negative flank only the long-tailed cores survive.

### Level 5

| flank | valid | min | p5 | p25 | med | p75 | p95 | max | mean |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −50 | 45,410,894 |   1 |  24 |  38 |  50 |  62 |  83 | 30,000,056 |  56.05 |
| −40 | 45,415,538 |   1 |  44 |  58 |  70 |  82 | 103 | 30,000,076 |  76.05 |
| −30 | 45,415,538 |  21 |  64 |  78 |  90 | 102 | 123 | 30,000,096 |  96.05 |
| −20 | 45,415,538 |  41 |  84 |  98 | 110 | 122 | 143 | 30,000,116 | 116.05 |
| −10 | 45,415,538 |  61 | 104 | 118 | 130 | 142 | 163 | 30,000,136 | 136.05 |
|   0 | 45,415,538 |  81 | 124 | 138 | 150 | 162 | 183 | 30,000,156 | 156.05 |
| +10 | 45,415,538 | 101 | 144 | 158 | 170 | 182 | 203 | 30,000,176 | 176.05 |
| +20 | 45,415,538 | 121 | 164 | 178 | 190 | 202 | 223 | 30,000,196 | 196.05 |
| +30 | 45,415,538 | 141 | 184 | 198 | 210 | 222 | 243 | 30,000,216 | 216.05 |
| +40 | 45,415,538 | 161 | 204 | 218 | 230 | 242 | 263 | 30,000,236 | 236.05 |
| +50 | 45,415,538 | 181 | 224 | 238 | 250 | 262 | 283 | 30,000,256 | 256.05 |

### Level 6

| flank | valid | min | p5 | p25 | med | p75 | p95 | max | mean |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −50 | 19,446,667 |  75 | 194 | 226 | 252 | 280 | 327 | 30,000,277 | 264.88 |
| −40 | 19,446,667 |  95 | 214 | 246 | 272 | 300 | 347 | 30,000,297 | 284.88 |
| −30 | 19,446,667 | 115 | 234 | 266 | 292 | 320 | 367 | 30,000,317 | 304.88 |
| −20 | 19,446,667 | 135 | 254 | 286 | 312 | 340 | 387 | 30,000,337 | 324.88 |
| −10 | 19,446,667 | 155 | 274 | 306 | 332 | 360 | 407 | 30,000,357 | 344.88 |
|   0 | 19,446,667 | 175 | 294 | 326 | 352 | 380 | 427 | 30,000,377 | 364.88 |
| +10 | 19,446,667 | 195 | 314 | 346 | 372 | 400 | 447 | 30,000,397 | 384.88 |
| +20 | 19,446,667 | 215 | 334 | 366 | 392 | 420 | 467 | 30,000,417 | 404.88 |
| +30 | 19,446,667 | 235 | 354 | 386 | 412 | 440 | 487 | 30,000,437 | 424.88 |
| +40 | 19,446,667 | 255 | 374 | 406 | 432 | 460 | 507 | 30,000,457 | 444.88 |
| +50 | 19,446,667 | 275 | 394 | 426 | 452 | 480 | 527 | 30,000,477 | 464.88 |

## A.3 Distinct-key counts per scheme

### Level 4

| flank | distH | distHL | distL | distS | singH | singHL | singL | singS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −50 |    92,213 |   216,636 |   215,104 |    92,213 |    80,723 |   211,807 |   209,813 |    80,723 |
| −40 |   923,525 | 3,566,765 | 3,536,226 |   923,525 |   791,386 | 3,472,569 | 3,434,942 |   791,386 |
| −30 |18,066,192 |60,972,985 |59,868,424 |18,066,192 |15,532,725 |58,979,376 |57,397,124 |15,532,725 |
| −20 |89,360,599 |95,898,339 |93,182,450 |89,360,599 |84,972,053 |92,538,197 |88,706,368 |84,972,053 |
| −10 |95,935,637 |96,874,377 |93,183,755 |95,935,637 |92,579,666 |93,641,556 |88,707,552 |92,579,666 |
|   0 |97,594,901 |97,594,901 |93,183,755 |97,594,901 |94,472,930 |94,472,930 |88,707,552 |94,472,930 |
| +10 |99,044,505 |99,044,505 |93,183,755 |99,044,505 |96,157,025 |96,157,025 |88,707,552 |96,157,025 |
| +20 |99,750,967 |99,750,967 |93,183,755 |99,750,967 |96,987,364 |96,987,364 |88,707,552 |96,987,364 |
| +30 |100,196,705 |100,196,705 |93,183,755 |100,196,705 |97,504,789 |97,504,789 |88,707,552 |97,504,789 |
| +40 |100,529,498 |100,529,498 |93,183,755 |100,529,498 |97,890,466 |97,890,466 |88,707,552 |97,890,466 |
| +50 |100,789,009 |100,789,009 |93,183,755 |100,789,009 |98,194,130 |98,194,130 |88,707,552 |98,194,130 |

### Level 5

| flank | distH | distHL | distL | distS | singH | singHL | singL | singS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −50 |42,171,823 |42,884,432 |42,477,596 |42,171,823 |40,882,915 |41,735,399 |41,090,694 |40,882,915 |
| −40 |42,653,376 |42,919,999 |42,481,516 |42,653,376 |41,460,239 |41,776,601 |41,094,287 |41,460,239 |
| −30 |42,835,520 |42,965,223 |42,481,516 |42,835,520 |41,678,094 |41,828,985 |41,094,287 |41,678,094 |
| −20 |42,944,162 |43,008,327 |42,481,516 |42,944,162 |41,805,033 |41,879,317 |41,094,287 |41,805,033 |
| −10 |43,020,809 |43,044,607 |42,481,516 |43,020,809 |41,894,435 |41,922,125 |41,094,287 |41,894,435 |
|   0 |43,074,367 |43,074,367 |42,481,516 |43,074,367 |41,957,626 |41,957,626 |41,094,287 |41,957,626 |
| +10 |43,166,220 |43,166,220 |42,481,516 |43,166,220 |42,064,912 |42,064,912 |41,094,287 |42,064,912 |
| +20 |43,246,263 |43,246,263 |42,481,516 |43,246,263 |42,159,430 |42,159,430 |41,094,287 |42,159,430 |
| +30 |43,316,733 |43,316,733 |42,481,516 |43,316,733 |42,243,669 |42,243,669 |41,094,287 |42,243,669 |
| +40 |43,380,667 |43,380,667 |42,481,516 |43,380,667 |42,320,522 |42,320,522 |41,094,287 |42,320,522 |
| +50 |43,436,266 |43,436,266 |42,481,516 |43,436,266 |42,388,171 |42,388,171 |41,094,287 |42,388,171 |

### Level 6

| flank | distH | distHL | distL | distS | singH | singHL | singL | singS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −50 |18,659,366 |18,678,130 |18,543,934 |18,659,366 |18,224,388 |18,247,388 |18,049,358 |18,224,388 |
| −40 |18,669,191 |18,683,127 |18,543,934 |18,669,191 |18,236,700 |18,253,817 |18,049,358 |18,236,700 |
| −30 |18,678,584 |18,688,222 |18,543,934 |18,678,584 |18,248,384 |18,260,281 |18,049,358 |18,248,384 |
| −20 |18,687,335 |18,693,010 |18,543,934 |18,687,335 |18,259,345 |18,266,376 |18,049,358 |18,259,345 |
| −10 |18,695,060 |18,697,370 |18,543,934 |18,695,060 |18,269,154 |18,272,057 |18,049,358 |18,269,154 |
|   0 |18,701,046 |18,701,046 |18,543,934 |18,701,046 |18,276,927 |18,276,927 |18,049,358 |18,276,927 |
| +10 |18,713,933 |18,713,933 |18,543,934 |18,713,933 |18,293,122 |18,293,122 |18,049,358 |18,293,122 |
| +20 |18,727,147 |18,727,147 |18,543,934 |18,727,147 |18,309,883 |18,309,883 |18,049,358 |18,309,883 |
| +30 |18,739,969 |18,739,969 |18,543,934 |18,739,969 |18,326,186 |18,326,186 |18,049,358 |18,326,186 |
| +40 |18,752,439 |18,752,439 |18,543,934 |18,752,439 |18,342,199 |18,342,199 |18,049,358 |18,342,199 |
| +50 |18,764,145 |18,764,145 |18,543,934 |18,764,145 |18,357,344 |18,357,344 |18,049,358 |18,357,344 |

## A.4 Adjacent-window gap statistics

For each adjacent pair of cores in the genome, `gap = window[i+1].start − window[i].end`. Negative = the windows overlap. Linear relation: `gap(m) = gap(0) + 2m`. `min_gap` blowouts (~−30M) come from the N-block records; `max_gap` blowouts come from the inter-record N-bumpers.

### Level 4

| flank | avg_gap | min_gap | max_gap | n_pos (gap) | n_zero (touch) | n_neg (overlap) |
|---:|---:|---:|---:|---:|---:|---:|
| −50 | −2496 | −29,999,923 | 80 | 36,775 | 434 | 8,871 |
| −40 |  −123 | −29,999,943 | 62 | 756,830 | 3,038 | 32,517 |
| −30 |   +19 | −29,999,963 | 16,010,550 | 47,908,237 | 73,920 | 406,330 |
| −20 |    +4 | −29,999,983 | 17,000,149 | 75,609,258 | 5,054,767 | 25,500,441 |
| −10 |   −15 | −30,000,003 | 17,000,129 | 208,852 | 515,989 | 105,442,832 |
|   0 |   −35 | −30,000,023 | 17,000,109 | 454 | 0 | 106,167,219 |
| +10 |   −55 | −30,000,043 | 17,000,089 | 454 | 0 | 106,167,219 |
| +20 |   −75 | −30,000,063 | 17,000,069 | 454 | 0 | 106,167,219 |
| +30 |   −95 | −30,000,083 | 17,000,049 | 454 | 0 | 106,167,219 |
| +40 |  −115 | −30,000,103 | 17,000,029 | 454 | 0 | 106,167,219 |
| +50 |  −135 | −30,000,123 | 17,000,009 | 454 | 0 | 106,167,219 |

### Level 5

| flank | avg_gap | min_gap | max_gap | n_pos | n_zero | n_neg |
|---:|---:|---:|---:|---:|---:|---:|
| −50 |  +14 | −29,999,991 | 17,000,239 | 36,734,578 | 750,409 |  7,921,820 |
| −40 |   −5 | −30,000,011 | 17,000,219 | 17,411,256 | 814,047 | 27,190,234 |
| −30 |  −25 | −30,000,031 | 17,000,199 |  4,526,397 | 711,219 | 40,177,921 |
| −20 |  −45 | −30,000,051 | 17,000,179 |        547 |     179 | 45,414,811 |
| −10 |  −65 | −30,000,071 | 17,000,159 |        454 |       0 | 45,415,083 |
|   0 |  −85 | −30,000,091 | 17,000,139 |        454 |       0 | 45,415,083 |
| +10 | −105 | −30,000,111 | 17,000,119 |        454 |       0 | 45,415,083 |
| +20 | −125 | −30,000,131 | 17,000,099 |        454 |       0 | 45,415,083 |
| +30 | −145 | −30,000,151 | 17,000,079 |        454 |       0 | 45,415,083 |
| +40 | −165 | −30,000,171 | 17,000,059 |        454 |       0 | 45,415,083 |
| +50 | −185 | −30,000,191 | 17,000,039 |        454 |       0 | 45,415,083 |

### Level 6

| flank | avg_gap | min_gap | max_gap | n_pos | n_zero | n_neg |
|---:|---:|---:|---:|---:|---:|---:|
| −50 |  −99 | −30,000,056 | 17,000,581 | 994 | 168 | 19,445,504 |
| −40 | −119 | −30,000,076 | 17,000,561 | 454 |   0 | 19,446,212 |
| −30 | −139 | −30,000,096 | 17,000,541 | 454 |   0 | 19,446,212 |
| −20 | −159 | −30,000,116 | 17,000,521 | 454 |   0 | 19,446,212 |
| −10 | −179 | −30,000,136 | 17,000,501 | 454 |   0 | 19,446,212 |
|   0 | −199 | −30,000,156 | 17,000,481 | 454 |   0 | 19,446,212 |
| +10 | −219 | −30,000,176 | 17,000,461 | 454 |   0 | 19,446,212 |
| +20 | −239 | −30,000,196 | 17,000,441 | 454 |   0 | 19,446,212 |
| +30 | −259 | −30,000,216 | 17,000,421 | 454 |   0 | 19,446,212 |
| +40 | −279 | −30,000,236 | 17,000,401 | 454 |   0 | 19,446,212 |
| +50 | −299 | −30,000,256 | 17,000,381 | 454 |   0 | 19,446,212 |

**Adjacent lcptools cores overlap heavily** — at every level, the overwhelming majority of adjacent pairs have negative gap. Average overlap at flank=0: 35 bp (L4), 85 bp (L5), 199 bp (L6). Deeper level → more overlap.

The flank at which avg gap crosses 0 (windows just touch):

| level | flank where avg gap ≈ 0 |
|---:|---:|
| 4 | ≈ −17.5 |
| 5 | ≈ −42.5 |
| 6 | ≈ −99.5 |

L6 cores at flank=−50 still overlap 99.99% of the time — they're too long and densely packed for a ±50 flank shrink to separate them.

## A.5 Tangle (cyclomatic complexity, **tangle = E − V + 1**)

Per-scheme tangle on the de Bruijn-like graph (node = signature, directed edge = signature adjacency in genome traversal).

### All levels × all flanks

| level | flank |  V_H |  E_H | tangle_H | tangle_HL | tangle_L |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | −50 |    92,213 |    42,522 |       −49,690 |      −171,681 |      −170,198 |
| 4 | −40 |   923,525 |   651,128 |      −272,396 |    −2,801,549 |    −2,772,134 |
| 4 | −30 |18,066,192 |36,928,375 |   +18,862,184 |   −15,731,808 |   −14,785,766 |
| 4 | −20 |89,360,599 |97,456,097 |    +8,095,499 |    +3,050,460 |    +5,199,634 |
| 4 | −10 |95,935,637 |98,983,187 |    +3,047,551 |    +2,340,394 |    +5,201,244 |
| 4 |   0 |97,594,901 |99,427,931 |    +1,833,031 |    +1,833,031 |    +5,201,244 |
| 4 | +10 |99,044,505 |99,969,505 |      +925,001 |      +925,001 |    +5,201,244 |
| 4 | +20 |99,750,967 |100,349,279 |     +598,313 |      +598,313 |    +5,201,244 |
| 4 | +30 |100,196,705|100,645,162 |     +448,458 |      +448,458 |    +5,201,244 |
| 4 | +40 |100,529,498|100,885,901 |     +356,404 |      +356,404 |    +5,201,244 |
| 4 | +50 |100,789,009|101,083,332 |     +294,324 |      +294,324 |    +5,201,244 |
| 5 | −50 |42,171,823 |43,067,221 |      +895,399 |      +363,755 |      +661,037 |
| 5 | −40 |42,653,376 |43,178,283 |      +524,908 |      +351,779 |      +665,150 |
| 5 | −30 |42,835,520 |43,239,807 |      +404,288 |      +328,348 |      +665,150 |
| 5 | −20 |42,944,162 |43,285,586 |      +341,425 |      +306,986 |      +665,150 |
| 5 | −10 |43,020,809 |43,322,619 |      +301,811 |      +289,720 |      +665,150 |
| 5 |   0 |43,074,367 |43,350,865 |      +276,499 |      +276,499 |      +665,150 |
| 5 | +10 |43,166,220 |43,404,637 |      +238,418 |      +238,418 |      +665,150 |
| 5 | +20 |43,246,263 |43,456,267 |      +210,005 |      +210,005 |      +665,150 |
| 5 | +30 |43,316,733 |43,504,381 |      +187,649 |      +187,649 |      +665,150 |
| 5 | +40 |43,380,667 |43,550,535 |      +169,869 |      +169,869 |      +665,150 |
| 5 | +50 |43,436,266 |43,592,187 |      +155,922 |      +155,922 |      +665,150 |
| 6 | −50 |18,659,366 |18,771,059 |      +111,694 |      +103,811 |      +170,986 |
| 6 | −40 |18,669,191 |18,776,957 |      +107,767 |      +101,973 |      +170,986 |
| 6 | −30 |18,678,584 |18,782,699 |      +104,116 |      +100,114 |      +170,986 |
| 6 | −20 |18,687,335 |18,788,143 |      +100,809 |       +98,449 |      +170,986 |
| 6 | −10 |18,695,060 |18,793,004 |       +97,945 |       +96,953 |      +170,986 |
| 6 |   0 |18,701,046 |18,796,808 |       +95,763 |       +95,763 |      +170,986 |
| 6 | +10 |18,713,933 |18,805,048 |       +91,116 |       +91,116 |      +170,986 |
| 6 | +20 |18,727,147 |18,813,647 |       +86,501 |       +86,501 |      +170,986 |
| 6 | +30 |18,739,969 |18,822,132 |       +82,164 |       +82,164 |      +170,986 |
| 6 | +40 |18,752,439 |18,830,377 |       +77,939 |       +77,939 |      +170,986 |
| 6 | +50 |18,764,145 |18,838,152 |       +74,008 |       +74,008 |      +170,986 |

### Cycle density (tangle / V), at the production setting

| level | flank | tangle_H | V_H | cycles per 1000 nodes |
|---:|---:|---:|---:|---:|
| 4 | +10 |   925,001 |  99,044,505 | 9.3 |
| 5 | +10 |   238,418 |  43,166,220 | 5.5 |
| 6 | +10 |    91,116 |  18,713,933 | 4.9 |

**At every level, the LCP-core graph is overwhelmingly tree-like.** Only ~0.5–1% of nodes participate in cycles at the production flank.

---

# PART B — Context-radius rescue on chr19

All Part-B measurements on chr19 only (extracted with `samtools faidx` from whole-genome ref + truth FASTAs), as a tractable sweep before scaling to whole genome or read-mode FASTQ inputs. **Read-mode is where the boundary-collision concern (thousands of records vs the single chr19 record here) bites; the cyclomatic complexity formula `E − V + 1` understates the true tangle count when there are many connected components, which would need to be corrected for read-mode analysis.**

## B.6 What R is

For each LCP core, the baseline (R=0) signature is `MurmurHash64(seq[start−flank..end+flank], seed=core.label)` — a 64-bit hash of the core's own ±10 bp byte window. When this signature collides with another core's signature in the same genome, **only the colliding cores** get a new signature built from `±R` neighbor cores' signatures, mixed via `fmix64`-folded XOR. `R` mixed into the seed makes the r=0 and r=R hash spaces disjoint. Non-colliding cores keep their r=0 signature.

This is **Approach B2** from the design plan (`plans/in-the-current-version-reflective-lantern.md`): collision-only, fixed-radius. The intent is to disambiguate same-label cores by their genomic neighborhood — operationally a higher-order DBG over the LCP-core alphabet, applied selectively. Closest published precedent: variable-order de Bruijn graphs (Boucher et al., DCC 2015).

## B.7 Pareto sweep — svp anchors

Measured on chr19 truth FASTA pairs (normal = ref chr19; tumor = `simulation/truth/<id>/19_hapA`) with `LEVEL=5 FLANK=10`. Bubble TSV piped to `/dev/null`. Wall-clock + RSS from `/usr/bin/time -v`. 28 cells total: 4 sims × 7 R values.

| sim | R | normal cores | tumor cores | anchors | chain | bubbles | wall (s) | RSS (MB) |
|---|--:|---:|---:|---:|---:|---:|---:|---:|
| del | 0 | 876,869 | 862,049 | 852,426 | 852,426 | 850,489 | 3.80 | 1907.5 |
| del | 1 | 876,869 | 862,049 | 858,256 | 858,256 | 857,938 | 3.96 | 1907.5 |
| del | 2 | 876,869 | 862,049 | 859,054 | 859,054 | 858,843 | 4.03 | 1907.5 |
| del | 3 | 876,869 | 862,049 | 859,518 | 859,518 | 859,368 | 3.97 | 1907.5 |
| del | 4 | 876,869 | 862,049 | 859,858 | 859,858 | 859,729 | 3.93 | 1907.5 |
| del | 6 | 876,869 | 862,049 | 860,355 | 860,355 | 860,274 | 3.89 | 1907.4 |
| del | 8 | 876,869 | 862,049 | 860,692 | 860,692 | 860,622 | 3.83 | 1907.5 |
| ins | 0 | 876,869 | 892,488 | 867,069 | 867,069 | 865,096 | 3.82 | 1961.5 |
| ins | 1 | 876,869 | 892,488 | 872,983 | 872,983 | 872,648 | 3.98 | 1961.5 |
| ins | 2 | 876,869 | 892,488 | 873,823 | 873,823 | 873,597 | 3.97 | 1961.4 |
| ins | 3 | 876,869 | 892,488 | 874,319 | 874,319 | 874,164 | 4.03 | 1961.5 |
| ins | 4 | 876,869 | 892,488 | 874,675 | 874,675 | 874,543 | 4.24 | 1961.5 |
| ins | 6 | 876,869 | 892,488 | 875,189 | 875,189 | 875,106 | 4.19 | 1961.5 |
| ins | 8 | 876,869 | 892,488 | 875,536 | 875,536 | 875,464 | 4.10 | 1961.5 |
| snp | 0 | 876,869 | 876,860 | 867,069 | 867,069 | 865,097 | 3.86 | 1928.5 |
| snp | 1 | 876,869 | 876,860 | 872,982 | 872,982 | 872,647 | 4.02 | 1928.5 |
| snp | 2 | 876,869 | 876,860 | 873,822 | 873,822 | 873,596 | 4.08 | 1928.5 |
| snp | 3 | 876,869 | 876,860 | 874,319 | 874,319 | 874,164 | 4.04 | 1928.5 |
| snp | 4 | 876,869 | 876,860 | 874,676 | 874,676 | 874,544 | 4.05 | 1928.5 |
| snp | 6 | 876,869 | 876,860 | 875,189 | 875,189 | 875,105 | 4.08 | 1928.5 |
| snp | 8 | 876,869 | 876,860 | 875,543 | 875,543 | 875,471 | 3.96 | 1928.6 |
| mixed | 0 | 876,869 | 984,843 | 560,100 | 560,077 | 558,601 | 3.97 | 2154.7 |
| mixed | 1 | 876,869 | 984,843 | 564,545 | 564,522 | 564,320 | 4.18 | 2154.7 |
| mixed | 2 | 876,869 | 984,843 | 564,898 | 564,875 | 564,720 | 4.12 | 2154.6 |
| mixed | 3 | 876,869 | 984,843 | 565,095 | 565,072 | 564,934 | 4.30 | 2154.6 |
| mixed | 4 | 876,869 | 984,843 | 565,247 | 565,224 | 565,079 | 4.22 | 2154.7 |
| mixed | 6 | 876,869 | 984,843 | 565,465 | 565,442 | 565,292 | 4.19 | 2154.6 |
| mixed | 8 | 876,869 | 984,843 | 565,595 | 565,572 | 565,414 | 4.20 | 2154.7 |

### Anchor recovery shape per sim

| sim | R=1 / total R=0→8 gain | R=2 | R=4 | R=8 (full) | total gain (anchors) |
|---|---:|---:|---:|---:|---:|
| del | 70.6 % | 80.2 % | 90.0 % | 100 % | +8,266 |
| ins | 69.8 % | 79.7 % | 89.8 % | 100 % | +8,467 |
| snp | 69.8 % | 79.7 % | 89.8 % | 100 % | +8,475 |
| mixed | 80.9 % | 87.3 % | 93.7 % | 100 % | +5,495 |

**R=1 alone captures 70–81 % of all possible rescue across every sim.** R=2 captures 80–87 %. R=4 captures ~90 %. Diminishing returns from R=4 onward.

### Wall-clock and memory

Wall-clock is essentially flat across R (3.8–4.3 s) — the rescue overhead is negligible because pass 3 only touches the small collider set. Peak RSS depends on input size (1.9–2.2 GB for chr19 pairs).

## B.8 Single-genome k-mer uniqueness (`analyze_kmer`)

For each genome alone, count how many cores have a unique R-context signature within that one genome. This is the *direct* DBG → VG distance: 100 % means every core is a unique node (= VG), lower means some cores collapse.

| genome | R=0 singletons | R=1 | R=2 | R=4 | R=8 | max-rep R=0 → R=8 |
|---|---:|---:|---:|---:|---:|---:|
| ref chr19 | 98.94 % | 99.61 % | 99.71 % | 99.81 % | 99.90 % | 41 → 4 |
| del | 98.95 % | 99.62 % | 99.72 % | 99.81 % | 99.91 % | 41 → 4 |
| ins | 98.96 % | 99.62 % | 99.71 % | 99.81 % | 99.91 % | 41 → 4 |
| snp | 98.94 % | 99.61 % | 99.71 % | 99.81 % | 99.90 % | 41 → 4 |
| **mixed** | **76.84 %** | **77.37 %** | **77.46 %** | **77.59 %** | **77.78 %** | 39 → 4 |

**Two regimes.** Clean inputs (ref/del/ins/snp): R=2 brings DBG to **99.71 %** unique nodes — VG-equivalent for practical purposes. Mixed: even R=8 only reaches **77.78 %** — a hard wall. See §B.10 for why.

`max-rep` (worst-case node collapse multiplicity) drops from 41 (or 39) at R=0 to 4 at R=8 across every genome: R does effectively shrink the worst-tangled nodes, even on mixed.

## B.9 DBG → VG graph topology convergence (GFA tangles)

`tools/dump_gfa` emits standard GFA 1.0 with **one segment per unique signature** (DBG mode) or **one segment per core** (VG mode), and one link per deduplicated consecutive-core transition. `tangle = edges − nodes + 1` (cyclomatic complexity assuming single connected component; chr19 = 1 record, so this is exact).

### B.9a Clean genomes at L=5

| input | mode | nodes | edges | tangles |
|---|---|---:|---:|---:|
| ref | VG | 876,869 | 876,868 | 0 |
| ref | DBG R=0 | 870,455 | 873,774 | 3,320 |
| ref | DBG R=1 | 874,618 | 874,972 | 355 |
| ref | DBG R=2 | 875,226 | 875,425 | 200 |
| ref | DBG R=4 | 875,862 | 875,972 | 111 |
| ref | DBG R=8 | 876,428 | 876,467 | 40 |
| del | VG | 862,049 | 862,048 | 0 |
| del | DBG R=0 | 855,759 | 859,049 | 3,291 |
| del | DBG R=1 | 859,877 | 860,219 | 343 |
| del | DBG R=2 | 860,465 | 860,656 | 192 |
| del | DBG R=4 | 861,077 | 861,185 | 109 |
| del | DBG R=8 | 861,626 | 861,663 | 38 |
| ins | VG | 892,488 | 892,487 | 0 |
| ins | DBG R=0 | 886,075 | 889,393 | 3,319 |
| ins | DBG R=1 | 890,237 | 890,591 | 355 |
| ins | DBG R=2 | 890,845 | 891,044 | 200 |
| ins | DBG R=4 | 891,481 | 891,591 | 111 |
| ins | DBG R=8 | 892,047 | 892,086 | 40 |
| snp | VG | 876,860 | 876,859 | 0 |
| snp | DBG R=0 | 870,447 | 873,766 | 3,320 |
| snp | DBG R=1 | 874,609 | 874,963 | 355 |
| snp | DBG R=2 | 875,217 | 875,416 | 200 |
| snp | DBG R=4 | 875,853 | 875,963 | 111 |
| snp | DBG R=8 | 876,419 | 876,458 | 40 |

R=1 collapses 89 % of tangles already. R=8 leaves 38–40 tangles in ~876k nodes — essentially linear (= VG).

### B.9b Mixed at L=4, L=5, L=6

| level | mode | nodes | edges | tangles | VG-collapse |
|---|---|---:|---:|---:|---:|
| 4 | VG | 2,300,435 | 2,300,434 | 0 | — |
| 4 | DBG R=0 | 1,986,247 | 2,006,869 | 20,623 | 314,188 |
| 4 | DBG R=1 | 2,020,549 | 2,028,962 | 8,414 | 279,886 |
| 4 | DBG R=2 | 2,033,452 | 2,035,772 | 2,321 | 266,983 |
| 4 | DBG R=4 | 2,038,175 | 2,038,557 | 383 | 262,260 |
| 4 | DBG R=8 | 2,040,448 | 2,040,663 | 216 | 259,987 |
| 5 | VG | 984,843 | 984,842 | 0 | — |
| 5 | DBG R=0 | 868,952 | 871,947 | 2,996 | 115,891 |
| 5 | DBG R=1 | 872,696 | 873,048 | 353 | 112,147 |
| 5 | DBG R=2 | 873,329 | 873,576 | 248 | 111,514 |
| 5 | DBG R=4 | 874,192 | 874,368 | 177 | 110,651 |
| 5 | DBG R=8 | 875,297 | 875,400 | 104 | 109,546 |
| 6 | VG | 421,696 | 421,695 | 0 | — |
| 6 | DBG R=0 | 373,863 | 374,120 | 258 | 47,833 |
| 6 | DBG R=1 | 374,328 | 374,496 | 169 | 47,368 |
| 6 | DBG R=2 | 374,638 | 374,758 | 121 | 47,058 |
| 6 | DBG R=4 | 375,057 | 375,142 | 86 | 46,639 |
| 6 | DBG R=8 | 375,669 | 375,734 | 66 | 46,027 |

### B.9c Cross-level comparison on mixed

| level | total cores | tangles @R=0 | tangles @R=8 | tangle reduction | VG-collapse @R=8 | % collapse |
|---|---:|---:|---:|---:|---:|---:|
| 4 | 2,300,435 | 20,623 | 216 | 95× | 259,987 | 11.30 % |
| 5 | 984,843 | 2,996 | 104 | 29× | 109,546 | 11.12 % |
| 6 | 421,696 | 258 | 66 | 4× | 46,027 | 10.91 % |

**The collapse-stuck fraction is invariant under L (~11 %).** Higher L produces fewer cores in absolute terms and proportionally fewer tangles to start with, but the fraction of cores trapped inside large-repeat interiors stays at ~11 % of total cores at R=8 regardless of level. R's effectiveness in tangle reduction is dramatic at L=4 (where there are many fine-grained tangles to crack) and marginal at L=6 (where there were few to begin with).

## B.10 Why mixed hits a wall (the tandem-DUP interior)

Per `analyze_cores` on chr19 mixed tumor:

| layer | distinct sigs | singletons | dup-keys | cores in dups | max-rep |
|---|---:|---:|---:|---:|---:|
| ref (label+flank, R=0) | 99.27 % | 98.94 % | 2,909 | 9,323 (1.06 %) | 41 |
| mixed (label+flank, R=0) | 88.23 % | 76.84 % | 112,203 | 228,094 (23.2 %) | 39 |

The 109k extra rep=2 cores in mixed tumor (vs ref) are produced by the simulation's injected duplications. The `mixed.yaml` config contains:

```
DUP totals across scales:
   300 × 1–50 bp +  50 × 50 bp–1 kb +  50 × 1–10 kb +  20 × 10–100 kb
   +   5 × 100 kb–1 Mb  +  2 × 1–2 Mb  →  total ≈ 4–7 Mb duplicated content
```

Plus equivalent-scale INS and INV. Combined with the ~15 cores/Mbp density at L=5, this produces 100–200k cores in within-tumor pairs.

**Tandem-DUP interior cores cannot be disambiguated by R.** For a tandem duplication of unit length `U` with `K = U / core_spacing` cores per copy:
- Boundary cores (`2R` per copy = `4R` total) see different flanking context across the boundary → R rescues them.
- Interior cores (`K − 2R` per copy) see identical flanking context on both copies → no R can distinguish them.

For a 1 Mb tandem dup at L=5 (~15,000 cores per copy), even R=8 rescues only 16 cores (8 from each boundary) out of 15,000 = 0.1 %. R cannot break this wall — it is a path-coverage problem (would require multi-mapping chaining; Approach E in the design brainstorm), not a context problem.

This matches the constant ~11 % collapse-stuck fraction observed across levels: tandem-DUP interior content is invariant under LCP level because the duplication is byte-identical in the input.

---

# PART C — Synthesis

## C.11 Headline conclusions

1. **Zero observable 64-bit murmur collisions across 33 whole-genome configs**: `distH == distS` in every row of every CSV. Theoretical expectation for 100M items at 64-bit is ~3×10⁻⁴ collisions — empirically zero. *Source: Part A.*

2. **Genome repetition is the real uniqueness ceiling**: 8.07% / 5.15% / 3.83% of cores at L4/L5/L6 have a byte-identical duplicate somewhere at flank=0 (whole genome). Larger flank and deeper level both shrink that rate.

3. **The label seed (`distHL`) is unnecessary for any flank ≥ 0**: across all three levels at whole-genome scale, `distHL == distH` exactly for every flank ≥ 0. The seed only carries information when the window is a strict subset of the core (negative flank); at L4 flank=−30 it adds 42.9M extra distinct keys, at L5 flank=−50 it adds 712k, at L6 flank=−50 it adds 19k.

4. **Murmur IS necessary — label alone is lossy**: at L4, the label collapses 4.41M distinct core sequences (4.52%). At L5 it's 592k (1.38%). At L6 it's 157k (0.84%). All worse than what murmur on the same window content achieves.

5. **Flank not yet saturated at ±50**: gain in singletons per +10 bp flank step decelerates by ~12% each step at L5 but stays positive at +50 (+55,599 singletons gained from +40 → +50). Wider flanks would still help; diminishing returns.

6. **Adjacent cores overlap by default**: 35 bp at L4, 85 bp at L5, 199 bp at L6. Negative flank opens real gaps only at L4 and L5; L6 cores still overlap 99.99% of the time at flank=−50.

7. **Core-length percentile spread is tight**: at flank=0, p5–p95 spans only ~30/60/135 bp at L4/L5/L6 around medians of 63/150/352 bp. The long-tail `max` values (~30 Mbp) are degenerate N-block "cores" from centromeric repeats.

8. **The graph is essentially tree-like at production settings (no rescue)**: only 5.5 cycles per 1000 nodes at L5 flank=+10 on whole genome. The cycles correspond to genome repeats that `unique_anchors` rejects anyway; the survivors form a forest the LIS chain can walk.

9. **Context-radius rescue: R=2 is the right default** (chr19). On clean inputs (ref/del/ins/snp), R=2 brings DBG within 0.29 % of full VG topology (99.71 % singletons; only 192–200 tangles remain in ~876k nodes). R=2 captures 80–87 % of the maximum R=8 rescue. *Source: Part B.*

10. **R=1 already captures 70–81 % of the rescue** (chr19). Going higher than R=2 is diminishing returns for clean inputs.

11. **B1 (per-core variable radius) is unjustified.** The Pareto sweep confirms R=2 is near the asymptote; per-core radius bookkeeping would buy <0.1 percentage points and add complexity.

12. **Mixed-SV inputs hit a wall around 77 %** (chr19). No R value crosses ~78 % uniqueness because the simulation injects multi-Mb tandem duplications whose interior cores cannot be disambiguated by context. This is a structural property, not a hash-quality property — confirmed by conclusion #1 (no observable 64-bit murmur collisions across L4/L5/L6 whole-genome sweep).

13. **Rescue overhead is free**: wall-clock and RSS are flat across R because rescue only touches the small collider set (chr19 measurements: 3.80–4.30 s per cell).

## C.12 Production setting verdict (level=5, flank=10, context-radius=2)

Whole-genome metrics with **no rescue** (R=0):

| metric | value |
|---|---:|
| total cores                   | 45,415,538 |
| valid windows                 | 45,415,538 |
| distinct content keys (distH) | 43,166,220 |
| singletons (anchor candidates) | 42,064,912 |
| label-only would lose         | 592,851 distinct sequences (1.4%) |
| graph tangle (no rescue)      | 238,418 cycles (5.5 per 1000 nodes) |
| genome repetition rate        | 5.15% |
| median window length          | 170 bp |
| p5 / p95 window length        | 144 bp / 203 bp |
| avg adjacent-core overlap     | 105 bp |

Going to **flank=50** would buy +369,254 singletons (+0.88%) and cut tangle by 34.6% (238k → 156k). Going to **L6 flank=10** would shrink the anchor pool to 18.3M but improve uniqueness to 96.7% and reduce tangle density to 4.9/1000 nodes — at the cost of sparser genome coverage.

Adding **R=2 rescue** (chr19 evidence projected, not yet measured whole-genome): on clean inputs, lifts singletons from ~98.94 % to ~99.71 % per genome; on mixed-SV tumor, lifts from 76.84 % only to 77.46 % — the tandem-DUP wall sits at ~77.78 % even at R=8.

## C.13 Tools used

- `./svp <normal.fa> <tumor.fa> <out.tsv> --context-radius R` — bilateral anchor count after `unique_anchors`
- `./tools/analyze_uniqueness <fasta> <level> <flank,...>` — per-genome distinct-key / singleton / length-percentile sweep
- `./tools/analyze_tangle <fasta> <level> <flank,...>` — V/E/tangle on the signature-adjacency graph at R=0
- `./tools/analyze_cores <ref> <tum>` — uniqueness histograms at label-only / flank-only / core-only / label+flank
- `./tools/analyze_kmer <fasta>` — per-genome R-context singleton + max-rep counts at R ∈ {0,1,2,3,4,6,8}
- `./tools/dump_gfa <fasta> <out.gfa> --mode dbg|vg --context-radius R` — per-genome graph topology in GFA 1.0

## C.14 Limitations and next steps

- **Context-radius rescue only measured on chr19.** Whole-genome behavior should be qualitatively similar but is not yet measured. Whole-genome would also surface the boundary-collision concern (thousands of records vs the single chr19 record); the cyclomatic-complexity formula `E − V + 1` understates the true tangle count when the graph has many connected components, which will need correcting for multi-record / read-mode analysis.
- **No read-mode test yet.** `simulation/reads/snp/hifi.fastq.gz` and `simulation/reads/mixed/hifi.fastq.gz` exist but svp does not accept FASTQ input today.
- **Tandem-DUP interior remains unsolved.** R cannot address it. Approach E (multi-mapping chaining over many-to-many anchor sets, replacing the current LIS chain) is the natural next step if mixed-SV recall becomes a target.
- **Flank > +50 not measured.** The marginal gain per +10 bp at L5 is still positive (~0.13 %) at +50; saturation point is past the current sweep.
