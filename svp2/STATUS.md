# svp2 — status as of 2026-04-22

Experimental modular SV-calling pipeline on branch `experimental/modular-pipeline`. Lives alongside the original `svp/` in the same repo; original is untouched.

## Purpose

Tumor-vs-normal structural-variant caller. Reads normal + tumor FASTA, emits a VCF of SVs. Completely different architecture from the original svp: the original is a monolithic C++ binary that does graph+bubble+SV calling in one process; svp2 splits the work into six independent stages that each read the previous stage's binary output, so you can swap one stage's algorithm without retouching the others.

## 6-stage flow

```
FASTA (ref + normal + tumor)
      │
      ▼  s1_lcp.cpp          LCP cores @ depth k (default 5)
  outputs/stage1_<sim>_d<k>/{ref,normal,tumor}.{lcp,seq}.bin + meta.json
      │
      ▼  s2_graph.cpp        2-colour graph + ref anchor map
  outputs/stage2_<sim>_d<k>/graph.bin, ref_anchors.bin
      │
      ▼  s3_bubbles.cpp      --extractor <name>
  outputs/stage3_<sim>_d<k>/<extractor>.bubbles.bin
      │
      ▼  s4_simple.py        two-pointer DEL/INS/(DUP-hint)/SUB
  outputs/stage4_.../<extractor>.{calls,unresolved}.bin
      │
      ▼  s5_align.py         --aligner <name> --threads N
  outputs/stage5_.../<extractor>.<aligner>.calls.bin
      │
      ▼  s6_finalize.py      INS→DUP via edlib HW + truvari-style dedupe
  outputs/vcf/<sim>_d<k>_<ext>_sv<min>_<aln>_dup<t>_dd<p>.vcf
```

Binary formats are defined once in `svp2/src/svp_io.{h,py}` (magic + version + records). Both C++ and Python read/write them.

## Extractors (stage 3)

| name | idea | status |
|---|---|---|
| `chain_rt` | ref∩tumor unique-LCP-core intersection, emit one bubble per consecutive anchor pair | **best**, O(n), finds all truth SVs as bubbles |
| `chain_3w` | ref∩normal∩tumor three-way unique intersection | essentially ties chain_rt |
| `bfs` | classic color-divergent BFS on graph (≈ original svp) | works, lower F1 |
| `bulge` | SPAdes-style 2-path bulge with hop/length budget | works, worst F1 |
| `superbubble` | Brankovic linear-time superbubble detection | **broken on cyclic LCP graph** (topo sort covers <1% of nodes); needs SCC condensation |

## Aligners (stage 5)

| name | backend | notes |
|---|---|---|
| `wfa` | pywfa (WFA2 bidirectional) | primary, bails to edlib above 30kbp |
| `wfa_rc` | pywfa + inner-strip + k-mer Jaccard on revcomp | adds INV detection; speculative INV for same-length bubbles (caveat below) |
| `edlib` | edlib NW path | fast fallback |
| `parasail` | parasail NW SIMD | optimal DP, slow above ~2kbp |
| `mappy` | minimap2 asm5 | auto-detects INV via strand flip |

## Benchmarks (truvari v5.4.0: `sizemin=30 pctseq=0.70 pctsize=0.70 refdist=500`)

All numbers are against **single-SV-type insilicoSV simulations** on chr19. Truth files live in `simulation/`.

| sim file | truth count (post-filter) | F1 | prec | rec | TP | FP | FN | best config |
|---|---|---|---|---|---|---|---|---|
| `deletion_simulation.vcf` | 155 | **0.997** | 1.000 | 0.994 | 154 | 0 | 1 | chain_rt + wfa |
| `insertion_simulation.vcf` | 156 | **0.833** | 0.833 | 0.833 | 130 | 26 | 26 | chain_rt + wfa_rc (regressed — see Known issues) |
| `duplication_simulation.vcf` | 154 | **0.842** | 0.853 | 0.831 | 128 | 22 | 26 | chain_rt + wfa |
| `inversion_simulation.vcf` | 165 | **0.890** | 0.952 | 0.836 | 138 | 7 | 27 | chain_rt + wfa_rc |
| `mixed_v2.vcf` | ~975 | **not benchmarked** | — | — | — | — | — | grid stalled before run completed |

Pre-regression run (before INV speculative-emit fix) had **INS F1=1.000** perfect.

Mixed, on a single pre-fix smoke test: `f1=0.784 p=0.755 r=0.814 tp=794 fp=257 fn=181` (chain_rt + wfa_rc) — treat as indicative only.

Mean F1 by **extractor** across the 52 completed grid configs:
```
chain_3w   0.872
chain_rt   0.802
bfs        0.697
bulge      0.487
```

Mean F1 by **aligner**:
```
wfa_rc   0.775
wfa      0.756
mappy    0.690
edlib    0.668
```

## Known issues / regressions

1. **INS regression 1.000 → 0.833**. In `s5_align.py::_cigar_tuples_wfa_rc`, same-length bubbles with inner ≥ 20bp get a speculative `__INV_INNER__` emit at confidence 0.15 even without revcomp-kmer evidence. This rescued INV from F1=0.03 to F1=0.89 but introduced 26 bogus INV calls on the pure-INS sim. Fix: gate the speculative path on a minimal k-mer-Jaccard rev-comp signal, or move speculative emit out of `wfa_rc` and into a separate `inv_speculative` aligner.

2. **DUP at 0.842 plateau**. chain_rt produces a bubble for 162/162 truth DUPs ≥30bp, but only ~128 get matched by truvari. Specific failure mode: truth DUP sv101 at 31088593 is "close" to our bubble at 31088681 (only 88bp apart and SVLEN=221 matches exactly) but the sequences of our INS alt and the truth DUP source have 58% divergence — the LCP anchors don't land where the tandem-duplication boundary does, so the bubble's inserted sequence isn't literally the DUP source. I ran out of time before sweeping the DUP threshold (0.55/0.65 candidates in the drafted grid).

3. **Mixed sim stall**. `bfs` extractor on `mixed_v2` was producing 93k spurious SNP/SUB calls that stage 6's O(N²) INS→DUP edlib loop choked on. Fixed three bugs (stage 5 SNP emission, stage 6 dedupe O(N²), 30kb cap on `ins_is_dup`), but the restarted grid was interrupted before it reached mixed_v2.

4. **superbubble extractor emits 0 bubbles**. The Brankovic topo sort only covers ~840 of 868k nodes because the colored graph contains many SCCs (expected for LCP with repeats). Needs SCC condensation before topo sort.

5. **INS→DUP conversion bails above 30kb**. Anything larger just stays as INS and becomes a FN when truth uses `<DUP>`. 3 calls in the INS sim, 1 in DUP.

## Grid state

- `svp2/outputs/metrics/grid.jsonl`: 52 completed configs out of 134 planned (v1 grid in `scripts/grid_search.py`).
- `svp2/scripts/grid_v2.py` drafted (~120 additional configs: DUP threshold sweeps, INV depth sweeps, refdist sweeps) — never executed.
- 109 intermediate VCFs cached in `svp2/outputs/vcf/`.

## Build & run

```bash
# Build the C++ stage binaries
make -C svp2 all          # → svp2/bin/{s1_lcp,s2_graph,s3_bubbles}

# One config end-to-end (bash wrapper, caches every stage)
svp2/scripts/run_pipeline.sh <sim> <extractor> <aligner>
#   sim       ∈ {del, ins, dup, inv, snp, mixed_v2}
#   extractor ∈ {chain_rt, chain_3w, bfs, bulge, superbubble}
#   aligner   ∈ {wfa, wfa_rc, edlib, parasail, mappy, none}

# Truvari benchmark (sizemin=30, returns summary JSON on stdout)
svp2/scripts/bench.sh <sim> <path/to/output.vcf>

# Leaderboard / aggregate
svp2/scripts/analyze_grid.py svp2/outputs/metrics/grid.jsonl

# Grid run (restart from where it left off — resumes via jsonl keys)
python3 svp2/scripts/grid_search.py --parallel 6 --workers-stage5 4

# Drafted v2 sweeps
python3 svp2/scripts/grid_v2.py --parallel 6 --workers-stage5 4
```

## Data

Pre-extracted chr19-only FASTAs live in `svp2/data/` (58–60 MB each, generated from `simulation/*.fa`). `run_pipeline.sh` uses these by name — no need to re-slice the 3 GB whole-genome FASTAs on every run.

## What to do next (in priority order)

1. Gate the speculative-INV emit on a minimal revcomp-kmer signal to recover INS F1=1.000 without giving up INV F1=0.89. `s5_align.py::_cigar_tuples_wfa_rc`.
2. Resume the v1 grid — 82 configs left, ~25 min at 6× parallelism. It already includes the DUP threshold sweep (0.55/0.65/0.80/0.90) that should raise DUP F1.
3. Run the drafted `grid_v2.py` (121 configs) for deeper dup_thr × pctseq × refdist × depth sweeps.
4. Benchmark mixed_v2 end-to-end — the stage 6 stall is fixed, so this should just work now.
5. Fix the `superbubble` extractor via SCC condensation, or delete it and replace with a sixth extractor worth having.

## File map

```
svp2/
  Makefile                 builds the 3 C++ stage binaries
  README.md                architecture + binary format spec
  STATUS.md                this file
  src/
    svp_io.{h,py}          binary format (magics SVSQ/SVGR/SVRA/SVBB/SVCL)
    svp_hash.h             FNV-1a 64-bit + reverse_complement
    svp_graph.{h,cpp}      2-colour graph + read/write + Bubble struct
    svp_bubble.{h,cpp}     the 5 extractors
    s1_lcp.cpp             stage 1 driver
    s2_graph.cpp           stage 2 driver
    s3_bubbles.cpp         stage 3 driver (--extractor)
    s4_simple.py           stage 4 — two-pointer
    s5_align.py            stage 5 — parallel aligners
    s6_finalize.py         stage 6 — merge + DUP + VCF
  scripts/
    run_pipeline.sh        single-config end-to-end (cached)
    bench.sh               truvari wrapper
    grid_search.py         v1 grid (134 configs)
    grid_v2.py             v2 grid (drafted, ~120 configs)
    analyze_grid.py        leaderboard from grid.jsonl
  data/                    chr19-only FASTAs (gitignored)
  outputs/                 all stage outputs (gitignored)
  bin/                     compiled binaries (gitignored)
```
