# Bubble-caller comparison: position-coverage recall on the chr19 LCP graph

Goal: maximize **bubble-coverage recall** — the fraction of ground-truth SVs whose position is covered by ≥1 bubble. This is the upstream ceiling for any SV caller built on the graph (a bubble must exist before anything can be typed/called). Type is ignored here; only position matters.

All callers run on the **same** combined colored LCP graph (chr19: ref record `19` + tumor `19_hapA`, level 5, ±10 bp flank, 64-bit signature; node colour 0x1=normal, 0x2=tumor, shared=both). Confirmed identical substrate: `vg-rewrite`/`gfa_collapse` = 1,067,263 nodes; `main` branch = 1,064,375 nodes, shared 670,558 / normal-only 197,570 / tumor-only 196,247 (<0.5% apart). GT = 1,702 SVs (427 DEL / 427 DUP / 421 INS / 427 INV) from the insilicoSV VCF. Coverage tolerance `W` swept; metric in `tools/bubble_recall.py`.

## Correction: the "558,601 svp bubbles" was a bug, not a flood

An earlier analysis reported svp (`vg-rewrite`) emitting 558,601 bubbles, 99.7% CPLX, and concluded "the caller floods 328×." **Retracted.** Confirmed from a clean rebuild: **556,884 / 558,601 rows (99.7%) have `n_start > n_end`** — backwards coordinates with negative declared lengths, producing a 4.6 GB garbage file. Root cause in `src/match.cpp:emit_bubbles`: consecutive kept anchors' cores **overlap** (~85 bp at L5), so `anchor[i].end > anchor[i+1].start`; `ext_r`/`ext_l` no-op and `nl = ns - ne` **unsigned-underflows to ~2⁶⁴**. svp's **real** output is **1,717 well-formed bubbles** (570 CPLX / 203 DEL / 558 INS / 386 SNP). Two-line fix: clamp `if (ne_end < ne_start) ne_end = ne_start;` (and the tumor side) before `ext_r`. (`main` branch already guards this correctly via its Phase-5 underflow guard.)

## Bubble-coverage recall (W=1000 bp), same metric, same GT

| caller (same LCP graph) | bubbles | recall | DEL | INS | DUP | INV | precision |
|---|--:|--:|--:|--:|--:|--:|--:|
| **Cortex-style 2-colour (permissive)** | 3,683 | **99.1%** | 100% | 100% | 96% | 100% | 68.1% |
| svp / vg-rewrite greedy (well-formed) | 1,717 | 97.4% | 100% | 100% | 90% | 100% | 78.8% |
| `main` — ALL VCF calls (bubbles + core-index pipelines) | 2,061 | 96.3% | 99% | 99% | 90% | 98% | 82.5% |
| Onodera strict superbubble | 1,736 | 84.4% | 76% | 69% | 94% | 98% | 86.8% |
| **`main` — BUBBLES only** | **1,650** | **79.7%** | **75%** | **66%** | 83% | 94% | 76.7% |

Recall vs W (coverage tolerance): main-bubbles 67.3→73.5→78.0→79.7→81.8% at W=0/100/500/1000/2000; Cortex 98.1→98.1→99.0→99.1→99.5%.

## Why `main`'s bubble coverage (79.7%) is ~19 pts below Cortex (99.1%)

It is **not the graph** (identical substrate). It is the bubble-detection/extraction algorithm. Measured size of covered vs **missed** SVs for `main` bubbles:

| type | covered n | covered med SVLEN | missed n | missed med SVLEN |
|---|--:|--:|--:|--:|
| DEL | 319 | 40 bp | 108 | **28 bp** |
| INS | 279 | 42 bp | 142 | **31 bp** |
| DUP | 356 | 30 bp | 71 | **3,526 bp** |
| INV | 403 | 39 bp | 24 | 40 bp |

So `main` bubbles miss **small DEL/INS** (median 28–31 bp, all below the ~85 bp core overlap) and **large DUP** (median 3.5 kb). Three causes, from `main`'s `ALGORITHM.md` / `sv.cpp`:

1. **Strict convergence requirement (Phase 5).** `findColoredPath` BFS-walks the normal and tumor branches and keeps a bubble only if both converge to the *same* shared end node; "invalid (open-ended) paths are discarded." The permissive Cortex decomposition instead marks **any** maximal run of single-colour nodes between shared anchors — no convergence test.
2. **Overlap-empty guard swallows small variants (Phase 5).** `main` sets the extracted sequence **empty** when "the start and end node overlaps fully cover the intermediate path." Since cores overlap ~85 bp, a DEL/INS smaller than the overlap (the missed median 28–31 bp) becomes empty → `detectSV` → `no_call` → no bubble. (This guard is the *correct* fix for the underflow bug above — but its side effect is dropping sub-overlap variants.)
3. **`skipped_in_inv = 166`** — bubbles inside large-INV regions are deliberately skipped (handed to the INV core-index pipeline).

`main` compensates with **separate non-bubble core-index pipelines** (`detectLargeDELs`=151, `detectLargeINS`=171, `detectLargeINVs`=98, `detectFlippedDUPs`=121), which lift total VCF coverage back to 96.3%. The existence of those four detectors is itself evidence the bubble layer under-covers.

## Tradeoff and the recall-maxing recipe

`main`'s strictness **buys precision** (bubbles 76.7% vs Cortex 68.1%) at the cost of recall. The permissive Cortex decomposition is the opposite end: max recall, lower precision (the "extra" bubbles are real small variants, filtered downstream by typing). To push bubble-coverage recall to ~99%:

- Detect a bubble as **any maximal run of single-colour nodes flanked by shared nodes** — drop the "both branches converge to the same shared end" requirement; keep open-ended/divergent regions.
- For sub-core-overlap variants, **emit the small INS/DEL from the flank breakpoint coordinates** instead of emptying to `no_call` (the breakpoint is well-defined even when the divergent core is overlap-covered).
- Keep large DUP on a **multiplicity / flipped-DUP** path — no graph topology recovers a pure tandem dup as a bubble (the intrinsic ~1% even Cortex leaves: 16 DUPs at 96.3% DUP recall).

This is the design the clean re-implementation targets: LCP cores → colored graph → **permissive** Cortex bubble decomposition, verified at 99.1% / 3,683 bubbles before any typing is added.

## Provenance

- Benchmark: `tools/bubble_recall.py` (position-coverage recall/precision, W sweep, miss analysis). GT: `outputs/gfa/chr19_gt_sv.tsv`.
- Cortex / Onodera intervals computed in `tools/gfa_viz.py` (+ recomputed in the benchmark) on `outputs/gfa/combined_chr19.gfa`.
- `main`-branch bubbles: built in isolated worktree `/tmp/svp-main` (branch `main`, htslib 1.23.1 + parasail), run `./svp --ref ref_chr19.fa --normal ref_chr19.fa --tumor tumor_chr19.fa --depth 5 --chr 19 --debug` → `main_chr19_bubbles.tsv` (1,650), `main_chr19.vcf` (2,061 calls). 26.0 s / 3.2 GB.
- svp/vg-rewrite well-formed bubbles: `outputs/gfa/svp_chr19_wellformed.tsv` (after excluding the underflow rows).
