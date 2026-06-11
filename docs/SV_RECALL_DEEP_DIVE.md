# SV calling deep dive + the strand=1 fix (2026-06-11)

Overnight investigation of three problems on `wgs_simit` (Sim-it Q30 HiFi; truth ≥50bp:
96 DEL / 52 DUP / 26 INV / 26 INS), and the one clean fix it produced. Full investigation
scaffolding (env-gated modeb/vcf instrumentation, analysis tools) lives on branches
`p1-detour-bubbles` / `p2-ins-del-accuracy` / `p3-dup-inv`. This file is the standalone
account on the working line.

## The fix that shipped: drop strand=1 bubbles

The caller does dual fw+RC LCP extraction per read, so every locus is represented twice:
a **strand=0** bubble and its **strand=1** RC mirror. The strand=0 arm coordinate maps
linearly to the forward genome; the strand=1 arm does **not** (its middle super-nodes can sit
at unrelated forward positions), so no POS formula lands it. Measured on the production
`bubble_pos.bin` + `modeb.bin` against the truth SNPs:

| | SNP bubbles | exact-POS hit (best of 6 formulas) |
|---|---:|---:|
| strand=0 | 15,786 | **98.1%** |
| strand=1 | 16,761 | **1.5%** |

strand=0 alone already accounts for ~all SNP recall (15,483 of 15,484 TP), so the strand=1
mirrors add ≈0 recall while producing ~all the FP. Dropping them (`vcf.cpp::apply_filters`,
one line `if (bp.strand == 1) return false;`):

| tier | TP | FP | recall | precision | **F1** |
|---|---:|---:|---:|---:|---:|
| SNP        | 15,484 (=) | 17,116 → **313** | 0.913 (=) | 0.48 → **0.98** | 0.625 → **0.945** |
| indel 1–49 | 850 → 838 | 1,131 → **32**  | 0.904 → 0.891 | 0.43 → **0.96** | 0.582 → **0.926** |
| ≥50bp SV   | 99 (=)     | 9 → **1**       | 0.495 (=) | 0.92 → **0.99** | 0.643 → **0.660** |

Recall is unchanged on SNP and SV and loses only 12 indels (strand=1-only). The ≥50bp
duplicate FPs (the truvari `Multi` calls) were strand=1 mirrors — they vanish too. This is
the cheap realisation of the "strand=1 POS" upstream lever: rather than *reposition* the
mirror (needs a per-super-node posmap walk), we *drop* it, since strand=0 already carries the
variant correctly. The residual loss (12 indels, ~258 strand=1-only SNPs) is recoverable
later by the proper per-super-node POS fix.

## Problem 1 — where `T_only_no_N` / `N_only_no_T` detour bubbles come from

Two populations, mixed:

**(A) Split halves of real large homozygous SVs.** Verified node-by-node at DEL `19:16995558`
(svlen 45903): the deletion is bracketed by a bucket=1 I-tag **junction** between the two
true flank anchors (covn=0, covt=28, span=45974≈svlen) **plus** bucket=2 D-tag **segments**
(covn=28, covt=0, n_resid≈svlen). They never merge into one bucket=0 bubble because an
**embedded c3 super-node inside the deleted segment** (a sub-sequence repeated elsewhere →
shared) breaks the c1 reference path so it closes at a different anchor than the junction.
This is why small SVs stay bucket=0 and large ones fragment.

**(B) Repeat-region / INV-locus artifacts.** Repeat-region LCP-determinism makes normal and
tumor pick different cores → spurious single-colour bubbles; INV loci emit a cluster (median
6 orphan bubbles within ±3kb of a truth INV vs 2 control).

**Classification is not solvable by local features.** Five discriminators measured
(local features, span-as-svlen, half-pairing, c1|c3 reconstruction, coverage-split) — all
precision 0.03–0.06. Real split-halves and repeat/INV phantoms are locally identical.

## Problem 2 — why INS/DEL recall caps at ~80/85%

`tools/trace_truth.py`: of 19 FN DELs, 1 has no bubble (upstream loss) and 18 are the
fragmented orphans above (correct size, dropped on `bucket!=0` + `pre+suf=0`). Letting them
through (modeb c1|c3 reconstruction + relaxed gate) lifts DEL 77→84 but FP 6→137 (F1
0.64→0.37); **155 of the 269 new FP sit on truth INVs, 39 on DUPs** — i.e. mis-typed real
structure, not noise. The `bucket=0` gate is effectively an INV/DUP-locus suppressor, which
is why orphan recovery is net-negative. So the remaining ≥50bp recall gap is structural
(fragmentation + untypeable INV/DUP), not a tunable. The FP side, by contrast, *was* fixable —
it was the strand=1 duplicates (above).

## Problem 3 — DUP and INV

**DUP:** coverage signature is clean (c3 bubbles inside truth DUPs have median covt/covn=2.00,
93% ≥1.5; the DUP junction bid 32474 gives exact boundaries) but coverage *runs* give 0
size-matched DUPs (boundaries too coarse for truvari `pctsize 0.7`) and only ~1 DUP forms a
clean junction bubble. **INV:** structurally invisible — 0 tumor-only (c2) sigs at INV loci
(dual fw+RC extraction erases orientation), so no byte-level detector can work; detectable
only as a locus cluster, with no recoverable svlen. Both remain upstream-blocked.

## Remaining upstream levers (ranked)

1. **per-super-node strand=1 POS** — recovers the 12 indels + ~258 SNPs the drop forgoes, and
   would let strand=1 SV junctions be used. (The drop is the 95%-of-benefit shortcut.)
2. **fragmentation by embedded c3 nodes** — needs an SCC-/posmap-aware walk in modeb so large
   SVs form one bucket=0 bubble (the naive c1|c3 mask wanders the giant c3 SCC: D 2180→1710,
   C 464→877).
3. **INV orientation** — needs single-strand sig provenance at read stage.
4. **RC-canonicalise sigs** — folds strand pairs, halving repeat phantoms (and is the
   structural version of the strand=1 drop).
