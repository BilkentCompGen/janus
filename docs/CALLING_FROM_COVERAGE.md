# Calling DEL and DUP from coverage alone — what actually works, measured

_2026-07-29, on the reference-free graphs (`hard`, `complex`: both strands, no blacklist, no vote).
Every number below comes from a scan of the WHOLE graph with no knowledge of the truth; the truth is
consulted only afterwards, to say what each candidate was. Earlier truth-conditioned measurements are
superseded — they measured recall and called it precision._

Tooling: `tools/gt_cores` (truth → core sigs), `scripts/gt_vs_graph.py` (join + threshold sweep),
`scripts/scan_dosage.py` (genome-wide scan), `scripts/scan_strat.py` (size-stratified attribution).

## The detector

```
dosage(node) = tcov / ncov                          per node, no reference, no coordinates
DEL candidate = a connected run of >= X nodes with dosage < 0.1 and ncov > 0
DUP candidate = a connected run of >= X nodes with dosage > 1.5, ncov inside the backbone band
backbone      = median ncov over nodes present in both samples (measured: 25 at 30x)
```

Connectivity is the graph's own edges. A candidate is scored against the truth by SIGNATURE.

## DEL — works, and the threshold falls out of the biology

Stratified by run length (`hard`, 87,986 candidates at X=3):

| run (nodes) | candidates | DEL ≥500bp | SMV (<50bp) | matches nothing |
|---|---|---|---|---|
| 3–4 | 83,741 | 1.0% | **93.8%** | 3.7% |
| 5–6 | 1,678 | 29.5% | 66.9% | 0.6% |
| 7–9 | 564 | **87.9%** | 5.1% | 0.4% |
| 10–14 | 584 | **99.3%** | 0.2% | 0.0% |
| 15+ | 1,419 | **100.0%** | 0.0% | 0.0% |

| X | calls | precision (as DEL) | precision (as any SV) | recall |
|---|---|---|---|---|
| 5 | 4,245 | 70.4% | 72.6% | 92.7% |
| **6** | 2,990 | **91.0%** | 92.9% | **92.3%** |
| **8** | 2,352 | **98.7%** | 99.7% | **91.4%** |
| 12 | 1,724 | 100.0% | 100.0% | 87.1% |

On `complex` the same rule gives **69.1% as DEL / 100.0% as any SV** at X=8 (recall 91.8%): the residual
is complex-grammar deletion blocks — real somatic loss, correctly located, whose event-level label is CPX.

**Why X≈6.** A SNV destroys the 2–3 reference cores it falls in and can never make a longer run. The
core step is ~77 bp, so 6 cores ≈ 460 bp — the threshold lands on the 500 bp boundary by itself, not by
tuning.

**Recall falls above X=12 for one boring reason:** size. A 500–1000 bp deletion has a median of 9 cores
per strand family, so raising X excludes it by construction. Nothing is being lost by the detector.

**Fragmentation is not a problem below 20 kb.** Per strand family, deletions of 500 bp–20 kb come out as
**one unbroken run** (median 1 run, longest/available = 1.00). Only >20 kb events fragment (median 4
runs, longest covering 59%). The "~14 fragments per event" figure in an earlier note was wrong: it
counted both strand families as fragments of one run.

## DUP — the ratio test alone is nearly useless; the fix is two-sided

`dosage > 1.5` alone selects **8.47% of the graph** (7.45 M nodes) and gives 14.6% precision at X=20.
The mask is not full of repeats — it is full of noise. Its nodes have **median ncov 18 against a
backbone of 25**, with median tcov 31: the ratio fires because ncov drifted DOWN, not because tcov went
up. At mean 25 and sd ~5 that pair is an ordinary Poisson fluctuation, which is why 8.5% of all nodes
qualify.

A real tandem duplication has ncov at ONE copy and tcov at TWO. Constraining ncov to the backbone band:

| rule | precision @X=8 | precision @X=20 | recall @X=20 |
|---|---|---|---|
| `dosage > 1.5` | 7.2% | 14.6% | 72.4% |
| + `ncov >= backbone/2` | 8.0% | 19.7% | 72.4% |
| + **`ncov` in [0.75, 1.5]×backbone** | **42.5%** | **83.3%** | 66.7% |

## The coverage-derived repeat mask: right diagnosis, wrong cure

A core occurring at C places collects ~C× the depth **in both samples** — copy number, not tumour excess —
so repeats can be flagged without any reference. Measured inside deletions ≥20 kb:

| | median ncov | 90th | fraction with ncov>100 |
|---|---|---|---|
| cores that BREAK the run (3.4% of the interior) | **96** (≈4 copies) | 18,228 | **49.2%** |
| cores that pass | 24 | 32 | 1.0% |

The identification is correct. **Using it to bridge runs is not.** Letting a run pass through repeat
nodes drops DEL precision from 98.7% to 91.7% (REP=1.75) and only recovers 1.7% recall at best (REP=10,
94.3%). A repeat core is a genome-wide hub: crossing it merges unrelated loci into one run and drags SMV
noise into the long bins. Restricting the bridge to low-degree repeats (degree ≤ 4) changes **nothing at
all** — the breakers are exclusively high-degree hubs.

So the repeat flag's real uses are elsewhere: excluding repeats from being anchors, and (untested)
sizing. It is not the fix for large-deletion fragmentation, and fragmentation was not costing anything.

## What is still open

- **DUP recall** stalls at ~67–77%; the two-sided depth band should be tightened around `tcov ≈ 2×
  backbone` rather than only bounding ncov.
- **Merging** candidates into events. At X=8 the 2,352 DEL calls cover 213 distinct deletions — the
  strand-mirror pair plus, for >20 kb, real fragments. A count of events cannot be claimed before this.
- **INV has no dosage signal at all** (measured: 1.04) and never will: with both strands inserted, an
  inversion's interior covers the mirror rail exactly as deeply as the forward one. An inversion is its
  two junctions, and separating an INV junction from a DEL/TRA junction needs the mirror-twin relation,
  which nothing currently records.
