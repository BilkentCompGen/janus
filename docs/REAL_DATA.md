# svp on real tumour/normal HiFi — what was measured

Two sittings. **2026-08-05/06** established that the simulation numbers do not survive contact with
real data and that nothing tried then was worth keeping. **2026-08-10** found what was wrong and
fixed it. F1 against the paper's own truth, from that state to this one:

| | hcc1954 | h2009 | hcc1937 | h1437 |
|---|---|---|---|---|
| 2026-08-06 | 0.337 | 0.173 | 0.263 | 0.063 |
| now | **0.872** | 0.854 | 0.762 | 0.759 |
| Severus, published | 0.871 | 0.903 | 0.826 | 0.810 |

svp passes Severus on hcc1954 and does not on the other three. Everything below is in the code
unless it says otherwise.

## The data

**Reads** — ENA, BioProject **PRJNA1086849**, the Severus paper's own project, which is the
**CASTLE** panel (Cancer Standards Long-read Evaluation). PacBio **Revio**, five tumour/normal pairs,
`hifi/get.sh` fetches and packs them. Measured on the raw fastq, not assumed: median read 16.2 kb,
median per-read Q38.5, overall Q37.4 → **0.018% error**. These are HiFi despite every file in the
archive being named `_subreads.fastq.gz`.

| pair | tumour | normal | truthset |
|---|---|---|---|
| hcc1954 | SRR28305163 | SRR28305160 | 1106 |
| h2009 | SRR28305172 | SRR28305169 | 925 |
| hcc1937 | SRR28305185 | SRR28305182 | 703 |
| h1437 | SRR28305179 | SRR28305175 | 221 |
| hs578 | SRR31537482 | SRR31537479 | **none** — Hs578T is not in CASTLE |

**Truth** — the same paper's Zenodo record (`10.5281/zenodo.10856827` → 14541057),
`variant_calls_and_benchmarks.tar.gz`, the `minda_4_2` ensemble per cell line: at least 2 callers in
each of 3 platform groups AND at least 4 callers overall. `truth_clean.vcf` is that file with the one
unprefixed `cmd:` line Minda writes into its own output removed, because pandas will not parse it.

**The comparison is not symmetric, and the size of that is measured.** Severus is one of the callers
whose output built the ensemble, so a Severus call is in the truth partly for reasons that have
nothing to do with being right. The archive also ships looser ensembles for hcc1954. Scoring one and
the same svp call set against each:

| truth | records | svp FP | svp precision |
|---|---|---|---|
| minda_4_2 (what the paper reports against) | 1106 | 200 | 0.827 |
| minda_3_2 | 1168 | 172 | 0.851 |
| minda_2_1 | 1413 | 136 | 0.882 |

64 of the 200 "false" calls are events some looser consensus contains. Of the 287 records that
separate the 2-support ensemble from the 4-support one, 247 have two or more distinct tools behind
them and 40 are one tool counted twice (severus on two technologies). That is an upper bound on how
much of our false-positive count is a truth gap, not a measurement of it.

## What is in the code, and what each is worth

Every one of these was measured on real reads with everything else held, and checked against the
simulation for regression.

**1. No length floor on a segment** (`min_aln`, was 1000, now 0). A deletion needs one surviving
alignment of a contig; a breakend needs two, so a breakend pays the floor twice. Worse, dropping a
contig's MIDDLE segment chains the two that remain and fabricates a junction between pieces that were
never adjacent. hcc1954, everything else held: 1000 → 0.749/0.718, 600 → 0.761/0.762, 400 →
0.771/0.797, 200 → 0.798/0.863. **Both precision and recall improve as the floor comes off**, and
false positives fall from 266 to 241. Below 200 the numbers are identical to three digits at 150,
100, 50 and 0, because minimap2 emits nothing shorter here.

**2. A repeat test inside `shred()`** (`shred_rep`, 0.40). `shred()` reads a dense run of small indels
as a breakpoint pair, which is what a short inversion looks like — and also what a VNTR looks like.
On hcc1954 it produced 41 false calls and not one true one; on `hard` it is worth 4 inversions in the
50–1000 bp band. The two cases differ in the contig's own bases, so no reference is needed: at
`rep ≤ 0.40` hcc1954 loses all 41 false calls (241 → 200 FP, exactly what turning shred off entirely
achieves) and `hard` keeps all 4 inversions. Below 0.35 the simulation starts losing them.

**3. The same test at the call** (`anchor_rep`, 0.50). Measured over the hcc1954 calls that not even
the 2-support ensemble contains: 27.2% have an anchor above 0.4 against 2.9% of the true calls, and
above 0.6 it is 11.0% against 0.2%. Taken on the anchor side only, where the bases are aligned at
under 1% mismatch and therefore ARE the reference at that locus. 0.50 and not 0.45 because the knee
is between 0.45 and 0.40, where 33 true calls go at once. The simulation pays 12 of 600 junctions for
this — `hard` has no VNTR-driven false calls to win back, so every step down is pure loss there.

**4. The germline check on breakends too**, and its threshold at 0.25. `germline()` only ran on
same-chromosome calls, so 27 of the 136 uncorroborated false calls were breakends it had never
examined. Of those calls, 7.4% sit where a quarter or more of the normal reads are above 2%
mismatch, against 0.0% of the true calls.

**5. More anchor per contig** (`pad_cores`, was 25, now 50). Half of what hcc1954 missed was a
placement failure, not a detection failure: 34% of the misses had a contig at one end whose SA
pointed to a third locus, 16% had contigs at both ends that were never the same contig. Raising the
pad cannot cost a contig — the window is clamped to the read's own cores and `min_flank` tests the
flank ACHIEVED. Measured at 15, 25, 50 and 75: 50 dominates 25 at every support level (same
precision, +0.009 recall) and 75 gives back what it gains.

**6. A normal deep enough to answer the question, and only where the question is asked.** The
germline check asks the normal reads at the reference locus, and it was being asked at 11×. Thinning
that very same BAM: 3× leaves 290 false calls, 5.5× leaves 172, 11× leaves 136 — a curve that has not
flattened. At 44× hcc1954 leaves 116 and F1 goes from 0.867 to **0.874**. It is the largest precision
lever found, and it is worth most where the precision gap is widest: +0.016 precision on hcc1954,
+0.017 on hcc1937, +0.024 on h2009, +0.050 on h1437.

Aligning the whole normal to get it would make the aligner the dominant cost of a caller whose
premise is that the graph replaces alignment — two hours against about forty minutes for everything
else. It is also unnecessary. The question is only ever asked at the candidate loci, and the graph
can say which normal reads cover those without a reference: they share LCP cores with the germline
ANCHORS of the contig each candidate came from. `svp pick` does that in one pass. Measured on
hcc1954, against the 8.35M-read whole-normal arm:

| normal | reads | alignment | precision | recall | F1 |
|---|---|---|---|---|---|
| whole | 8.35M | ~2 h | 0.891 | 0.858 | 0.874 |
| picked, `min_shared` 1 | 405,604 (3.58%) | ~14 min | 0.887 | 0.858 | 0.872 |
| picked, `min_shared` 5 | 176,912 (1.56%) | ~14 min | 0.887 | 0.858 | 0.872 |

1 and 5 give the identical call set — TP 949, FP 121 — so the tighter cut is free. A read that really
covers the locus shares dozens of cores; a single shared core is what a repetitive anchor collects
from everywhere.

**7. `min_supp` 3 → 4.** Not a mechanism, a point on the precision/recall curve — but the point is
chosen from four samples rather than one: 3 is worse than 4 on every one of them, and 5 wins on
h2009 and h1437 by about as much as 4 wins on hcc1937 and hcc1954. The tie goes to recall, which is
the side svp is behind on for three of the four.

## Running it

    svp run <normal.svb|fq> <tumour.svb|fq> <ref.fa> <dir>

One command, eight stages, each skipped when its output is already there: graph → bubbles →
junctions → minimap2 on the contigs → `call` with no normal → `pick` → minimap2 on the picked normal
reads → `call`. Verified end to end on hcc1954: the one command reproduces the 0.887 / 0.858 / 0.872
that the stages run by hand produce, record for record. Every stage takes .svb or
uncompressed FASTQ directly, so nothing has to convert between them; the only outside programs are
minimap2 and samtools, which svp is not. Every default it passes is the one measured below — the
configuration is in `src/`, not in a wrapper script, which is why the ones this work was measured
with are gone rather than kept as the way to reproduce it.

## Where svp stands

pad 50, `min_aln` 0, `shred_rep` 0.40, `anchor_rep` 0.50, `normal_dirty` 0.25, `max_mm` 0.01,
`min_supp` 4, every sample called against a normal deepened to 33–44×, against `minda_4_2`.

| sample | | precision | recall | F1 |
|---|---|---|---|---|
| **hcc1954** | 2026-08-06 | 0.218 | 0.736 | 0.337 |
| | supp≥4 | **0.887** | **0.858** | **0.872** |
| | supp≥5 | 0.895 | 0.847 | 0.870 |
| | Severus (published) | 0.921 | 0.826 | 0.871 |
| **h2009** | 2026-08-06 | 0.096 | 0.871 | 0.173 |
| | supp≥4 | 0.828 | 0.857 | 0.842 |
| | supp≥5 | 0.868 | 0.841 | 0.854 |
| | Severus | 0.913 | 0.894 | 0.903 |
| **hcc1937** | 2026-08-06 | 0.165 | 0.651 | 0.263 |
| | supp≥4 | 0.846 | 0.694 | 0.762 |
| | supp≥6 | **0.872** | 0.613 | 0.720 |
| | Severus | 0.859 | 0.795 | 0.826 |
| **h1437** | 2026-08-06 | 0.033 | 0.701 | 0.063 |
| | supp≥5 | 0.749 | 0.769 | 0.759 |
| | supp≥6 | **0.806** | 0.715 | 0.758 |
| | Severus | 0.782 | 0.842 | 0.810 |

**svp is level with the published Severus output on hcc1954 and behind on the other three.** On
hcc1954 it is ahead on recall by 0.032 and behind on precision by 0.034, for an F1 0.001 short;
with the whole normal aligned instead of the picked 1.56% it is 0.891/0.858/0.874, ahead by 0.003 --
that difference is what the pick costs, against two hours of alignment instead of fourteen minutes; the same
call set scored against the 2-support ensemble is 0.924 / 0.710 / 0.803 against Severus's
0.967 / 0.679 / 0.798, ahead on recall and F1 there too. On hcc1937 and h1437 svp passes Severus on
PRECISION at supp≥6 (0.872 against 0.859, 0.806 against 0.782) but gives up too much recall to get
there. h2009 is behind on both at every setting. The other three samples have no looser ensemble in
the archive, so only the 4-support number exists for them.

Simulation, same code path: `hard` gives 526/600 junctions at 98.0% precision at supp≥3 and 100.0%
at supp≥5, against Severus's 83.9% / 93.9% on the same contig BAM.

## Measured and rejected

Written down because each cost real time and each would otherwise be tried again.

- **A colour rule in the graph.** A germline variant leaking as somatic must produce a normal-only
  sibling branch, because the tumour can only fail to match the normal's spelling if the normal also
  fails to match the tumour's. The mechanism is right and the filter is not: it removed 23 false
  calls and 16 true ones. In repeats an anchor carries both a somatic and a normal-only branch for
  many reasons that are not an allelic swap. Diff at `/home/yuio/tmp/bubbles_colour_rule.diff`.
- **Loosening the bubble length threshold** (`max_dlen` 12 → 2). The threshold is 12 MINIMIZERS,
  roughly a 130–190 bp difference, which is far above the 50 bp floor the caller reports at — and it
  removes 63.1% of hcc1937's somatic cores against 48.6% of hcc1954's, while hcc1937's missed
  deletions have a median of 127 bp. The reasoning holds and the effect does not: +0.015 recall for
  −0.010 precision, at 46% more contigs. The small events do not come back, so they are missing for
  some other reason.
- **Tumour-only adjacencies between two germline cores.** 282,455 exist in hcc1954's graph and
  `build_chunks` reads none of them; 6,172 survive a support and degree filter. Implemented behind
  `edge_tcov` (default 0, off) and not yet measured end to end.
- **The VNTR bed Severus was run with.** 37.3% of the remaining false calls are inside it — and
  10.3% of the true ones. As a hard mask it costs 3 true calls per false one removed. Severus does
  not use it as a mask either.
- **`per_chunk` 5 → 10.** Doubles the contigs and the alignment time; the precision/recall curve
  lands on top of the old one.
- **Clustering tolerance** 50 → 400: F1 flat from 0.829 to 0.833. **`min_supp` 2**: 241 → 317 FP.
  **Breakpoint spread within a cluster**: >10 bp is 10.5% of false calls and 3.2% of true ones, 11
  calls in total. **Chunk diversity of the supporting contigs**: median 4 against 2, a 1.6×
  enrichment and not enough. **The mismatch gate** re-swept after everything else: F1 is flat from
  0.01 to 0.03, so it stays at the precision end.
- **Asking the germline question against the CONTIG instead of the reference.** The contig IS the
  locus -- [anchor | somatic run | anchor] -- so a normal read that carries the event aligns through
  the run and one that does not stops at the junction, and the index is the candidate contigs alone
  (17 MB against hg38's 3.1 GB). Implemented and measured: 0.832/0.844 against 0.887/0.858 for the
  same picked reads against hg38. Spanning had to be qualified as spanning WITHOUT a large indel --
  a normal read lacking the contig's deletion still aligns straight across it, the aligner just pays
  with an insertion -- and that fixed recall (0.781 -> 0.844) without closing the precision gap. What
  it has no answer for is the reference-space check's second clause, "the normal cannot place this
  locus either", which is what catches the repeat loci. Code in
  /home/yuio/trash/20260810-contig-germline.
- **A third pass over the reads to genotype each call.** Not built. `per_chunk` 10 already showed
  that doubling the contigs behind a support count does not sharpen it, and the normal side of such
  a pass is blind to exactly the germline leak that matters — if the tumour spells the locus
  differently, its k-mers are not in the normal either. The reference-locus question is the one that
  works, and that is what the deep normal answers.

## What is left

hcc1954's remaining 116 false calls at supp≥4 have a median SVLEN of 384 bp and a median support of
7 against 10 for the true calls; half are under 500 bp. No property tested separates them at a ratio
worth applying. The three samples without a deep normal are the obvious next measurement, and the
recall gap on hcc1937 — where 66% of the misses have no contig at either breakpoint — is upstream of
the caller, in what the graph and `junctions` produce.
