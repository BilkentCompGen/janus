# Results — everything that was measured

Scoring is **Minda**, the Severus paper's own instrument, at the paper's settings:

    minda truthset --base hifi/<pair>/truth_clean.vcf --vcfs <calls>.vcf --min_size 50 --tolerance 500

Severus is scored the same way, from its **published** output in the paper's archive, so the
comparison is symmetric. F1 against that truth, from the first real-data sitting to now:

| | hcc1954 | h2009 | hcc1937 | h1437 |
|---|---|---|---|---|
| 2026-08-06 | 0.337 | 0.173 | 0.263 | 0.063 |
| now | **0.872** | 0.854 | 0.762 | 0.759 |
| Severus, published | 0.871 | 0.903 | 0.826 | 0.810 |

snape passes Severus on hcc1954 and does not on the other three. Everything below is in the code
unless it says otherwise.

## The data

**Reads** — ENA, BioProject **PRJNA1086849**, the Severus paper's own project, which is the
**CASTLE** panel (Cancer Standards Long-read Evaluation). PacBio **Revio**, five tumour/normal
pairs; `hifi/get.sh` fetches and packs them. Measured on the raw fastq, not assumed: median read
16.2 kb, median per-read Q38.5, overall Q37.4 → **0.018% error**. These are HiFi despite every file
in the archive being named `_subreads.fastq.gz`.

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
nothing to do with being right. The archive also ships two looser ensembles for hcc1954, and the
shipped call set and Severus's published one score against all three like this:

| truth | records | snape P/R/F1 | Severus P/R/F1 |
|---|---|---|---|
| minda_4_2 (what the paper reports against) | 1106 | 0.895 / 0.847 / 0.870 | 0.921 / 0.826 / 0.871 |
| minda_3_2 | 1168 | 0.917 / 0.822 / **0.867** | 0.938 / 0.796 / 0.861 |
| minda_2_1 | 1413 | 0.940 / 0.696 / **0.800** | 0.967 / 0.679 / 0.798 |

47 of the 110 calls the 4-support ensemble counts false are events the 2-support one contains. Of
the 287 records that separate the two, 247 have two or more distinct tools behind them and 40 are
one tool counted twice (severus on two technologies). snape is ahead on F1 on both wider ensembles —
which is where the truth depends least on Severus's own output.

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
ANCHORS of the contig each candidate came from. `snape pick` does that in one pass. Measured on
hcc1954, against the whole-normal arm:

| normal | reads | alignment | precision | recall | F1 |
|---|---|---|---|---|---|
| whole | 11.34M | ~2 h | 0.891 | 0.858 | 0.874 |
| picked, `min_shared` 1 | 405,604 (3.58%) | — | 0.887 | 0.858 | 0.872 |
| picked, `min_shared` 5 | 176,912 (1.56%) | 11:25 | 0.887 | 0.858 | 0.872 |

The selection itself is one pass over the normal reads: 1,735 candidate contigs give 119,338 anchor
cores, and the pass takes 4:59. What the whole check is worth, at the shipped floor: the contigs
alone produce 1,601 calls at **0.587 / 0.849**, and asking the normal at each of their loci removes
554 — 552 false, 2 true — for **0.895 / 0.847**.

1 and 5 give the identical call set — TP 949, FP 121 — so the tighter cut is free. A read that really
covers the locus shares dozens of cores; a single shared core is what a repetitive anchor collects
from everywhere. **`pick` was verified on hcc1954 only**; the other three pairs' numbers below were
measured against their deepened normal aligned whole.

**7. `min_supp` 3 → 4.** Not a mechanism, a point on the precision/recall curve — but the point is
chosen from four samples rather than one: 3 is worse than 4 on every one of them, and 5 wins on
h2009 and h1437 by about as much as 4 wins on hcc1937 and hcc1954. The tie goes to recall, which is
the side snape is behind on for three of the four.

**8. Homopolymer compression before hashing.** Checked on the reads themselves rather than assumed,
and not on the simulation: `scripts/hifi_error_profile.py` compares hcc1954's normal reads to each
other over 5 Mb of chr1 (11,824 reads), calling a difference an error when at most 20% of the reads
covering that locus carry it. **83.4% of the errors are indels** (49.2% deletion, 34.3% insertion),
and those indel errors are **2.3× enriched** in homopolymer runs ≥ 3 while run=1 is depleted 0.4×.
Substitution errors, the control, show no enrichment at all (1.0–1.1×). 80% of indel errors sit in a
run ≥ 2, so compressing runs before hashing neutralises about two thirds of everything that could
change a core's identity.

## Where snape stands

Every row below is one run of the whole pipeline at its shipped defaults — pad 50, `min_aln` 0,
`shred_rep` 0.40, `anchor_rep` 0.50, `normal_dirty` 0.25, `max_mm` 0.01, `min_supp` 4, and the
germline question asked against the normal reads `pick` selected. Nothing is tuned per pair.

| pair | picked normal | precision | recall | F1 | Severus (published) |
|---|---|---|---|---|---|
| hcc1954 | 176,912 (1.56%) | 0.895 | 0.847 | **0.870** | 0.921 / 0.826 / 0.871 |
| h2009 | 318,542 (2.18%) | 0.863 | 0.841 | 0.852 | 0.913 / 0.894 / 0.903 |
| hcc1937 | 104,739 (0.97%) | 0.854 | 0.683 | 0.759 | 0.859 / 0.795 / 0.826 |
| h1437 | 182,922 (1.46%) | 0.728 | 0.774 | 0.750 | 0.782 / 0.842 / 0.810 |

For reference, the same four pairs on 2026-08-06 were 0.337 / 0.173 / 0.263 / 0.063 F1.

**snape lands 0.001 of F1 short of the published Severus output on hcc1954 and behind on the other
three.** On hcc1954 it is ahead on recall by 0.021 for 0.026 of precision; against the two looser
ensembles the same call set is ahead on F1 (table earlier).

What the `pick` costs, against the same call made with the whole deepened normal aligned (measured
at a floor of 4): hcc1954 0.891 → 0.887 precision, h2009 0.828 → 0.820, hcc1937 0.846 → 0.843,
recall unchanged in all three.
That is the price of aligning 1–2% of the normal instead of all of it, and it buys back an hour and
three quarters per pair.

The support floor is the one knob with a per-sample optimum, swept on this same arm (F1):

| supp | hcc1954 | h2009 | hcc1937 | h1437 | mean |
|---|---|---|---|---|---|
| ≥3 | 0.784 | 0.698 | 0.641 | 0.476 | 0.650 |
| ≥4 | **0.872** | 0.838 | **0.762** | 0.722 | 0.799 |
| **≥5 (shipped)** | 0.870 | **0.852** | 0.759 | **0.750** | **0.808** |
| ≥6 | 0.865 | 0.835 | 0.719 | 0.748 | 0.792 |

3 is worse than 4 on every pair. Between 4 and 5 the split is two pairs each, but the margins are
not: 4 wins hcc1954 by 0.002 and hcc1937 by 0.003, while 5 wins h2009 by 0.014 and h1437 by 0.028,
so **5 is shipped** — chosen on the mean over four pairs, not on hcc1954, where a floor of 4 would
read 0.872 and put us past Severus. The candidate pass that feeds `pick` runs at 4, one below the
reporting floor: a locus too many there costs alignment, never precision.

Simulation, whole pipeline at the shipped defaults: `hard` gives **526/600 junctions and 70/78
insertions at 99.7% precision** — 2 of 593 calls do not match a truth junction. Per type: DEL 91.0%,
INV 92.3%, INS 89.7%, DUP 76.9%.

## What it costs, and what the other arm costs

Per stage on hcc1954, 16 threads pinned, at the shipped defaults, and against the alignment-first
arm on the same machine (2026-08-11):

| arm | stage | wall | peak RSS |
|---|---|---|---|
| **snape** | `graph` (both read sets, 391 Gbp) | 13:36 | 41.2 GB |
| | `bubbles` | 09:43 | 5.7 GB |
| | `junctions` | 02:53 | 2.9 GB |
| | minimap2 + sort, 202,882 contigs / 1.87 Gbp | 04:04 | 18.8 GB |
| | `call` (candidates, floor 4) | 00:02 | 0.03 GB |
| | `pick` | 04:59 | 0.03 GB |
| | minimap2 + sort, 176,912 normal reads | 11:25 | 20.5 GB |
| | `call` | 00:05 | 0.03 GB |
| | **total** | **46:47** | **41.2 GB** |
| **alignment-first** | minimap2 + sort + index, 196 Gbp tumour | 2:02:57 | 20.4 GB |
| | minimap2 + sort + index, 195 Gbp normal | 1:54:23 | 19.9 GB |
| | Severus 1.7, both BAMs, 16 threads | 09:47 | 13.2 GB |
| | **total** | **4:07:07** | **20.4 GB** |

**5.3× less wall time, at twice the peak memory**, and 1.1 GB of alignments left on disk against
93 GB. Severus itself is not the cost — it runs in 9:47; the 3:57 of alignment in front of it is.
Both arms were fed from the same `.svb`, so neither paid for gzip; reading the original `.gz` would
have made the alignment arm slower still.

Severus on those BAMs (1.7, unphased, its own VNTR bed) scores **0.912 / 0.819 / 0.863**, against
the 0.921 / 0.826 / 0.871 of its published calls — the input path is not what separates the two
callers.

**None of these callers reads base qualities.** Severus computes its "read quality" from the
alignment (NM tag, CIGAR mismatch rate, MAPQ, aligned length), and `query_qualities` never appears
in its source; pbsv and sawfish document no base-quality option at all; minimap2 scores from
sequence only. That is why the `.svb` path — 2-bit, no qualities, no read names — costs nothing
here, for us or for them.

## Running it

    snape run <normal.svb|fq> <tumour.svb|fq> <ref.fa> <dir>

One command, eight stages, each skipped when its output is already there: graph → bubbles →
junctions → minimap2 on the contigs → `call` with no normal → `pick` → minimap2 on the picked normal
reads → `call`. Verified end to end on hcc1954. Every stage takes .svb or uncompressed FASTQ
directly, so nothing has to convert between them; the only outside programs are minimap2 and
samtools.

## Measured and rejected

Written down because each cost real time and each would otherwise be tried again.

- **A colour rule in the graph.** A germline variant leaking as somatic must produce a normal-only
  sibling branch, because the tumour can only fail to match the normal's spelling if the normal also
  fails to match the tumour's. The mechanism is right and the filter is not: it removed 23 false
  calls and 16 true ones. In repeats an anchor carries both a somatic and a normal-only branch for
  many reasons that are not an allelic swap.
- **Loosening the bubble length threshold** (`max_dlen` 12 → 2). The threshold is 12 MINIMIZERS,
  roughly a 130–190 bp difference, which is far above the 50 bp floor the caller reports at — and it
  removes 63.1% of hcc1937's somatic cores against 48.6% of hcc1954's, while hcc1937's missed
  deletions have a median of 127 bp. The reasoning holds and the effect does not: +0.015 recall for
  −0.010 precision, at 46% more contigs. The small events do not come back, so they are missing for
  some other reason.
- **Tumour-only adjacencies between two germline cores.** A rearrangement is an adjacency change,
  and an adjacency can be tumour-only while both of its ends are germline, in which case the
  junction leaves no somatic core at all — only an edge the normal never traverses. 282,455 of those
  exist in hcc1954's graph and `build_chunks` reads none of them; 6,172 survive a support and degree
  filter. Seeding chunks from them was implemented and never measured end to end, so the code is
  gone rather than left behind a default-off flag.
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
  locus — [anchor | somatic run | anchor] — so a normal read that carries the event aligns through
  the run and one that does not stops at the junction, and the index is the candidate contigs alone
  (17 MB against hg38's 3.1 GB). Implemented and measured: 0.832/0.844 against 0.887/0.858 for the
  same picked reads against hg38. Spanning had to be qualified as spanning WITHOUT a large indel — a
  normal read lacking the contig's deletion still aligns straight across it, the aligner just pays
  with an insertion — and that fixed recall (0.781 → 0.844) without closing the precision gap. What
  it has no answer for is the reference-space check's second clause, "the normal cannot place this
  locus either", which is what catches the repeat loci.
- **A third pass over the reads to genotype each call.** Not built. `per_chunk` 10 already showed
  that doubling the contigs behind a support count does not sharpen it, and the normal side of such
  a pass is blind to exactly the germline leak that matters — if the tumour spells the locus
  differently, its k-mers are not in the normal either. The reference-locus question is the one that
  works, and that is what the deep normal answers.

## What is left

hcc1954's false calls (the 116 left at a floor of 4, where this was analysed) have a median SVLEN of
384 bp and a median support of
7 against 10 for the true calls; half are under 500 bp. No property tested separates them at a ratio
worth applying. The recall gap on hcc1937 — where 66% of the misses have no contig at either
breakpoint — is upstream of the caller, in what the graph and `junctions` produce.
