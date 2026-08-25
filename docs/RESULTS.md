# Results — everything that was measured

Scoring is **Minda**, the Severus paper's own instrument, at the paper's settings:

    minda truthset --base hifi/<pair>/truth_clean.vcf --vcfs <calls>.vcf --min_size 50 --tolerance 500

Severus is scored the same way, from its **published** output in the paper's archive, so the
comparison is symmetric. F1 against that truth, from the first real-data sitting to now:

| | hcc1954 | h2009 | hcc1937 | h1437 |
|---|---|---|---|---|
| 2026-08-06 | 0.337 | 0.173 | 0.263 | 0.063 |
| before the base-span bubble rule | 0.870 | 0.852 | 0.759 | 0.750 |
| **now** | **0.874** | **0.865** | **0.765** | **0.754** |
| Severus, published | 0.871 | 0.903 | 0.826 | 0.810 |

On the shipped ensemble snape passes Severus on hcc1954 and does not on the other three. On ensembles
rebuilt without each compared caller's own votes (see "Leave-one-caller-out" above) it is ahead on
hcc1954 (+0.018, p=4.8e-07), tied on h2009 (p=0.80) and h1437 (p=0.12), and behind on hcc1937
(-0.027, p=1.0e-07). Everything below is in the code unless it says otherwise.

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
`variant_calls_and_benchmarks.tar.gz`, the `minda_4_2` ensemble per cell line. The rule is in the
truth VCF's own `##cmd:` line — **support in at least 2 of the 3 platform groups AND at least 4 of
the 11 caller runs overall**. (This file used to say "at least 2 callers in each of 3 platform
groups", a stricter rule that only 25.8–69.7% of the records satisfy — 176 of hcc1954's 1106 records
have support in only two of the three groups, so they cannot have two in each of three. Note the
weaker point too: the rule above cannot be *confirmed* from the records, because every record
satisfies both clauses — none has fewer than four runs and none touches only one group. The `##cmd:`
line is present in two of the four truth files and states the same conditions in both.)
`truth_clean.vcf` is that file with the one unprefixed `cmd:` line Minda writes into its own output
removed, because pandas will not parse it.

**Leave-one-caller-out.** Four of the callers compared against are members of that ensemble, and a
Severus vote is in 93.9–97.3% of the records. `loo_truth.py` re-applies the rule with one caller's
runs deleted, which gives a truth that caller did not help build; snape is in no ensemble, so for a
snape-vs-X comparison both are scored on `truth_no_X` and neither side helped build the event set.
Doing that costs **every** ensemble member F1 on average over the four panels — severus −0.032,
nanomonsv −0.026, sniffles −0.021, savana −0.018 — which is what being inside an ensemble is worth,
and is the check that the construction is not chosen because it helps: snape is a member of none.
Against the debiased truths snape is ahead of nanomonsv by 0.041, savana by 0.031 and sniffles by
0.321, and behind severus by 0.018; on hcc1954 specifically the Severus comparison flips from −0.001
to **+0.018** (paired McNemar p = 4.8e-07).

(An earlier version of this paragraph said nanomonsv and savana *gained* 0.023 and 0.015. That came
from a bug: `loo_truth.py` had the SUPP_VEC letters for those two callers swapped, so each was being
scored against the truth the OTHER one had helped build. A caller gaining by the removal of its own
votes is impossible — deleting the records your votes carried can only remove events you found — and
that impossibility is what exposed it.)

**The scoring tolerance is not what decides it.** 500 bp is ten times the size floor and applies to
both breakpoints, so the ordering was re-checked at `--tolerance 100`. Both callers lose absolute F1
— snape 0.874 → 0.835 on hcc1954, Severus 0.871 → 0.831 — and the ordering is unchanged on every
panel, with snape's relative position slightly better at the tighter tolerance on all four:
+0.004 vs +0.003 (hcc1954), −0.032 vs −0.038 (h2009), −0.053 vs −0.060 (hcc1937), −0.052 vs −0.057
(h1437). Breakpoint precision is not where these callers differ.

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

**0. The bubble filter's length test, in bases rather than in minimizers** (`max_dbp`, was
`max_dlen` 12 minimizers, now 50 bp). A minimizer at k=21 w=22 over homopolymer-compressed bytes is
one sample per **16.6 raw bases** — measured, median over 175,753 bubbles on hcc1954, from the
`SNAPE_BUBBLE_TSV` dump. So `max_dlen = 12` was a ~200 bp test standing in front of a caller that
reports from 50 bp, and that is where the small insertions and deletions were going.

An arm's length in bases comes off its cores: consecutive level-5 cores of one read overlap by
**122.5 bp** (measured over 1.97e9 pairs on hcc1954 and 2.02e9 on h2009, and now re-measured on every
run rather than hard-coded — the simulated reads give 122.8), so an arm of n cores spans
`sum(len) - 122.5*(n-1)`. The threshold is not swept: it is `call`'s own `min_sv`. An arm with a core
no read gave bytes for has no measurable span and is kept rather than judged.

Calibrated offline against the truth loci — the dump plus `sigprobe`'s per-locus somatic signatures
score a candidate rule in two seconds instead of a 4-minute re-run. At equal recovery of
currently-missed events the base rule keeps far fewer cores than any setting of the minimizer rule:
h2009, 59 of 61 missed insertions for 145k cores against 294k; 27 of 27 missed deletions for 145k
against 547k.

End to end, each panel run twice from the same graph.bin with the same binary:

| panel | shipped P/R/F1 | new P/R/F1 | INS fn | DEL fn | BND fn | DUP fn |
|---|---|---|---|---|---|---|
| hcc1954 | 0.895/0.847/0.870 | 0.894/0.854/**0.874** | 12 → 10 | 26 → 20 | 118 = 118 | 13 = 13 |
| h2009 | 0.863/0.841/0.852 | 0.854/0.876/**0.865** | 63 → 33 | 48 → 46 | 34 = 34 | 2 = 2 |
| hcc1937 | 0.854/0.683/0.759 | 0.847/0.698/**0.765** | 55 → 50 | 77 → 71 | 80 = 80 | 11 = 11 |
| h1437 | 0.728/0.774/0.750 | 0.727/0.783/**0.754** | 8 → 8 | 16 → 14 | 23 = 23 | 3 = 3 |

**Breakend and duplication false negatives are identical on all four panels**, so the change moves
insertions and deletions and nothing else. Pooled: INS fn 138 → 101, DEL fn 167 → 151. Insertion
recall at 50–100 bp goes 49.6% → 71.3% (Severus 82.6%) and at 100–200 bp 77.8% → 83.7%; bands at and
above 200 bp do not move. Contigs cost +21% (h2009, 5.08 → 6.20 Gbp).

Simulation regression: junctions 526 → 527 of 600, insertions 70 → 73 of 78, precision unchanged at
99.7%. `hard` has no VNTRs, so a filter that had merely been loosened would show there as precision
lost for no recall gained.

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

The table above is the pre-fix build, which is where the ablation was run. On the shipped build the
same pair gives 1,909 candidates on 1,848 contigs, 128,219 anchor cores and 201,853 reads (1.78%),
and the pass takes 2:01. What the whole check is worth, re-measured on the shipped build at the
shipped floor: the contigs alone produce 1,692 calls at **0.560 / 0.857** (TP 948, FP 744), and
asking the normal at each of their loci removes 635 — 632 false, 3 true — for **0.894 / 0.854**
(TP 945, FP 112, 1,057 records). An earlier revision of this section quoted 1,601 → 554 → 0.895 /
0.847 from the pre-fix build and it did not reconcile with the 1,057 calls the shipped build
actually emits; the numbers above were re-run on 2026-08-13 and do.

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
| hcc1954 | 201,853 (1.78%) | 0.894 | 0.854 | **0.874** | 0.921 / 0.826 / 0.871 |
| h2009 | 346,852 (2.58%) | 0.854 | 0.876 | 0.865 | 0.913 / 0.894 / 0.903 |
| hcc1937 | 126,787 (1.17%) | 0.847 | 0.698 | 0.765 | 0.859 / 0.795 / 0.826 |
| h1437 | 189,710 (1.52%) | 0.727 | 0.783 | 0.754 | 0.782 / 0.842 / 0.810 |

For reference, the same four pairs on 2026-08-06 were 0.337 / 0.173 / 0.263 / 0.063 F1.

**snape is 0.003 of F1 ahead of the published Severus output on hcc1954 and behind on the other
three.** On hcc1954 it is ahead on recall by 0.028 for 0.027 of precision. These margins are small
relative to what this benchmark resolves: on ensembles rebuilt without each compared caller's own
votes, snape is ahead on hcc1954 (+0.018, paired McNemar p=4.8e-07), TIED on h2009 (p=0.80) and
h1437 (p=0.12), and behind on hcc1937 (-0.027, p=1.0e-07). The pipeline's own run-to-run spread IS measured, and it does not
reach the output: three complete replicates on hcc1954 from one graph.bin with one binary gave
identical chunk counts, 222,387 contigs each with 222,348 names in common (39 differ, 0.018%), and
**identical call sets** -- 1,692 records in all three, none appearing in only one. Those three
stopped before the germline check; two *complete* runs, germline check included, then picked 201,853
and 212,754 normal reads from the same 1,909 candidates -- a 5.4% difference -- and still produced
**1,057 identical calls each**. Replicated on all four pairs, ordered by cap saturation
(the share of populated chunks already holding five contigs, i.e. how often the race gets to choose):

| pair | cap saturation | contigs | names differing | calls | calls differing |
|---|---|---|---|---|---|
| hcc1954 | 20.7% | 222,387 | 39 (0.018%) | 1,057 | **0** |
| hcc1937 | 24.9% | 191,793 | 40 (0.021%) | 948 | **0** |
| h1437 | 50.3% | 298,768 | 151 (0.051%) | 1,073 | **0** |
| h2009 | 58.5% | 666,659 | 214 (0.032%) | 2,759 | **0** |

The intermediate representation moves; the call set does not, on any of the four. It does **not**
move in proportion to saturation — h1437 diverges most while h2009 is the most saturated. A fit on
the first three panels predicted 0.029% for h1437 and it came in at 0.051%; that prediction was
written down before the run, and there is no measured explanation for the miss. `junctions` does race for its per-chunk contig budget across threads,
but exchanging one of the five reads that report a junction does not change whether five reports
exist, and the germline check is a locus property that survives a 5% change in which normal reads it
is handed.

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

Simulation, whole pipeline at the shipped defaults. Two insilicoSV configurations on the same
reference, same binary, same scorer, same ±100 bp tolerance:

| | `hard` | `complex` |
|---|---|---|
| SVs ≥ 50 bp | 522 | 1,000 |
| event grammars | a handful | 23, incl. nested/composite (`delINVdel`, `dupINVdup`, `dDUP_iDEL`, `INV_rTRA`) |
| truth junctions | 600 | 1,602 |
| junctions recovered | **527 (87.8%)** | **1,358 (84.8%)** |
| insertions recovered | 73/78 (93.6%) | 123/130 (94.6%) |
| precision | 595/597 = **99.7%** | 1,444/1,464 = **98.6%** |

Per type on `hard`: DEL 90.6%, INV 92.3%, DUP 78.2%. On `complex` the same shape appears with 2.7×
the truth and far more resolution — **every weak class is a duplication**: DUP 64.6%, dDUP and
INV_dDUP 59.3%, INV_DUP 51.9%, and DUP at 50–1,000 bp 50.0%, against DEL 90.4% and INV 95.4%.
`complex` was run 2026-08-13 09:03 and nothing is tuned for it; `hard` is also the drop rule's
regression test, so it is not independent of the change and `complex` is.

(Before the change `hard` gave 526/600 and 70/78 at the same 99.7%.)

## Held out: two benchmarks outside CASTLE

The whole night's tuning happened on CASTLE, so at 04:55 the parameter vector and the binary's
sha256 were frozen and written to `/home/yuio/tmp/sv-callers/PARAMETER_FREEZE.md`, **before** either
of these samples' reads existed on this box in a runnable form. Nothing was changed afterwards on the
basis of either result. Both are PacBio Revio HiFi tumour/normal pairs, packed with `snape pack` and
run with `snape run` at the shipped defaults.

| benchmark | truth | events | P | R | F1 |
|---|---|---|---|---|---|
| HCC1395 / HCC1395BL | SVScope ensemble (Zenodo 18688406) | 1,118 | 0.802 | 0.810 | **0.806** |
| COLO829 / COLO829BL | Zenodo 7515830 v4, PCR/capture/Bionano validated | 62 | 0.563* | 0.726 | 0.634 |
| COLO829, Paulin 2025 stable subset | same, restricted to records stable over four stocks | 50 | 0.513* | **0.820** | 0.631 |

Depths measured from the packed bases: hcc1395 58x tumour / 41x normal, colo829 63x / 66x.

**The 32 PCR-validated records inside HCC1395's truth are the sharpest thing measured tonight.**
Scored as a truth set they are meaningless (every call outside them counts false); the number that
means something is recall, and minda matches **7 of the 32** — deletions 1/12, insertions 6/20.

What those 32 are is the point: **all 32 are TDscope calls and 26 of 32 fall inside a TRF tandem
array**. TDscope is tandem-repeat-specific, so this is not a random sample of the 1,118 — it is a
tandem-repeat event set, and the misses span 50–1,035 bp with a median of 100. That is largely the
population
`anchor_rep = 0.50` exists to suppress, and E09 measured the trade: relaxing it to 0.60 recovers nine
insertions on h2009 and takes false positives from 138 to 345. So: on orthogonally validated events
inside tandem repeats snape recovers a quarter, by a trade it makes deliberately. A platform artefact
was ruled out — 18 of the 32 followed hifi+ont-supported calls and 14 ont-only, and snape gets 4 of
each. It does not explain all of it: 10 of the 32 also carry a general-purpose caller's vote and
snape misses 6 of those 10.

*COLO829's precision is a LOWER BOUND and must not be compared with the CASTLE precisions. Its truth
is a validated subset — a list of events somebody confirmed by PCR, capture or optical mapping — not
a complete catalogue of the cell line's somatic SVs, so a call absent from it is not thereby false.
snape emits 80 calls against 62 truth records. Recall is what this pair supports, and there are two
because the cell line drifts: 0.726 against all 62, 0.820 against the 50 Paulin et al. found stable
across four replicate stocks.

COLO829's 62 events are also a small denominator — a 95% interval on a recall from 62 trials is about
±11 points — so it can show a gross failure and cannot resolve a small difference. Its value is that
its truth was built by orthogonal validation from platforms that do not include these reads, so no
caller in any comparison here helped construct it.

## Where the events die, and what the graph itself is worth

Two measurements that re-use work already on disk, and both should be repeated whenever a stage
changes.

**The ceiling partition.** Everything after minimap2 is `call`, which runs in seconds, so running it
on one aligned contig set at nested settings says exactly how many truth events the contigs contain
and what each filter removes. `ceiling.sh` + `ceiling_report.py` do this in 66 seconds. h2009,
shipped rule, of 572 truth insertions:

| where an insertion is lost | count |
|---|---|
| never reaches the aligned contigs | 29 |
| the support floor `min_supp = 5` | 25 |
| the repeat guards | 9 |
| the germline check | 0 |

The four sum to the panel's 63 insertion false negatives. Of the 25 the floor discards, **eleven sit
at support exactly 4** while `junctions` emits at most 5 contigs per chunk — but only 11 of the 32
floor-failing INS/DEL events are in a chunk that actually reached the cap, so raising `per_chunk`
cannot reach the rest (see "Measured and rejected").

**What the core graph contributes.** Run the identical caller — same code, same thresholds — over
whole-read alignments instead of junction contigs, with the germline check reading the WHOLE normal
alignment rather than the selected 1.6%:

| hcc1954 | precision | recall | F1 | TP | FP |
|---|---|---|---|---|---|
| whole reads, no graph | 0.635 | 0.873 | 0.735 | 966 | 556 |
| the graph pipeline | 0.894 | 0.854 | 0.874 | 945 | 112 |

**The graph is worth +0.259 of precision and costs −0.019 of recall — as an UPPER BOUND**: 444 false calls removed, 21
true ones lost. For scale, the germline check against the normal is worth more still (+0.334,
0.560 → 0.894), so the graph is the second-largest precision instrument in the pipeline, not the
largest. It is the largest single accuracy instrument in the pipeline, larger than any
threshold in `call` — and it is not free. Only 23 insertions survive in the whole-read arm against
530 in the graph pipeline.

That is a bound and not a matched comparison: `min_supp = 5` counts distinct **contigs, capped at five per chunk** on the graph side and distinct **reads at 60× depth, uncapped** on the whole-read side, so part of the precision gap is the cap rather than the graph. Two further asymmetries run the same way — the whole-read arm gets the whole normal for its germline check, and the two arms cluster on different units. `experiments/10_whole_read_ablation.md` recorded this in its Limits; the round-9 referee (2026-08-13) caught that the paper and this file had not. A matched-semantics re-run would turn the bound into a measurement and has not been done.

## What it costs, and what the other arm costs

Per stage on hcc1954, 16 threads pinned, machine otherwise idle, at the shipped defaults. The snape
arm was re-measured end to end on 2026-08-13 with the shipped binary; the alignment-first arm is the
2026-08-11 measurement on the same machine and the same `.svb` inputs, and involves no snape code.

| arm | stage | wall | peak RSS |
|---|---|---|---|
| **snape** | `graph` (both read sets, 391 Gbp) | 10:21 | 32.8 GB |
| | `bubbles` | 04:02 | 2.9 GB |
| | `junctions` | 02:06 | 2.9 GB |
| | minimap2 + sort + index, 222,387 contigs / 2.04 Gbp | 04:00 | 19.1 GB |
| | `call` (candidates, floor 4) | 00:02 | 0.04 GB |
| | `pick` | 02:01 | 0.12 GB |
| | minimap2 + sort + index, 212,754 normal reads / 3.79 Gbp | 16:25 | 20.9 GB |
| | `call` | 00:05 | 0.04 GB |
| | **total** | **39:01** | **32.8 GB** |
| **alignment-first** | minimap2 + sort + index, 196 Gbp tumour | 2:02:57 | 20.4 GB |
| | minimap2 + sort + index, 195 Gbp normal | 1:54:23 | 19.9 GB |
| | Severus 1.7, both BAMs, 16 threads | 09:47 | 13.2 GB |
| | **total** | **4:07:07** | **20.4 GB** |

**6.3× less wall time, at 1.6× the peak memory**, and 1.1 GB of alignments left on disk against
93 GB. Severus itself is not the cost — it runs in 9:47; the 3:57 of alignment in front of it is.
Both arms were fed from the same `.svb`, so neither paid for gzip; reading the original `.gz` would
have made the alignment arm slower still.

Against the 2026-08-11 table this replaces, three things moved and one prediction failed. `graph`
13:36 → 10:21, `bubbles` 9:43 → 4:02 and `pick` 4:59 → 2:01, all from a reader/worker split committed
on 08-12 that has nothing to do with the small-event work — **the tool was faster than its own
documentation for a day and a half**. The prediction that 10% more contigs would cost the contig
alignment ~25 s failed: that stage is 3:59.77 against 4:04, i.e. unchanged within noise. The one real
increase is the picked-normal alignment, 11:25 → 16:24, of which `pick` selecting ~20% more reads
accounts for about half; the rest is unexplained and was not chased.

`pick`'s output is the one quantity that varies between runs: 201,853 reads (1.78%) in one complete
run and 212,754 (1.88%) in another, from the same 1,909 candidates, because the ~39 contigs in
222,387 that differ run to run carry different anchor sequence. **The final call set does not vary
with it — 1,057 calls in both, identical in position and allele.**

Severus on those BAMs (1.7, unphased, its own VNTR bed) scores **0.912 / 0.819 / 0.863**, against
the 0.921 / 0.826 / 0.871 of its published calls — the input path is not what separates the two
callers.

**None of these callers reads base qualities.** Severus computes its "read quality" from the
alignment (NM tag, CIGAR mismatch rate, MAPQ, aligned length), and `query_qualities` never appears
in its source; pbsv and sawfish document no base-quality option at all; minimap2 scores from
sequence only. That is why the `.svb` path — 2-bit, no qualities, no read names — costs nothing
here, for us or for them.

## Running it

    snape run --samples <sheet.tsv> <ref.fa> <dir>

The sheet is `name  role  path`, one row per read set, exactly one `normal` and at least one
`tumor`. One command, eight stages, each skipped when its output is already there: graph → bubbles →
junctions → minimap2 on each line's contigs → `call` with no normal → `pick` → minimap2 on the
picked normal reads → `call`. Verified end to end on hcc1954 (one tumour) and on a synthetic
1-normal/2-tumour case, where each line's deletion appears in that line's column and nowhere else.
Every stage takes .svb or uncompressed FASTQ directly, so nothing has to convert between them; the
only outside programs are minimap2 and samtools.

## Measured and rejected

Written down because each cost real time and each would otherwise be tried again.

- **Tumour-only adjacencies between two germline cores.** This entry used to say the idea was
  "implemented and never measured end to end". It has now been measured and the seed is not there:
  `sigprobe` over 1,462 truth loci on four panels finds such an edge at 0.7-3.0% of loci, no more
  often at missed events than at found ones (Fisher p = 0.7 INS, 0.79 DEL), and **not one** of the
  141 loci where no somatic core forms carries one. The geometry explains it: for a junction to leave
  no novel core, its breakpoint has to land on a core boundary, and a core is 217 bp. The loci that
  produce no somatic core are large deletions instead, median 3.7 kb, where there is no recall gap.
- **Keeping the bubble filter's cores fuses chunks and starves `per_chunk`.** It does not.
  `chunkstat` reproduces `build_chunks` over hcc1954's graph; with NO filter at all the component
  size distribution is essentially unchanged (median 1 -> 2 cores, q0.999 10 -> 8, max 3479 either
  way) and the share of somatic cores in a component of >=10 cores FALLS from 6.9% to 3.6%.
- **Raising `per_chunk`**, re-examined with a mechanism rather than a sweep. `junctions` writes each
  contig's chunk id into its FASTA header and `call` writes a supporting contig into every record, so
  the chunk behind each discarded event is recoverable. Of the 32 truth INS/DEL events that fail
  `min_supp = 5` on h2009, only 11 sit in a chunk that reached the cap of 5; the other 21 are in
  chunks fewer than five reads carried at all. That is why "5 -> 10 lands on top of the old curve".
- **Node support as a bubble-drop criterion.** Never dropping a bubble the tumour traversed >= T
  times is dominated by the rule it was meant to improve on: at its best setting it recovers all 61
  of h2009's currently-missed insertions for 558,714 surviving cores, while simply lowering the OLD
  minimizer rule to `dlen <= 2` recovers the same 61 for 547,180. Length is the discriminator; depth
  is not.
- **Tightening the colinearity clause** from 20 minimizers (~330 bp) towards the 50 bp floor, which
  is the only protection a balanced junction has once the length clause is in bases. Calibrated
  against all five SV types at 1,917 hcc1954 breakpoint loci: it recovers 1 breakend, 1 duplication
  and 6 deletion loci for **3.2x** the surviving cores. Left at 20.
- **Breaking a near-tie in `tandem()`** towards a fixed frame, so an insertion in a tandem array does
  not split its support across two clusters `len` apart. h2009 at the shipped guards: TP 812 either
  way, one extra false call.
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
