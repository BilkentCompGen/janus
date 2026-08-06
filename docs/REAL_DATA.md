# svp on real tumour/normal HiFi — what was measured, 2026-08-05/06

The simulation numbers do not survive contact with real data. This records the measurements so the
next attempt starts from evidence instead of from the `hard` result. **Nothing here is in the code:**
every intervention below was tried, measured, and reverted. The working tree at `b34dbcd` is the
state these experiments started from and the state they were rolled back to.

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
`variant_calls_and_benchmarks.tar.gz`, the `minda_4_2` ensemble per cell line. It is a consensus
call set, not experimental truth: built by `minda ensemble --min_size 50 --tolerance 500 --conditions
"[[[caller_names[:4], caller_names[4:8], caller_names[8:11]], '>=', 2],'&', [caller_names,'>=', 4]]"`
— at least 2 callers in each of 3 platform groups AND at least 4 callers overall. Copied to
`hifi/<sample>/truth_minda_4_2.vcf`; `truth_clean.vcf` beside it is the same file with the one
unprefixed `cmd:` line Minda writes into its own output removed, because pandas will not parse it.
Coordinates are **chr-prefixed hg38**.

**Scoring** — Minda, the paper's instrument, at the paper's settings:
`minda truthset --base <truth> --vcfs <calls> --min_size 50 --tolerance 500`.

## Where svp stands

Full pipeline (`graph → bubbles → junctions → minimap2 → call`), against the paper's own Severus
output on the same truth:

| sample | | TP | FN | FP | precision | recall | F1 |
|---|---|---|---|---|---|---|---|
| hcc1937 | Severus (reads) | 559 | 144 | 92 | 0.859 | 0.795 | 0.826 |
| | svp | 458 | 245 | 2323 | 0.165 | 0.651 | 0.263 |
| hcc1954 | Severus (reads) | 914 | 192 | 78 | 0.921 | 0.826 | 0.871 |
| | svp | 814 | 292 | 2913 | 0.218 | 0.736 | 0.337 |
| h2009 | Severus (reads) | 827 | 98 | 79 | 0.913 | 0.894 | 0.903 |
| | svp | 806 | 119 | 7611 | 0.096 | 0.871 | 0.173 |
| h1437 | Severus (reads) | 186 | 35 | 52 | 0.782 | 0.842 | 0.810 |
| | svp | 155 | 66 | 4556 | 0.033 | 0.701 | 0.063 |

**Recall is competitive; precision is not.** On h2009 recall is 0.871 against Severus's 0.894.

Timings, per sample, 16 cores: graph 12–14 min (peak RSS 34–39 GB), bubbles ~10 min, junctions
4–6 min, minimap2 ~2 min, call < 1 s. About 30 minutes end to end — the OOM that was feared did not
happen, a real pair is 2.2× the simulation's read volume but peaks below the simulation's 52 GB.

## Splitting the blame: same contigs, two callers

Severus run on **our** contig BAM, so the only variable is the caller:

| sample | Severus / our contigs | svp call / our contigs | Severus / reads |
|---|---|---|---|
| hcc1937 | 0.305 prec, 0.738 rec | 0.165, 0.651 | 0.859, 0.795 |
| hcc1954 | 0.432, 0.882 | 0.218, 0.736 | 0.921, 0.826 |
| h2009 | 0.232, 0.878 | 0.096, 0.871 | 0.913, 0.894 |
| h1437 | 0.099, 0.801 | 0.033, 0.701 | 0.782, 0.842 |

Two facts, both consistent across four samples:

1. **Our caller is worse than Severus on identical input** — roughly half the precision at lower
   recall. Something is missing in `svp call`, not only upstream.
2. **The contigs cost Severus itself ~0.5 precision.** Severus drops from 0.86–0.92 on reads to
   0.10–0.43 on our contigs. That is the larger half of the gap and it is upstream of the caller.
3. On hcc1954 Severus recovers **more** truth from our contigs (0.882) than from the reads (0.826) —
   the contigs carry more signal than the read pileup, and drown it.

## What the false positives are

Measured on hcc1954 (814 TP / 2913 FP), each against the TP set as control:

| property | TP | FP |
|---|---|---|
| median alignment mismatch rate at the locus | **0.27%** | **15.62%** |
| tandem-repeat score of ±150 bp (best period 1–6; random DNA ≈ 0.25) | 0.326 | 0.466 |
| repeat score > 0.5 | **0.4%** | **41.7%** |
| supporting contigs, median | 9 | 5 |
| supporting contigs ≤ 4 | 3.1% | 43.5% |
| median SVLEN | 63,583 bp | **172 bp** |

Best-period breakdown of the repetitive FP loci: p2 31%, p5 37%, p4 11%, p6 5%, p3 4%, **p1 0%** —
the homopolymer collapse in `extract()` is already doing its job, and what is left is period ≥ 2.

Type composition, hcc1954: we emit 1232 DEL / 1150 INS against Severus's 145 / 33 and a truth of
186 / 20. FPs are 74% DEL+INS. Note hcc1954's truth is dominated by large events (median SVLEN
497 kb, only 10% under 500 bp) so small-event FPs there are partly unmeasurable — but **h2009's truth
is 46% under 500 bp and its precision is the worst of the four (0.096)**, so truth coverage does not
explain it.

## Interventions tried, all reverted

| # | intervention | result |
|---|---|---|
| 1 | `min_vaf` on chunks (already in `junctions`, never used) | hcc1954 F1 0.337 → 0.399 at 0.05, then **down**: 0.385 at 0.10, 0.271 at 0.20. Weak lever, dominated by everything else. |
| 2 | tandem-repeat compression, `SVP_HPC_PERIOD` = 2 / 6 in `extract()` | **negative.** hcc1937 F1 0.390 → 0.383 → 0.378. Nodes fall 7% but contigs **triple** (168k → 497k): a core spans more raw bases in a repeat, so tumour/normal differences spread over a wider window and produce *more* tumour-only cores. Sequence cost is small (hg38 collapses 29.8% at period 1, 33.6% at 3, 34.0% at 6) — the idea fails on its own terms, not on cost. |
| 3 | mismatch-rate gate in `call` (`NM` / aligned length) | **the one large effect.** Fixed cut at 1%: hcc1954 F1 0.335 → **0.692** (beats Severus-on-our-contigs, 0.580); h1437 0.063 → **0.483**; hcc1937 0.263 → 0.390 (0.426 at 2%); h2009 0.172 → 0.258. Costs real recall. |
| 3b | Severus's own adaptive version of the same gate (`background_mm_rat` + `background_mm_hist`: 95th percentile, fired only in windows where clean and divergent alignments coexist) | ported faithfully, **underperforms the fixed cut**: hcc1954 F1 0.387 vs 0.692. Its 95th percentile lands at **10.6%** on our contigs and flags only 433 of 206,285 windows — the rule assumes a read population, and our contig population is far dirtier. |
| 4 | carry the graph's read counts into the caller — `junctions` writes `>read_n_c<chunk>_t<tumour reads>_d<depth>`, `call` gates on reads instead of counting contigs | **negative, and it was my "fundamental defect" claim.** At matched recall it is *worse* than counting contigs (hcc1937, no mm gate: contigs ≥3 → 0.165 prec / 0.651 rec; reads ≥3 → 0.063 / 0.649). Only ≥20 reads gains anything (F1 0.263 → 0.315) and it buys that with recall. |

The motivation for #4 still stands as an observation even though the fix did nothing: `per_chunk=5`
caps contigs per chunk and 15,122 of hcc1937's chunks sit exactly on that cap, so a contig count
cannot separate a junction carried by 60 reads from one carried by 3. It simply turns out that read
count does not separate true from false either.

## The open question

**1,071 false positives survive at ≥ 20 supporting tumour reads.** These are not thin-evidence
artefacts: the sequence is in 20+ tumour reads and absent from a 60× normal, so sampling cannot
explain them. Either they are real somatic events the 4-caller consensus does not contain, or the
exact-string somatic test (`ncov == 0` on a 217 bp core) calls a difference where the normal merely
carries a different length variant of the same repeat.

The experiment that separates those two, not yet run: take the somatic core signatures behind the
well-supported FPs and **scan the normal reads for them** — `svp bubbles`' `collect()` already does
exactly this kind of scan. Found in the normal reads ⇒ the graph's `ncov` is leaking (min-mult floor
or the bloom cascade) and the somatic definition is the defect. Not found ⇒ the differences are real
and the ceiling is the truth set, which changes what "precision" on this benchmark even means.

## Reproducing / picking this up

- Batch: `hifi/run_svp.sh` (five pairs, strictly sequential, every stage skipped if its output exists).
- The experiment scripts are `hifi/vaf_test.sh`, `hifi/hpc_test.sh`, `hifi/depth_test.sh`.
- The reverted diff for interventions 2–4 is at `/home/yuio/tmp/svp-realdata-experiments.diff`.
- Outputs under `/home/yuio/tmp/jx/<sample>/`, graphs under `outputs/svp/<sample>/`.
