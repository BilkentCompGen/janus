<!-- Deep-research report (workflow w192pbfvx, 2026-06-08).
Headline: Target ~200 somatic SVs/genome (a typical metastatic solid-tumor median): ~96 DEL, ~52 DUP, ~26 INV, ~26 INS, 0 translocations (out of scope for svp). -->

# Realistic per-tumor-genome somatic SV counts for one simulated WGS tumor (insilicoSV)

## 1. The methodological pitfall: never divide a cohort union by sample count

The single most common error in sizing a simulated tumor is taking a published **cohort total of unique SVs** and dividing it by the number of donors. That `union / N` quotient is **not** a per-genome count, because somatic (and germline) SVs are overwhelmingly *private*: each new genome added to a cohort contributes mostly novel breakpoint junctions, so the discovery union keeps growing with cohort size and composition rather than converging to a per-genome value.

The germline data proves this decisively. gnomAD-SV's discovery union grew from **433,371 SVs across 14,891 genomes** to **1,199,117 SVs across 63,046 genomes** — a 2.77x jump in the union for a 4.23x jump in sample count — while the **per-genome median barely moved: 7,439 -> 11,844 SVs/genome** ([Collins et al. 2020, PMID 32461652](https://pubmed.ncbi.nlm.nih.gov/32461652/); [gnomAD-SV v4, 2023](https://gnomad.broadinstitute.org/news/2023-11-v4-structural-variants/)). Dividing the v2 union by its sample count gives 433,371 / 14,891 ≈ 29 — off from the true per-genome median of ~7,400 by ~250x. The union is a property of the *cohort*, not of any genome in it.

The same caution applies on the somatic side. PCAWG is a cohort-aggregate union: **288,416 somatic SVs across 2,583 white-listed donors** ([PCAWG flagship, Nature 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7025898/)). Unlike gnomAD's *discovery union of distinct sites*, this PCAWG number is the total count of SV calls summed over donors, so the per-donor **mean** 288,416 / 2,583 ≈ 112 is a defensible average burden — *not* a union/N artifact. (The existing `simulation/configs/wgs.yaml` correctly derives 112 this way, as a per-donor mean alongside ~16,950 SNVs and ~936 indels.) The problem with 112 is therefore not that it is mathematically wrong, but that it is (a) a **mean**, which is inflated and unstable for the right-skewed SV distribution, and (b) drawn from **primary** tumors, which carry roughly half the SV burden of metastatic tumors.

**The correct statistic is the per-genome MEDIAN number of SV events, stratified by tumor type, with a range.** A median is robust to the heavy right tail (chromothriptic and hypermutated genomes); type stratification matters because the between-type spread is ~7-fold and a single pooled number misrepresents any individual tumor.

## 2. The anchor: per-genome median, by cancer type

The cleanest per-genome somatic anchor is the **Hartwig metastatic pan-cancer cohort**: a **median of 193 SVs per tumour** across **2,520 tumour-normal pairs** (cohort total 653,452 SVs) ([Priestley et al., Nature 2019, PMC6872491](https://pmc.ncbi.nlm.nih.gov/articles/PMC6872491/)). This is explicitly reported as a per-tumour median, not an aggregate (the naive pooled mean would be 653,452 / 2,520 ≈ 259, larger than the median exactly because the distribution is right-skewed).

Per-cancer-type medians span a **~7-fold range (≈56–412)**, which is the quantitative reason a single pooled number is misleading ([Priestley et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6872491/)):

| Profile | Cancer types | Median SVs/genome |
|---|---|---|
| SV-quiet (low) | neuroendocrine tumours (56), kidney (71) | **~56–71** |
| Typical solid tumor | cohort median | **193** |
| SV-driven (high) | oesophageal (372), ovarian (412) | **~372–412** |

PCAWG *primary* tumors run lower (per-donor mean ~112), consistent with metastatic SV burden being roughly 2.5x that of primary. Independent single-genome truth sets bracket the same range: **COLO829 = 68 somatic SVs**, **HCC1395 ≈ 1,300–1,788** — both genuine single-genome counts, useful if you later want to reproduce a specific benchmark rather than a generic median.

### Recommended target: ~200 SVs (a typical SV-active metastatic solid tumor)

Round the cohort median 193 -> **200** for clean per-type splits. This is the defensible "typical" point: neither an SV-quiet type (~56–71) nor an SV-driven catastrophe type (~372–412). If you instead want to stress the caller with an SV-high profile, **scale every count below by ~2x to land at ~400** (ovarian/oesophageal-like).

| Burden tier | Target total SVs | Matching profile |
|---|---|---|
| Low | ~60 | NET / kidney / CLL-like |
| **Typical (recommended)** | **~200** | **pan-cancer metastatic median** |
| High | ~400 | ovarian / oesophageal |

## 3. Type composition (per genome, ~200 SVs)

Deletions are the single most common somatic SV subtype at **33% of all SVs** in the Hartwig cohort, and are "the most prevalent in every cancer type except stomach and oesophageal tumours" ([Priestley et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6872491/)). Beyond that, there is **no universal type ordering**: the PCAWG SV companion ([Li et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32025012/)) finds SV classes "assort unevenly across tumour types and patients," describing the spectrum as sixteen signatures of structural variation rather than a fixed rank order. The DEL > DUP > INV > INS ordering used below is therefore a **modeling choice for a representative SV-active solid tumor** (DEL clearly dominant per Hartwig; tandem-DUP the conventional #2 via the Tandem Duplicator Phenotype; INV and INS in the residual) — not a claimed universal law, and it does not hold for translocation-enriched types like stomach/oesophageal.

**Translocations (BND/TRA) are deliberately set to 0.** The svp caller treats the genome as a flat concatenated buffer with no inter-chromosomal edges, so TRA/BND events are out of scope (the existing config already omits them). In the real spectrum translocations are large — ~24% of all PCAWG somatic SVs — so the remaining DEL/DUP/INV/INS are re-normalized to 100% of the simulated set. Re-normalizing Hartwig's DEL≈33%-of-full-spectrum to the translocation-free set gives DEL ≈ 33 / (100 − 24) ≈ 43%; I round to:

| Type | % of simulated SVs | Count (target 200) | Notes |
|---|---|---|---|
| DEL (deletion) | 48% | 96 | most common subtype (33% of full spectrum, ~43% translocation-free); small DEL <10 kb most enriched in metastases |
| DUP (tandem duplication) | 26% | 52 | conventional #2 (Tandem Duplicator Phenotype) |
| INV (inversion) | 13% | 26 | residual; reciprocal inversions uncommon |
| INS (insertion: mobile-element + novel) | 13% | 26 | residual; MEI-dominated |
| BND / translocation | 0% | 0 | **out of scope — svp has no inter-chromosomal edges** |
| **Total** | 100% | **200** | |

DEL is nudged slightly above the raw ~43% because it is the most prevalent class in nearly every cancer type and small deletions are the metastasis-enriched class. INV/INS each take the residual ~13%; INS spans mobile-element insertions (the dominant biological INS class) and novel/templated sequence.

## 4. Size distributions (per type)

Both DEL and tandem-DUP have **multimodal** size distributions ([Li et al. 2020, Nature](https://www.nature.com/articles/s41586-019-1913-9)). insilicoSV draws sizes **uniformly within each `length_ranges` bin**, so a multimodal shape is approximated by **stacking several bins with counts weighted toward the real modes** — standard practice for SV simulators (VISOR / insilicoSV / BAMSurgeon all bin-and-weight). Concrete bp modes come from the 560-breast-cancer rearrangement signatures ([Nik-Zainal et al. 2016, PMC4910866](https://pmc.ncbi.nlm.nih.gov/articles/PMC4910866/)), the Tandem Duplicator Phenotype span classes ([Menghi et al. 2018, PMC6481635](https://pmc.ncbi.nlm.nih.gov/articles/PMC6481635/)), and gnomAD-SV mobile-element peaks ([Collins/Brand et al. 2019](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674.full.pdf)).

**DEL — 96 events.** Min 50 bp, median ≈ 1–2 kb, max ~1 Mb. Bimodal: a dominant <10 kb mode plus a long tail; deletion rearrangement signature RS5 peaks <100 kb.

| Size range | Count | Rationale |
|---|---|---|
| 50 bp – 1 kb | 38 | small-DEL mode, the metastasis-enriched class |
| 1 kb – 10 kb | 34 | second small mode (PCAWG splits DEL at 10 kb) |
| 10 kb – 100 kb | 18 | Nik-Zainal deletion signature RS5 (<100 kb) |
| 100 kb – 1 Mb | 6 | long tail |

**DUP (tandem) — 52 events.** Min 50 bp, median ≈ 10–15 kb, max ~2 Mb. Three TDP span classes — class 1 median **11 kb** (1.64–51 kb), class 2 median **231 kb** (51–622 kb), class 3 median **1.7 Mb** (622 kb–6.2 Mb) ([Menghi 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6481635/)); signatures peak <10 kb (RS3, BRCA1-type) and >100 kb (RS1).

| Size range | Count | Rationale |
|---|---|---|
| 50 bp – 1 kb | 12 | short-TD mode |
| 1 kb – 51 kb | 22 | TDP class 1, modal ~11 kb (BRCA1-type small TD) |
| 51 kb – 622 kb | 13 | TDP class 2, modal ~231 kb |
| 622 kb – 2 Mb | 5 | TDP class 3, modal ~1.7 Mb (CDK12-type large TD) |

**INV — 26 events.** Min 50 bp, median ≈ 1–5 kb, max ~1 Mb. No strong canonical modes; spread broadly, weighted small.

| Size range | Count |
|---|---|
| 50 bp – 1 kb | 9 |
| 1 kb – 10 kb | 10 |
| 10 kb – 100 kb | 5 |
| 100 kb – 1 Mb | 2 |

**INS — 26 events.** Min 50 bp, median ≈ 300 bp, max ~10 kb. Dominated by mobile-element peaks: **Alu ~300 bp, SVA ~2.1 kb, full-length LINE-1 ~6 kb** ([gnomAD-SV](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674.full.pdf)), plus short novel/templated insertions.

| Size range | Count | Rationale |
|---|---|---|
| 50 – 350 bp | 13 | Alu peak ~300 bp (most common MEI) + short templated |
| 1.5 kb – 2.5 kb | 4 | SVA peak ~2.1 kb |
| 5 kb – 7 kb | 5 | full-length LINE-1 peak ~6 kb |
| 1 kb – 10 kb | 4 | residual novel-sequence insertions |

## 5. insilicoSV config (drop-in replacement for the `>= 50 bp` block of `simulation/configs/wgs.yaml`)

Keep the existing header (`reference`, `filter_small_chr: 10000000`, `min_intersv_dist: 200`, `haploid: True`, `max_tries`, `random_seed`) and the existing SNP/small-indel layer (the ~16,950 SNV and 470+470 indel counts are separate per-genome figures, not affected here). Replace only the large-SV `variant_sets` entries:

```yaml
    # ===== somatic SVs (>= 50 bp), ~200 total =====
    # Per-genome MEDIAN, representative SV-active metastatic solid tumor.
    # Anchor: Hartwig metastatic median 193 SVs/tumour (Priestley 2019,
    #   PMC6872491); rounded to 200. Per-type medians span ~56-412.
    # Composition: DEL most common (33% full / ~43% translocation-free,
    #   Hartwig); DEL>DUP>INV>INS is a modeling choice, not a universal
    #   law (Li 2020 finds SV classes assort unevenly across tumor types).
    # Sizes: Nik-Zainal 2016 RS1/RS3/RS5; TDP span classes (Menghi 2018);
    #   gnomAD-SV MEI peaks Alu~300bp/SVA~2.1kb/L1~6kb. DEL+DUP multimodal.
    # Translocations omitted: svp treats the genome as a flat buffer
    #   (no inter-chromosomal edges).

    # --- DEL = 96 (bimodal <10kb + long tail) ---
    - {type: "DEL", number: 38, length_ranges: [[50, 1000]]}
    - {type: "DEL", number: 34, length_ranges: [[1000, 10000]]}
    - {type: "DEL", number: 18, length_ranges: [[10000, 100000]]}
    - {type: "DEL", number: 6,  length_ranges: [[100000, 1000000]]}

    # --- DUP (tandem) = 52 (3 TDP span classes ~11kb/~231kb/~1.7Mb) ---
    - {type: "DUP", number: 12, length_ranges: [[50, 1000]]}
    - {type: "DUP", number: 22, length_ranges: [[1000, 51000]]}
    - {type: "DUP", number: 13, length_ranges: [[51000, 622000]]}
    - {type: "DUP", number: 5,  length_ranges: [[622000, 2000000]]}

    # --- INV = 26 ---
    - {type: "INV", number: 9,  length_ranges: [[50, 1000]]}
    - {type: "INV", number: 10, length_ranges: [[1000, 10000]]}
    - {type: "INV", number: 5,  length_ranges: [[10000, 100000]]}
    - {type: "INV", number: 2,  length_ranges: [[100000, 1000000]]}

    # --- INS = 26 (MEI peaks: Alu ~300bp, SVA ~2.1kb, L1 ~6kb + novel) ---
    - {type: "INS", number: 13, length_ranges: [[50, 350]]}
    - {type: "INS", number: 4,  length_ranges: [[1500, 2500]]}
    - {type: "INS", number: 5,  length_ranges: [[5000, 7000]]}
    - {type: "INS", number: 4,  length_ranges: [[1000, 10000]]}
```

Config file path: `/home/yuio/svp/simulation/configs/wgs.yaml`.

## 6. Caveats

- **The 50 bp minimum is an arbitrary, technology-driven convention — not a biological SV/indel boundary.** gnomAD-SV defines SVs as rearrangements ≥50 bp, and Manta (default `minScoredVariantSize=50`) and Severus (`--min-sv-size 50`) default to it ([gnomAD-SV 2019](https://gnomad.broadinstitute.org/news/2019-03-structural-variants/); [Manta config](https://github.com/Illumina/manta/blob/master/src/python/bin/configManta.py.ini); [Severus](https://github.com/KolmogorovLab/Severus)). But the standard SV-methods review states the 50 bp cutoff "was selected for historical and technical reasons related to the limitations of short-read sequencing, not because a 50-bp threshold is biologically significant in any way" ([Mahmoud et al. 2019, Genome Biology](https://link.springer.com/article/10.1186/s13059-019-1828-7)). Long-read callers — the relevant class for svp's HiFi inputs — routinely place 30–49 bp events in the SV layer **by default**: Sniffles2 `--minsvlen 35`, cuteSV minimum 30, pbsv 20. So if you later push svp toward long-read-native benchmarking, consider lowering the SV-layer minimum from 50 to ~30–35 bp; with a 50 bp floor here, 1–49 bp events stay in the small-indel layer purely by short-read convention.
- **insilicoSV samples sizes uniformly within each bin**, so these stacked bins *approximate* the true multimodal/log-normal shapes by weighting counts toward the real modes — they are not draws from the empirical length distribution. Adequate for caller benchmarking; if exact mode fidelity matters later, narrow each bin around its mode.
- **This is a non-catastrophic genome:** no chromothripsis cluster, no whole-genome doubling, no complex templated-insertion chains. About 29% of real tumors have chromothripsis (tens-to-hundreds of clustered breakpoints in one event) and a few hypermutated genomes carry thousands of SVs — those live in the right tail and need a separate, explicitly "complex/SV-high" config (scale to ~400+ and add a clustered DEL/INV burst).
- **To match a specific benchmark instead of a generic median,** use a real single-genome truth set directly: **COLO829** (68 somatic SVs: 38 DEL, 7 DUP, 7 INV, 3 INS, 13 TRA) or **HCC1395** (~1,788 SVs). Both are single-genome counts, not cohort aggregates.
