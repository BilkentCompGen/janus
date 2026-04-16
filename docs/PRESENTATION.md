# SVP: Alignment-Free Somatic Structural Variant Detection via Locally Consistent Parsing

**Progress Presentation**

---

## Slide 1 — What Are Structural Variants?

Structural variants (SVs) are genomic rearrangements ≥50bp that alter chromosome structure:

| Type | What happens | Example |
|------|-------------|---------|
| **Deletion (DEL)** | Segment removed | `ABC` → `AC` |
| **Insertion (INS)** | Novel sequence added | `AC` → `AXC` |
| **Duplication (DUP)** | Segment copied in tandem | `ABC` → `ABBC` |
| **Inversion (INV)** | Segment reversed in place | `ABC` → `A(rc B)C` |

SVs account for more divergent bases between any two human genomes than SNPs despite being far fewer in number (~4,400 germline SVs vs. ~4M SNPs per individual).

*Abel et al., Nature 2020; Collins et al., Nature 2020*

---

## Slide 2 — SVs in Cancer

Somatic SVs are a hallmark of cancer genomes. They drive tumorigenesis by disrupting tumor suppressors, creating oncogenic fusions, and amplifying driver genes.

**Scale (PCAWG, 2,658 tumors across 38 cancer types):**
- Median ~100-300 somatic SVs per tumor (varies widely by cancer type)
- Breast/ovarian cancers: 5-6 driver SVs per tumor on average
- Deletions most common, followed by tandem duplications, then inversions
- Over half of all SV breakpoints occur within complex rearrangement clusters

**Additionally:** ~10,000 somatic SNVs (~3/Mb for moderate TMB) and ~2,000 somatic short indels per tumor genome.

*PCAWG Consortium, Nature 2020; Li et al., Nature 2020*

---

## Slide 3 — How Other SV Callers Work

All major SV callers require **aligned reads (BAM files)** as input:

| Caller | Approach | Key Signals |
|--------|----------|-------------|
| **Delly** | Split-read + discordant pairs | Abnormal insert sizes, clipped reads |
| **Manta** | Discordant pairs + local assembly | Graph-based candidate generation |
| **GRIDSS** | Positional de Bruijn graph assembly | Soft-clipped, discordant, unmapped reads |
| **SvABA** | String graph local assembly | Gapped/clipped/discordant reads |
| **LUMPY** | Probabilistic integration | Multiple alignment signals |

**Common pipeline:** FASTQ → Aligner (BWA/minimap2) → BAM → SV Caller → VCF

**Limitations:** Alignment artifacts, reference bias, missed complex rearrangements, high computational cost of alignment.

*Rausch et al., Bioinformatics 2012; Chen et al., Bioinformatics 2016; Cameron et al., Genome Research 2017; Wala et al., Genome Research 2018; Layer et al., Genome Biology 2014*

---

## Slide 4 — SVP: A Different Approach

SVP skips alignment entirely. It uses **Locally Consistent Parsing (LCP)** to decompose genomes into deterministic subsequences ("cores"), builds a colored graph, and finds structural differences as divergent paths.

**Current stage:** Works on FASTA genome assemblies.
**Pipeline:** Reference FASTA + Normal FASTA + Tumor FASTA → **SVP** → VCF

**Roadmap:** Extend to work directly on **raw long reads (FASTQ)** — PacBio / ONT — calling SVs without any alignment step at all. The LCP framework is designed for this: parsing is alignment-free by construction.

**Advantages:**
- No alignment step — no alignment artifacts or reference bias
- Detects all SV types (DEL, INS, DUP, INV) + SNPs in a single pass
- Linear-time parsing: O(N) per genome
- Natural path to direct long-read SV calling without BAM intermediates

---

## Slide 5 — Algorithm Step 1: LCP Parsing

Each input genome is parsed into **LCP cores** — deterministic, locally defined subsequences.

**How:** Scan the sequence left-to-right. At each position, detect local extrema patterns in the alphabet ordering:
- **Local minimum** (LMIN): `a > b < c` — a valley
- **Local maximum** (LMAX): `a < b > c` — a peak
- **Run-interrupted** (RINT): flanking a run of identical characters
- **Spanning sequences** (SSEQ): fill gaps between cores

**Result:** ~N/1.5 cores per genome after level-1 parsing.

**Deepening:** Iteratively compress via Deterministic Coin Tossing (DCT) + re-parsing. Each level reduces core count by ~1.5×. After d levels: ~N/1.5^d cores.

**Complexity:** O(N) time and space, run once per genome (×3: ref, normal, tumor).

---

## Slide 6 — Algorithm Step 2: Colored Graph Construction

Cores are hashed and inserted into a graph. Each node represents a unique core sequence.

- **Normal genome cores** get color bit `0x1`
- **Tumor genome cores** get color bit `0x2`
- Cores appearing in both get **shared color** `0x3`

Consecutive cores are linked with directed edges. The graph encodes the structure of all genomes simultaneously.

**Shared regions:** Most of the genome is identical between normal and tumor → shared nodes.
**Divergent regions:** Structural variants create runs of single-color nodes between shared anchors.

---

## Slide 7 — Algorithm Step 3: Bubble Detection & Extraction

A **bubble** is a divergent subgraph between two shared anchor nodes:

```
         ── normal-only path ──
shared → ── ── ── ── ── ── ── → shared
         ── tumor-only path  ──
```

**Detection:** Scan nodes for shared→single-color transitions (bubble start) and single-color→shared transitions (bubble end). O(V + E).

**Extraction:** Breadth-first search from each bubble start along normal-only and tumor-only edges. Sequences are reconstructed using **KMP-based overlap merging** (handles overlapping LCP core boundaries).

Each bubble yields a (normal sequence, tumor sequence) pair — the raw material for SV calling.

---

## Slide 8 — Algorithm Step 4: Reference Anchoring

Bubble positions must be mapped back to reference genome coordinates.

**Method:** Each core has a hash. A reference index maps core hashes to their reference positions. If a bubble's start node has a unique reference position → direct anchor (O(1)).

**Fallback:** Walk backward through predecessor nodes in the graph until a uniquely-anchored core is found, then compute the offset. Forward walk via the end node as a last resort.

**Typical case:** O(1) — most bubble starts anchor directly.

---

## Slide 9 — Algorithm Step 5: SV Calling

For each anchored bubble, classify the variant:

1. **Quick scan:** Compare normal and tumor sequences from both ends inward. Identifies SNPs, clean indels, and flags complex cases. O(L).

2. **Complex variants:** SIMD-accelerated **Needleman-Wunsch** global alignment + CIGAR parsing. Groups operations into variant blocks separated by anchor matches. O(L²).

3. **Inversion detection** (two pipelines):
   - **Bubble-based:** Smith-Waterman local alignment of tumor sequence vs. reverse complement of normal sequence. Normalized score ≥0.90 → INV.
   - **Core-index:** For large inversions. Classifies tumor cores against forward-ref and reverse-complement-ref hash sets. Runs of reverse-complement-only cores → INV call.

4. **Duplication detection** (two mechanisms):
   - **Flipped DUP:** Color pattern matching on bubble boundary nodes — detects the distinctive signature of tandem duplications in the graph.
   - **INS→DUP reclassification:** Smith-Waterman of inserted sequence vs. adjacent reference. Score ≥0.70 → reclassify as duplication.

---

## Slide 10 — Overall Complexity

| Phase | Time | Space |
|-------|------|-------|
| LCP parsing (×3 genomes) | O(N) | O(N) |
| Graph construction | O(C_d) | O(V + E) |
| Reference indexing | O(C_d) | O(C_d) |
| Bubble detection | O(V + E) | O(B) |
| Path extraction | O(Σ path lengths) | O(L_max) |
| SV calling (simple) | O(L) per bubble | O(1) |
| SV calling (alignment) | O(L²) per bubble | O(L²) |

**Overall: O(N + B·L²)** where N = genome length, B = number of bubbles, L = avg bubble sequence length.

**Dominant cost:** Needleman-Wunsch alignment for complex variants (mitigated by SIMD vectorization).

---

## Slide 11 — Simulation Setup

**Somatic simulation** (GRCh38 reference → tumor genome):
- 10,000 somatic SNPs (~3/Mb, moderate tumor mutational burden)
- 2,000 somatic short indels (1-49bp)
- **300 somatic SVs:**
  - 120 DEL (50bp-1Mb)
  - 70 DUP (50bp-1Mb)
  - 60 INS (50bp-100kb)
  - 50 INV (50bp-1Mb)

**Genome-wide** across all 22 autosomes (~3.1 Gb). Reference genome serves as the normal. Generated with insilicoSV.

*SV type ratios and size distribution based on PCAWG (Li et al., Nature 2020). Somatic SNV rate based on moderate-TMB solid tumors (~3/Mb).*

---

## Slide 12 — Published Accuracy of Other Callers (for context)

All require aligned BAM input. Numbers from published benchmarks on simulated tumor-normal WGS at 60× coverage.

**Overall sensitivity and precision (all SV types combined):**

| Caller | Sensitivity | Precision | Notes |
|--------|------------|-----------|-------|
| **Manta** | ~91% | >90% | Best standalone performer; best for INS (79%) |
| **LUMPY** | ~77% | >90% | Good DEL recall; highest call count |
| **GRIDSS** | ~75% | >90% | Best breakpoint precision (±2bp); high memory (38 GB) |
| **Delly** | ~70% | >90% | Good INV recall (>80%); slow runtime |
| **SvABA** | ~63% | >90% | Best for large SVs (>1Mb); assembly-based |

**By SV type (key observations):**
- **DEL:** All callers perform best here. Manta F1=74%, Precision=94%, Recall=61% on GIAB HG002.
- **INV:** Manta, GRIDSS, Delly recover >80%; others struggle.
- **INS:** Only Manta achieves meaningful detection (79%); others <30%.
- **DUP:** All callers F1 < 0.2 — tandem duplications remain challenging.

**Runtime and memory (30× WGS, single sample):**

| Caller | Wall time | Peak memory |
|--------|-----------|-------------|
| **Manta** | ~49 min | 2.8 GB |
| **Delly** | ~54 min | 0.9 GB |
| **LUMPY** | ~29 min | 5.7 GB |
| **GRIDSS** | ~8.4 hr | 38.5 GB |

*Sources: Kosugi et al., Genome Biology 2019; Cameron et al., BMC Genomics 2024; Smolka et al., PMC 2025; Kuzniar et al., GigaScience 2020*

Note: Direct comparison with SVP is not apples-to-apples — these callers work from aligned reads (BAM), while SVP works from assembled genomes (FASTA). SVP's runtime includes genome parsing but not the upstream alignment step that BAM-based callers require.

---

## Slide 13 — SVP Runtime Performance

AMD Ryzen 9 9950X, single-threaded. 22 autosomes (~3.1 Gb), LCP depth=5.

| Metric | Value |
|--------|-------|
| **Total CPU time** | 1,347 s (~22 min) |
| **Wall time (8 cores)** | ~3 min |
| **Peak memory (per chr)** | 11.4 GB (chr2) |

Runtime and memory scale linearly with chromosome size. Per-chromosome processing is trivially parallelizable. No alignment step needed.

---

## Slide 14 — SVP Accuracy by SV Type (Truvari, SVs ≥50bp)

300 somatic SVs + 12,000 somatic SNPs/indels across 22 autosomes. Evaluated with Truvari.

| Type | Truth | TP | FP | FN | Precision | Recall | F1 |
|------|-------|----|----|----|-----------|--------|----|
| DEL | 120 | 85 | 137 | 35 | 38.3% | 70.8% | 49.7% |
| INS | 60 | 44 | 178 | 16 | 19.8% | 73.3% | 31.2% |
| DUP | 70 | 50 | 197 | 20 | 20.2% | 71.4% | 31.5% |
| INV | 50 | 44 | 8 | 6 | 84.6% | 88.0% | 86.3% |
| **ALL** | **300** | **250** | **539** | **50** | **31.7%** | **83.3%** | **45.9%** |

**Recall is strong** — SVP finds 83% of true somatic SVs across all types and sizes.

**INV detection is the standout** — 88% recall, 85% precision (F1=86%).


---

## Slide 15 — Accuracy by SV Size

Precision improves dramatically with SV size — the FP problem is concentrated in small SVs (50-500bp).

| Min size | Truth | TP | FP | FN | Precision | Recall | F1 |
|----------|-------|----|----|-----|-----------|--------|----|
| ≥50bp | 300 | 250 | 539 | 50 | 31.7% | 83.3% | 45.9% |
| ≥100bp | 292 | 242 | 260 | 50 | 48.2% | 82.9% | 61.0% |
| ≥500bp | 226 | 196 | 48 | 30 | 80.3% | 86.7% | 83.4% |
| ≥1kb | 221 | 191 | 44 | 30 | 81.3% | 86.4% | 83.8% |
| ≥5kb | 178 | 150 | 20 | 28 | 88.2% | 84.3% | 86.2% |
| ≥10kb | 115 | 91 | 18 | 24 | 83.5% | 79.1% | 81.3% |

At ≥500bp: **80% precision, 87% recall, F1=83%**. At ≥5kb: **88% precision, F1=86%**.

---

## Slide 16 — Accuracy by Chromosome (SVs ≥50bp)

| Chr | Truth | TP | FP | FN | Prec | Rec | F1 | | Chr | Truth | TP | FP | FN | Prec | Rec | F1 |
|-----|-------|----|----|-----|------|-----|-----|---|-----|-------|----|----|-----|------|-----|-----|
| 1 | 15 | 11 | 15 | 4 | 42.3% | 73.3% | 53.7% | | 12 | 12 | 10 | 20 | 2 | 33.3% | 83.3% | 47.6% |
| 2 | 14 | 11 | 12 | 3 | 47.8% | 78.6% | 59.5% | | 13 | 12 | 12 | 17 | 0 | 41.4% | 100% | 58.5% |
| 3 | 11 | 10 | 21 | 1 | 32.3% | 90.9% | 47.6% | | 14 | 13 | 12 | 12 | 1 | 50.0% | 92.3% | 64.9% |
| 4 | 13 | 11 | 22 | 2 | 33.3% | 84.6% | 47.8% | | 15 | 24 | 18 | 23 | 6 | 43.9% | 75.0% | 55.4% |
| 5 | 9 | 8 | 23 | 1 | 25.8% | 88.9% | 40.0% | | 16 | 12 | 11 | 36 | 1 | 23.4% | 91.7% | 37.3% |
| 6 | 19 | 17 | 19 | 2 | 47.2% | 89.5% | 61.8% | | 17 | 12 | 10 | 18 | 2 | 35.7% | 83.3% | 50.0% |
| 7 | 16 | 12 | 17 | 4 | 41.4% | 75.0% | 53.3% | | 18 | 13 | 9 | 20 | 4 | 31.0% | 69.2% | 42.9% |
| 8 | 20 | 16 | 25 | 4 | 39.0% | 80.0% | 52.5% | | 19 | 11 | 11 | 18 | 0 | 37.9% | 100% | 55.0% |
| 9 | 11 | 9 | 102 | 2 | 8.1% | 81.8% | 14.8% | | 20 | 8 | 7 | 26 | 1 | 21.2% | 87.5% | 34.1% |
| 10 | 13 | 12 | 17 | 1 | 41.4% | 92.3% | 57.1% | | 21 | 13 | 9 | 28 | 4 | 24.3% | 69.2% | 36.0% |
| 11 | 19 | 16 | 21 | 3 | 43.2% | 84.2% | 57.1% | | 22 | 10 | 8 | 27 | 2 | 22.9% | 80.0% | 35.6% |

Recall is consistently 69-100% across chromosomes. chr9 shows a FP spike from complex/repetitive regions.

---

## Slide 17 — Known Limitations

- **FP rate:** Somatic SNPs and short indels near SV breakpoints can produce spurious SV-sized calls. Filtering and post-processing can reduce this.
- **DUP precision:** Tandem duplications remain the hardest SV type — a known challenge across all callers (published benchmarks show F1 < 0.2 for most BAM-based tools as well).
- **Germline divergence:** When the normal genome diverges significantly from the reference, coordinate anchoring accuracy degrades. Active work on germline-aware anchoring.
- **FASTA-only input:** Currently requires assembled genomes; long-read FASTQ support is the next milestone.

---

## Slide 18 — Current Status & Next Steps

**Done:**
- Alignment-free SV calling from FASTA genomes using LCP
- Detects DEL, INS, DUP, INV, and SNPs in a single pass
- Linear-time genome parsing, quadratic only for complex variant alignment
- Genome-wide benchmarking: 83% recall, F1=83% at ≥500bp, INV F1=86%
- ~22 min single-core / ~3 min with 8 cores for full human genome

**Next steps:**
- Improve precision via call filtering and breakpoint refinement
- Extend LCP parsing to work directly on raw long reads (FASTQ)
- Call somatic SVs from PacBio / ONT reads without alignment
- Benchmark against Delly, Manta, GRIDSS on matched datasets
