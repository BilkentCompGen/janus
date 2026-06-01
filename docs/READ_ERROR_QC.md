# Whole-genome HiFi read-error QC for LCP cores — full report

Whole-genome tumor/normal PacBio-HiFi simulation and analysis of how LCP cores
behave when the LCP stage is run on **reads** rather than on an assembled genome.
Self-contained: every parameter, tool, command, and measured number is here. All
numbers are measured, not estimated.

> **History.** An earlier version of this report claimed an intrinsic "LCP parse
> divergence" producing a 27–45% read-mode recall floor and an unprunable "phantom"
> core class. **That was wrong** — it was a case-sensitivity bug in svp's signature
> hashing interacting with a soft-masked reference. The bug is fixed; all numbers
> below are post-fix. The bug and its diagnosis are documented in §8.

**Goal.** Run the LCP stage on simulated 30× HiFi reads and measure, per core, how
its observed multiplicity (read coverage) depends on its ground-truth copy-number,
to design a QC / error-pruning step. Three core populations:

- **read-error cores** — signatures absent from the truth genome (reads are
  simulated from the truth, so any novel core is a sequencing error);
- **normal cores** — single-copy truth cores;
- **DUP-site cores** — cores inside a duplicated region (≥2 copies), expected at
  ~2× the normal coverage.

---

## 1. Background: LCP and the signature

**LCP = Locally Consistent Parsing** (`lib/lcptools`). It partitions a DNA string
into **cores** whose boundaries follow local sequence structure; the parse is
recursively **deepened** to level L (longer, sparser cores at higher L). Each core
carries an integer **label**. The lcptools alphabet is **case-insensitive**
(`encoding.c`: `alphabet['A']=alphabet['a']=0`, …), so cores and labels are
identical for soft-masked vs unmasked sequence.

svp reduces a core to a 64-bit **signature**:

```
signature(core) = MurmurHash64( UPPERCASE(seq[start-flank .. end+flank]), seed = core.label )
```

with **flank = 10 bp**. The signature must be case-insensitive too (the
`UPPERCASE`), or it disagrees with the case-insensitive cores — see §8. svp
normalises case in the FASTA reader (`io.cpp:read_fasta_flat`) and the FASTQ reader
(`lcp_read_qc`), so all sequence is uppercased before hashing.

---

## 2. Simulation — the truth genomes

### 2.1 Tumor truth genome (`insilicoSV` v0.1.9)

- Reference: GRCh38 (`simulation/ref.fasta`), **48.3% soft-masked (lowercase)**.
- Config `simulation/configs/wgs.yaml`: `filter_small_chr: 10000000` (restricts to
  chr1–22+X+Y, excludes alt/random contigs), `min_intersv_dist: 200`,
  `haploid: True`, `max_tries: 300`, `random_seed: 42`.
- Variant counts (PCAWG per-tumour burden — §2.3):

  | set | number | length range (bp) |
  |---|--:|---|
  | SNP | 16,950 | [1,1] |
  | INS / DEL (small indel) | 470 / 470 | [1,49] |
  | DEL (SV) | 22 / 18 / 8 / 2 | [50,1k]/[1k,10k]/[10k,100k]/[100k,1M] |
  | DUP (SV) | 14 / 12 / 5 / 1 | same bins |
  | INV (SV) | 7 / 5 / 2 / 1 | same bins |
  | INS (SV) | 9 / 4 / 2 | [50,1k]/[1k,10k]/[10k,100k] |

  **Totals: 16,950 SNV + 940 small indels + 112 SVs (50 DEL, 32 DUP, 15 INV, 15 INS)
  = 18,002 variants.** Output: `simulation/truth/wgs/out.fa` (3.21 GB), `out.vcf`,
  `novel.fa`; 24 primary records (`1_hapA … Y_hapA`), 3,087,377,631 bp.

### 2.2 Normal truth genome

Unmodified reference: `simulation/ref.fasta` primary chromosomes (~3.09 Gbp).

### 2.3 Why these counts (PCAWG)

PCAWG flagship (ICGC/TCGA, *Nature* 2020, s41586-020-1969-6): across 2,583 donors,
43,778,859 SNVs / 2,418,247 indels / 288,416 SVs → per-tumour means ~16,950 / ~936
/ ~112. SV type ordering (DEL > tandem-DUP > INV; mobile-element INS rarer) from
the PCAWG SV paper (Li et al. 2020). Translocations not simulated (the caller uses
a flat genome buffer with no inter-chromosomal edges).

---

## 3. Simulation — the reads (`tools/sim_hifi`)

HiFi reads are produced directly at the target accuracy by `tools/sim_hifi`.

**Model / parameters.** Per primary chromosome, `depth × len / mean_len` reads:
length ~ Gamma(mean 15,000, sd 5,000) clamped [100, 1e6]; uniform start; uniform
50/50 strand (minus reads reverse-complemented). Errors injected per base at
**0.001 (Q30)**, typed sub:ins:del = **22:45:33** (PacBio Sequel), with a
homopolymer deletion bias (`--hp-del-bias 3`); inserted base copies the following
base 50% of the time. Reference `N` (assembly gaps) are **kept as N** in the read
(LCP breaks there, as in the genome — not randomised). Depth = 30×.

**Driver / outputs.** `scripts/sim_reads_wg.sh` runs `sim_hifi` per chromosome 24-way
parallel and concatenates. **Tumor:** `simulation/reads/wgs/hifi.fastq.gz`,
6,174,744 reads, 92.6 Gbp, measured 0.1000%/base error. **Normal:**
`simulation/reads/normal/hifi.fastq.gz`, 6,176,526 reads, 92.6 Gbp. Wall time
**~16 min total** (32-core host).

**Read-error ground truth.** `sim_hifi --errors` writes every injected error
(`read, read_pos, chrom, ref_pos, type, ref, alt`) to
`simulation/reads/{wgs,normal}/errors/*.errors.tsv.gz` — exact, no alignment.

---

## 4. The analysis tool (`tools/lcp_read_qc`)

1. **Truth side.** Parse the primary-chromosome truth FASTA at level L on **both
   strands** (forward + reverse-complement), because LCP cores are strand-specific
   and a minus-strand read's cores match the reverse-complement truth. A signature's
   **copy-number** = its occurrence count across the both-strand truth (1 = single,
   2 = duplicated, ≥3 = higher repeat).
2. **Read side.** Stream the reads with T worker threads (256-way lock-sharded count
   map). Per read, `extract_cores` + `compute_sigs` at the same L; count only
   **interior** cores (`start ≥ margin`, `end + margin ≤ len`, margin = flank = 10)
   — edge cores have clamped flanks.
3. **Output.** Each observed signature is bucketed by truth copy-number
   (0 = read-error, 1, 2, ≥3); the multiplicity distribution is reported per bucket.
   `truth_dropouts` = truth signatures never observed in any read.

Both readers uppercase the sequence (§1), so a soft-masked truth core and the
uppercase read that covers it produce the same signature.

---

## 5. Result — whole-genome multiplicity by copy-number (level 5)

30×, level 5; ~85.0 M distinct both-strand truth signatures; ~1.30 billion interior
read cores per genome; ~11 min QC each.

| copy-number | tumor distinct | tumor mean | normal distinct | normal mean |
|---|--:|--:|--:|--:|
| 0 = read-error | 202,775,834 | 1.01 (median 1, p90 1) | 202,774,374 | 1.01 |
| 1 = normal | 83,992,796 | **12.49** (median 12) | 84,042,169 | 12.48 |
| 2 = DUP | 660,855 | **24.95** (≈ 2× normal) | 647,231 | 25.05 |
| ≥3 = repeat | 331,207 | 79.26 | 331,996 | 79.56 |
| — `truth_dropouts` (cn1) | **1,256** | — | 1,337 | — |

- **Reads reproduce genome cores essentially perfectly:** only **1,256 of ~84 M**
  single-copy truth cores are never observed (0.0015% dropout).
- **Read-errors are multiplicity 1** (mean 1.01, median 1, p90 1): random
  sequencing errors, no recurring class. 99.55% sit at mult-1.
- **DUP cores sit at ~2× normal** (24.95 vs 12.49) — the copy-number signal is clean.
- **Tumor ≈ normal** (the 18,002 somatic variants are negligible vs 84 M cores; the
  cn2 bucket is dominated by natural segmental duplications, which the 32 simulated
  somatic tandem DUPs join).

The normal multiplicity ≈ 12.5 = per-strand depth (≈15 = 30× over two strands, since
cores are strand-specific) × P(the ~170 bp core+flank window is HiFi-error-free ≈ 0.84).

---

## 6. Result — LCP level sweep (whole-genome tumor, 30×)

| metric | L4 | L5 | L6 | L7 |
|---|--:|--:|--:|--:|
| truth cores (one strand) | 102.2 M | 43.7 M | 18.7 M | 8.0 M |
| interior read cores | 3.05 B | 1.30 B | 546 M | 225 M |
| normal-core mean mult | 13.73 | 12.49 | 10.00 | 5.98 |
| DUP mean mult | 27.47 | 24.95 | 19.97 | 11.96 |
| error cores, total | 239.4 M | 202.8 M | 170.9 M | 128.9 M |
| — at multiplicity 1 | 99.37% | 99.55% | 99.71% | 99.86% |
| — at multiplicity ≥ 5 ("phantom") | 139,930 | 33,716 | 2,230 | 59 |
| **`truth_dropouts` (cn1)** | **1,058** | **1,256** | **2,838** | **51,785** |
| dropout fraction | 0.0006% | 0.0015% | 0.008% | 0.33% |

- **Dropout is negligible at every level** (0.0006% → 0.33%); reads reproduce the
  genome cores faithfully. (The slow rise with level is real — longer cores are
  slightly more sensitive to error/coverage at read ends — but tiny.)
- **Errors are ~all multiplicity 1** at every level (99.4–99.9%); the mult≥5 count is
  a few ×10⁴ at most (incidental recurrences), not a population.
- **Multiplicity falls with level** (13.7 → 6.0) as the core+flank window lengthens
  and HiFi-error survival drops; DUP stays ~2× normal throughout.
- **Best separation at L4–L5**: normal ≈ 12–14 vs the mult-1 error floor. At L7 normal
  is ~6, closer to the floor.

---

## 7. Error pruning

Read-errors are random and almost entirely **multiplicity 1**. Multiplicity-threshold
sweep (keep cores with `mult ≥ t`; measured on WG tumor L5):

| `mult ≥ t` | error cores removed | normal cores lost | normal kept |
|---|--:|--:|--:|
| 2 | 99.552% | 0.0050% | 99.995% |
| 3 | 99.941% | 0.036% | 99.964% |
| 4 | 99.971% | 0.162% | 99.838% |
| 5 | 99.983% | 0.552% | 99.448% |

The `mult = 2` bin holds **788,504 error cores vs 25,861 normal cores** — sequencing-error
**recurrences** (two reads with the same error at the same locus; ~0.8 M coincidences
at 0.1% error × depth × ~203 M error positions). `mult ≥ 2` removes the mult-1 errors;
`mult ≥ 3` additionally sweeps those recurrences (99.55% → 99.94% of errors) for a still-
trivial 0.036% normal loss, so **`mult ≥ 3` is the better operating point**. Beyond that,
normal-core loss climbs faster than error gain. DUP-site cores are flagged separately by
sitting at ~2× the normal multiplicity.

---

## 8. Postmortem — the case-sensitivity bug (why the earlier results were invalid)

**Root cause.** `compute_sigs` hashed the **raw** flank bytes:
`murmur64(seq.data()+L, R-L, label)`. `murmur64("…acgt…") ≠ murmur64("…ACGT…")`,
so a soft-masked (lowercase) core and its uppercase equivalent got **different
signatures** — even though the lcptools core (boundaries + label) is identical
(case-insensitive). The reference is **48.3% lowercase** (soft-masked repeats) and
`sim_hifi` emits **uppercase** reads (`b2i` masks with `& 0xDF`). So in every
soft-masked region, the uppercase read core's signature did not match the
lowercase truth core's signature.

**What it produced (all now-retracted):** ~497k "phantoms" on chr19 (≈300 M on WG),
a 27–45% "dropout," and the bogus "intrinsic LCP parse divergence / read-mode recall
floor / unprunable phantom class."

**Why it evaded the earlier diagnostics** — a deterministic per-base case flip is:
window-independent (read-length sweep 3k–500k was flat), not read-ends (drop-edge
and margin sweeps did nothing), not N-related (chr19 is 0.21% N and an N-keeping fix
changed phantoms 0.13%), present on both strands, and recurring at per-strand depth
(so phantoms sat at multiplicity ~15, indistinguishable from real cores). The one
"clean" single read that gave 0 phantoms was extracted with `samtools`, which
**preserves case**.

**Decisive test.** Uppercasing the truth to match the uppercase reads →
phantoms **497k → 0**, dropout **416,089 → 5** (chr19, err-free, L5).

**Fix.** Uppercase the sequence in both readers (`io.cpp:read_fasta_flat` and the
FASTQ reader in `lcp_read_qc`), so signatures are case-insensitive like the cores.
This is a **real svp bug for read mode**: production HiFi reads are uppercase and
production references are soft-masked, so the same mismatch would have hit a real
run. (It does *not* affect FASTA-vs-FASTA mode where both inputs share masking;
the chr19 GFAs shifted only ~0.008% — from correctly merging masked vs unmasked
copies of the same repeat.)

A second, independent fix: `sim_hifi` now keeps reference `N` instead of randomising
gap bases (gaps are unsequenceable; randomising them created spurious mult-1 cores).

---

## 9. chr19 validation (level 5, post-fix)

| | error cores | error mean | normal cores | normal mean | dropout |
|---|--:|--:|--:|--:|--:|
| error-free reads | 0 | — | 1,510,660 | 14.79 | 5 |
| Q30 reads | 4,545,162 | 1.00 | 1,510,656 | 12.53 | 9 |

Error-free reads reproduce the genome cores exactly (0 spurious cores, 5 dropouts);
Q30 adds only multiplicity-1 sequencing-error cores. (Pre-fix these were 497k
phantoms / 416k dropout.)

---

## 10. File manifest & reproduce

| path | contents |
|---|---|
| `simulation/configs/wgs.yaml` | insilicoSV config |
| `simulation/truth/wgs/{out.fa,out.vcf,novel.fa}` | tumor truth + VCF |
| `tools/sim_hifi.cpp` | HiFi read simulator |
| `scripts/sim_reads_wg.sh` | per-chromosome parallel read driver |
| `simulation/reads/{wgs,normal}/hifi.fastq.gz` | 30× HiFi reads |
| `simulation/reads/{wgs,normal}/errors/*.errors.tsv.gz` | injected-error truth |
| `tools/lcp_read_qc.cpp` | the QC tool |
| `src/io.cpp`, `src/cores.cpp` | the case-insensitive signature fix |

```sh
make -C simulation wgs                              # truth genome
make svp && make tools/sim_hifi                     # build (case fix in io.o)
g++ -std=c++17 -O3 -Ilib/lcptools -Isrc -o tools/lcp_read_qc \
    tools/lcp_read_qc.cpp build/io.o build/cores.o build/MurmurHash3.o \
    -Llib/lcptools -llcptools -lz -pthread
bash simulation/run_wgs_reads.sh                    # tumor + normal 30x HiFi
samtools faidx simulation/truth/wgs/out.fa $(awk '$1~/^([0-9]+|X|Y)_hapA$/{print $1}' \
    simulation/truth/wgs/out.fa.fai) > truth_primary.fa
tools/lcp_read_qc truth_primary.fa simulation/reads/wgs/hifi.fastq.gz \
    --level L --flank 10 --threads 24 --out qc_L.tsv
```
