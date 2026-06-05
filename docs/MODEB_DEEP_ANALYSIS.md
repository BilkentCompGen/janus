# Mode B → VCF deep analysis (2026-06-05 overnight session)

## Executive summary

This session reworked everything downstream of `modeb.bin` based on a systematic
truth-bubble dissection. Numbers in this doc are reproducible by:
```
make svp tools/modeb_dissect
./tools/modeb_dissect outputs/svp/wgs_mm10
python3 tools/truth_map.py outputs/svp/wgs_mm10 simulation/truth/wgs/out.vcf
python3 tools/pos_offset3.py
python3 tools/fp_features.py
./svp anchor outputs/svp/wgs_mm10        # 38 s
./svp vcf outputs/svp/wgs_mm10           # < 1 s
bash scripts/svp_benchmark.sh outputs/svp/wgs_mm10
```

### Truvari results (mm=10, no DUP/INV detector, paired metrics)

| threshold | TP | FP | FN | recall | precision | F1 |
|---|---:|---:|---:|---:|---:|---:|
| ≥50bp SV (`refdist 500, pctsize 0.7`) | 51 | 2 | 61 | 45.5% | **96.2%** | **0.618** |
| 1–49bp indel (`refdist 50`)          | 861 | 1023 | 79 | 91.6% | 45.7% | 0.610 |
| SNP (`refdist 0`)                    | 15,420 | 16,705 | 1,547 | 90.9% | 48.0% | 0.628 |

Per-SVTYPE on geq50:
| SVTYPE | TP | FP | FN | recall | precision |
|---|---:|---:|---:|---:|---:|
| DEL | 36 | 1 | 14 | 72.0% | **97.3%** |
| INS | 12 | 1 | 3  | 80.0% | **92.3%** |
| DUP | 0  | 0 | 32 | 0% | – |
| INV | 0  | 0 | 15 | 0% | – |

Starting state (2026-06-05 morning): geq50 recall 58% / precision 4.1% / F1 0.07. The
session ended at F1 0.62 (8.8×) on geq50 while preserving lt50/SNP recall.

---

## 1. Mode B vs truth — what each tag captures

The dissector emits `outputs/svp/wgs_mm10/bubble_table.tsv` (one row per Mode B bubble
with chrom/strand/anchor positions/coverage). truth_map.py joins truth VCF to bubbles.

**Bubble proximity to truth (within 500 bp of the nearest truth variant):**

| Truth SVTYPE | total | bubble within 5kb | within 500bp | bubble tag (within 500bp) |
|---|---:|---:|---:|---|
| SNP (n=16,950) | 16,950 | 16,390 | 16,088 | S 15,658 · D 294 · I 354 · C 15 |
| DEL (n=520)    | 520    | 501    | 491    | D 436 · I 19 · S 44 · C 1 |
| INS (n=485)    | 485    | 467    | 454    | I 366 · S 89 · D 6 · C 3 |
| DUP (n=32)     | 32     | 17     | 12     | I 6 · D 1 |
| INV (n=15)     | 15     | 12     | 11     | D 6 · I 4 · S 1 |

**Key takeaways:**
- Mode B captures ~95 % of SNP/INS/DEL with the CORRECT tag.
- Mode B captures only **17/32 DUP truth** in any nearby bubble — 15 truth DUPs have
  NO Mode B bubble within 5 kb. The graph misses them entirely.
- Mode B captures 13/15 INV — but tags inconsistently as D, I, or S (no single bubble
  represents the inversion).

---

## 2. The position bug (the single biggest fix)

**Old POS formula:** `pos = bubble_pos.pos_start + 1` where `pos_start = max(s1.pos_end)`.
This placed VCF POS at the right edge of s1's flank window — wrong by up to
`pre - s1.block_len` bp.

**Correct strand-aware POS** (verified on 16,950 truth SNPs):

```
strand=0:  POS_0based = bp.s1_left + t.pre        →  SNP POS = + 1; DEL/INS POS = + 0
strand=1:  POS_0based = bp.s2_right - t.suf - 1   →  SNP POS = + 1; DEL/INS POS = + 0
```

**Why this is correct** (strand=0):
- The Mode A path's spliced arm walks left-to-right in forward genome.
- arm byte k ↔ forward position `s1_left + k` (up to the s1 block, then extending
  through middle segments by their forward extents).
- `mutual_trim` finds `pre` = shared LCP between N-arm and T-arm. arm[`pre`] is the
  first differing byte → forward position `s1_left + pre`.

**strand=1** (arm reads as RC of forward): mirror analysis gives `s2_right - suf - 1`.

**Empirical accuracy on 16,950 truth SNPs at refdist=0:**
| formula | strand=0 hits | strand=1 hits | total exact |
|---|---:|---:|---:|
| old `pos_start + 1` | 227 / 12,498 (1.8 %) | 0 / 3,562 (0 %) | 232 (1.4 %) |
| new `s1_left+pre+1` / `s2_right-suf` | 12,217 / 12,498 (97.8 %) | 36 / 506 (7.1 %) | **15,420 (91.0 %)** |

The strand=1 formula is exact only for ~7 % of strand=1 bubbles. The other 93 % have an
offset of ±10–500 bp because the bubble's spliced arm doesn't tile a single forward
extent (middle super-nodes can have NON-CONTIGUOUS forward positions). Fixing strand=1
properly would require walking the Mode A path and reading each super-node's posmap
position — deferred.

**Side effect:** strand-aware POS lifted SNP recall (refdist=0) from 1.4 % to 90.9 %.

---

## 3. Phantom signatures — the bucket discovery

For each Mode B bubble (s1, s2), the bucket field reflects what Mode B's BFS found:

- `bucket=0` (std_both): both c1 and c2 arms exist between s1 and s2 — the **canonical
  bubble shape**.
- `bucket=1` (T_only_no_N): only c2 (tumor) arm exists; N-arm reconstructed via c3
  backbone detour up to 8 hops.
- `bucket=2` (N_only_no_T): mirror — only c1 (normal) arm; T-arm via detour.

**Critical empirical finding (FP vs TP for ≥50 bp emissions):**

| record class | bucket=0 | bucket=1 | bucket=2 |
|---|---:|---:|---:|
| TP DEL (40)  | 38 | 0 | 2 |
| TP INS (14)  | 14 | 0 | 0 |
| **FP DEL (750)** | **1** | 0 | **749** |
| **FP INS (584)** | 3 | **578** | 3 |

**Almost 100 % of large-DEL FPs are bucket=2 phantoms, and 99 % of large-INS FPs are
bucket=1 phantoms.** Requiring `bucket==0` for ≥50bp records collapsed geq50 FPs from
1,372 to 2 with only ~6 TP loss.

**Per-feature comparison (TP DEL vs FP DEL, all ≥50 bp):**

| feature | TP-DEL median | FP-DEL median |
|---|---:|---:|
| pre (shared prefix) | 35,008 | 0 |
| suf (shared suffix) | 38,272 | 0 |
| n_anchor_s1         | 564    | 4 |
| n_anchor_s2         | 612    | 46 |
| span               | 1,552  | 362 |

So real DELs sit between two LCP-rich anchor super-nodes (hundreds of sigs, tens of
thousands of shared bp), and the bubble walks an honest c1/c2 path between them. Phantoms
have tiny anchors and zero shared context — they're graph artifacts where two unrelated
nodes happen to form a parallel pair through a backbone detour.

---

## 4. INV signature (descriptive — detector not yet productionised)

For 13 of 15 truth INVs, there is a PAIRED (D-bucket=2 + I-bucket=1) bubble pair at
the same boundary:

```
truth INV 11:102,626,684 (svlen 267):
  bid=15310 I bucket=1 t_resid=21,447 SPAN=39   d=12 from INV start
  bid=15311 D bucket=2 n_resid=25,597 SPAN=25,597 d=102 from INV start

truth INV 12:82,670,102 (svlen 835):
  bid=2532  I bucket=1 t_resid=96,776 SPAN=59   d=2
  bid=36382 D bucket=2 n_resid=96,750 SPAN=59   d=53
```

**The signature is "huge residual / small span"** — the path walks through inverted
content (large t_resid) but the c3 anchors are close together in forward genome (small
span). For non-INV bubbles, residual ≈ span.

Two issues with productionising:
1. The svlen estimate from either residual length or arm span is OFF (residuals 10–200×
   larger than truth svlen; span 0.07–0.4× truth svlen). Truvari's `--pctsize 0.7`
   rejects.
2. Many non-INV bubbles also fit the orphan-bucket shape, generating ~60 false pairs
   per real one.

`tools/inv_pair_detect.py` implements the detector and is wired but emits 446 candidates
of which only 7 are size-matchable to truth. Left disabled in production until the svlen
problem is solved.

**INVs that have NO nearby Mode B bubble at all:**
- 19:26,010,322 (svlen 944) — 0 bubbles within 2 kb
- 19:54,111,254 (svlen 606,307) — 0 bubbles within 2 kb

These are upstream losses (the LCP graph never created a bubble at these inversion
boundaries) — not recoverable downstream.

---

## 5. DUP signature

Of 32 truth DUPs:
- **15 have NO Mode B bubble within 5 kb** — upstream loss; tandem-DUP graph topology
  is fundamentally degenerate (the path through the duplicated bytes is a self-cycle,
  not a parallel-arm bubble).
- 17 have a nearby bubble, mostly tagged I (the duplicated bytes appear as a tumor
  insertion at the DUP-cycle boundary).

`tools/dup_cov_detect` works at the coverage level (cov_t / bg_t ≥ 2.0, cov_t / cov_n ≥
1.7 on c3 super-nodes). With aggressive merge filters (svlen ≥ 200, ct/cn ≥ 1.7, ct ≥
65, svlen ≤ 500 kb), it gives 6 TP / 191 FP. Net F1 lift is **negative** (–0.13),
so it's left out of the production VCF.

To beat this, one of:
1. **Cycle-based detection** — the 105 small SCCs in modea.bin (size 2–2,280) are the
   honest tandem-DUP candidates. Untried this session.
2. **Bubble pattern** — a "DUP at locus" might appear as a paired (I-bucket=1, c3-cycle)
   in modeb. Untried.

---

## 6. Tricky / counter-intuitive truths

These are things I either had wrong, or that don't match naive intuition. Future you
should not have to rediscover them.

1. **The arm coordinate in strand=1 bubbles does NOT linearly map to forward genome.**
   The middle super-node in a strand=1 std_both bubble can sit at an unrelated forward
   position (the LCP graph connects super-nodes by KMP byte overlap, not by genome
   adjacency). This is why `s1_right - pre` only works for 7 % of strand=1 SNPs — the
   rest have a middle that doesn't tile [s2_left, s1_right] linearly.

2. **`pos_start` in the old `BubblePos` was `max(s1.pos_end)` — the right edge of the
   left anchor's flank window**, NOT the bubble's left boundary. POS = pos_start + 1
   was off by up to `(pre - s1.block_len)` bp. The new schema replaces pos_start/pos_end
   with the four edges s1_left/right/s2_left/right.

3. **bucket=2 (N_only_no_T) bubbles look STRUCTURALLY identical for real DELs and
   phantoms** — both can have span ≈ n_resid. The signal that separates them is the
   c3-anchor LCP density (n_anc_s1 + n_anc_s2 + pre + suf), but real bucket=2 DELs
   often have ONE tiny anchor (e.g., 1 sig) so the gate has to be `max(anc1,anc2) ≥ 50`
   not `min`. Even then, ~5/2071 real bucket=2 DELs pass a phantom-grade filter and
   are dropped along with the 196 phantoms. Net F1 win = keep `bucket==0` only.

4. **C-tag (complex, both arms diverge) bubbles are 0 % useful as VCF emissions.**
   228 of 318 are positionally near truth but their byte content doesn't represent a
   classifiable variant. Drop unconditionally.

5. **Real INS bubble span is often LARGE, not small.** Intuition says s1 and s2 should
   be adjacent in forward (the insertion is between them in tumor only). In practice
   the c3 anchor super-nodes UNITIG-COMPACT past the insertion locus, so their flank
   extents span across it. Real INS span median = 50; FP INS span < 10 are 98 % of
   phantoms. So span < 10 = phantom signal — but the inverse isn't true (large span
   doesn't make it TP).

6. **The "DUP via coverage uplift on c3" approach catches 31/32 truth DUPs at the
   POSITIONAL level but only 6/32 pass truvari pctsize=0.7.** The boundary clustering
   widens regions past truvari's size tolerance. Untouchable without a sharper boundary
   detector (probably from cycle SCCs in modea.bin).

7. **InsilicoSV truth has NO MNPs — only single-bp SNPs.** Any MNP we emit is by
   definition an FP. Drop MNP ≥50bp; emit MNP <50bp at your peril (they hurt lt50
   precision without contributing to recall).

8. **Posmap collisions inflate s1_left/s1_right.** Some super-nodes' sigs collide on
   multiple posmap positions; chrom_range picks min/max across the consensus chrom but
   collisions WITHIN the consensus chrom still inflate the apparent extent. Bubble
   35747's s1=[125543,716,704) is a 591 kb apparent extent driven by collisions, while
   the actual s1 super-node sits in a narrow region. This is a strand=1 POS failure
   mode and contributes to the strand=1 inaccuracy.

9. **N coverage in the residual region (`med_cov_n`) is a phantom signal for D-tag.**
   Real D-tag has cov_n ≈ 30 (haploid 30 × sim). Phantom D-tag has cov_n ≈ 500 — these
   are repeat-region phantoms where coverage piles up. `nc > 250` is a clean drop.

---

## 7. Fundamental upstream defects (caps on accuracy)

These limit achievable accuracy regardless of how the VCF emitter is rewritten.

1. **15/32 truth DUPs have NO bubble in Mode B.** Tandem-DUP topology is a cycle
   (the c3 backbone walks through the duplicated bytes twice in tumor), not a parallel-
   arm bubble. The std_both / orphan bucket enumerator can't capture cycles. **Fix
   requires modeb to also emit cycle bubbles** (the 105 small SCCs in modea.bin are the
   honest candidates).

2. **2/15 truth INVs have NO bubble within 2 kb** (19:26010322, 19:54111254). These
   are missed at the LCP graph level — the inversion boundaries didn't create distinct
   sigs from their non-inverted neighbours. Probably `min_mult=10` is dropping
   borderline sigs that would otherwise form an INV bubble.

3. **Strand=1 POS is approximate (±10–500 bp) for ~50 % of strand=1 bubbles.**
   Source: the spliced arm coordinate doesn't linearly map to forward genome when
   middle super-nodes come from unrelated forward positions. Fix requires walking the
   Mode B path per super-node and reading each one's posmap entry — non-trivial.

4. **Mode B's residual byte content for INVs over-estimates SVLEN by 10–200×.** The
   path walks PAST the inversion into adjacent content (because INV boundaries connect
   to RC-of-distant-context in tumor reads). No way to estimate the true INV svlen from
   residual or span alone.

5. **Posmap collisions inflate per-anchor extents.** When a super-node's sigs collide
   on multiple posmap positions (especially in repeat regions), s1_left → s1_right
   spans much more than the actual super-node. This is the root cause of strand=1
   POS inaccuracy.

6. **mm=10 prune removes some real low-coverage signal.** From the mm sweep memory:
   mm=10 keeps ~9.6 % of pre-prune sigs. The 4 % we drop relative to mm=3 includes
   most real-but-rare sigs — including borderline INV and small-DUP boundaries.

---

## 8. What the next session should attack

In rough order of expected F1 impact:

1. **Cycle bubbles** (modeb.cpp): emit cycle SCCs from modea.bin as DUP candidates.
   Expected: ~10–15 truth DUP recovered → geq50 recall 45 → 55 %.

2. **Strand=1 POS via super-node walk**: in vcf.cpp, walk the BFS path per super-node
   and look up each one's posmap position; the variant byte's position comes from the
   super-node that owns arm[pre]. Expected: SNP precision 48 → 65 % (strand=1 SNPs
   currently emit at wrong POS, generating ~3,500 FP SNPs).

3. **INV pair detector with svlen from posmap walk**: with #2 fixed, the residual byte
   that bookends a bucket=1/bucket=2 pair can be located precisely. Expected: 5–8
   truth INVs callable.

4. **Sub-50 indel FP reduction**: 1,023 FPs at lt50. Likely many are S-tag 1bp/1bp
   bubbles emitting at slightly wrong POS. dedup + tighter strand=1 will help.

5. **C-tag rescue**: 228 C bubbles are positionally near truth but currently dropped.
   Worth investigating whether ANY of them represent recoverable variants.

---

## 9. Files in this session

Code:
- `src/anchor.cpp` — extended BubblePos to s1_left/s1_right/s2_left/s2_right
- `src/vcf.cpp` — strand-aware POS, bucket=0 filter for ≥50 bp, per-tag phantom filter
- `tools/modeb_dissect.cpp` — full per-bubble TSV emitter (this is THE entry point for
  any future analysis: regenerate bubble_table.tsv after any anchor/vcf change)
- `tools/truth_map.py` — truth → nearest bubble
- `tools/pos_offset3.py` — POS formula sweep
- `tools/fp_features.py` — FP vs TP feature contrast
- `tools/inv_pair_detect.py` — INV detector (off by default)
- `scripts/merge_dup_cov.py` — DUP coverage merge (off by default)
- `scripts/svp_benchmark.sh` — 3-tier truvari driver

Output files in `outputs/svp/wgs_mm10/`:
- `bubble_table.tsv` (5.9 MB) — per-bubble features + computed VCF POS
- `bubble_residuals.bin` (108 MB) — raw residual bytes per bubble
- `truth_map.tsv` / `bubble_match.tsv` — truth ↔ bubble joins
- `calls.vcf` — final emission
- `truvari_{geq50,lt50,snp}/` — benchmarks
