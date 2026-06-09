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

## 6½. Findings from the cycle / loop / orphan-bucket follow-up push

After the initial pass closed at F1 0.618 on geq50, this section addresses three
questions the user raised:

### (a) Does Mode B "filter" the graph?

YES. Verified by re-running modeb with the actual counts:
```
[modeb] enumerated=37,654 (std_both=33,177 T_only_no_N=2,899 N_only_no_T=1,578)
        direct=1,550,621 cycle=4 multi_path=0 trunc=0 no_path=0
```
- **`direct=1.5M`** — c3 super-nodes connected directly to another c3 with no c1/c2
  branch (pure backbone segment — no variant). Dropping these is correct.
- **`cycle=4`** — bubbles whose c1/c2 DFS branch revisited a node. Tiny number because
  `dfs_branch` only follows c1/c2 nodes; c3-only cycles are never even attempted.
- **Mode B does NOT enumerate cycles through c3 super-nodes at all.** The
  enumerator's outer loop is `for each c3 s1: for each non-c3 out-neighbour x:
  dfs_branch(s1, x ...)`. Pure c3→c3 cycles (which is what tandem DUPs would form)
  are invisible to modeb. **This is the upstream blocker for DUP detection.**

### (b) Truth DUPs in the Mode A super-graph

`tools/truth_dup_locate` indexes each truth DUP's posmap region back to the super-nodes
that overlap it. Result (all 32 truth DUPs):

| metric | value |
|---|---|
| avg super-nodes per truth DUP | 15.1 |
| total c3 | 476 (98 %) |
| total c1 | 8 |
| total c2 | 0 |
| in giant SCC | 484 of 484 (100 %) |
| in small SCC | 0 |

**Every truth DUP's super-nodes are in the giant SCC** (96.7 % of all c3 super-nodes
are in this one SCC). The earlier `tools/dup_cycle_detect` looked only at SCCs of size
2-50 (84 candidates, 0 truth covered) — they're the wrong place to look.

### (c) Truth INVs in the Mode A super-graph

Same probe with `--svtype INV` for the 15 truth INVs:

| metric | value |
|---|---|
| avg super-nodes per truth INV | 37 |
| total c3 | 491 |
| total c1 | 66 |
| total c2 | **0** |

**Truth INVs have ZERO c2 super-nodes.** This kills every byte-level INV detector
that assumes tumor-only c2 sigs exist at the inversion. Mechanism (confirmed): the
caller does **dual fw+RC LCP extraction per read** (per `project_svp_read_pipeline`).
Tumor reads at an INV locus contain the inverted sequence X. Forward extraction sees
X. RC extraction sees RC(X) = the un-inverted sequence X'. Normal reads at the same
locus contain X'. Forward extraction sees X', RC extraction sees X. So **both samples
produce both X and X' at the INV locus** — no tumor-only sigs exist. INV is invisible
to c1/c2 colouring.

### (d) Phantom-bucket=2 origin — full descent

Four independent probes, building each on the previous:

**(d.1) Path topology** — `tools/phantom_descent` replays modeb's BFS for every
bucket=2 D-tag ≥50 bp bubble (1,111 total: 745 real + 366 phantom) and aggregates
features over the actual N-arm path (s1 → c1 middle → s2):

| feature | REAL (745) | PHANTOM (366) |
|---|---:|---:|
| path_len (median/p10/p90) | 3 / 3 / 3 | 3 / 3 / 3 |
| n_c1 in interior | 1 | 1 |
| n_c2 in interior | 0 | 0 |
| n_total_chain_ids in path | 2 / 1 / 3 | 2 / 1 / 9 |
| n_chain_ids with t=0 | EQUAL to total | EQUAL to total |
| n_chain_ids with t≥1 | 0 | 0 |
| median normal-count (n) | 25 | 24 |
| median tumor-count (t) | 0 | 0 |
| pre + suf | 0 / 0 / 0 | 0 / 0 / 0 |

**Real and phantom bucket=2 D-tag bubbles have IDENTICAL bubble-level structure**:
3-node path (s1, single c1 middle, s2), every chain_id (= original LCP node) in the
c1 super-node has tumor count t=0 and normal count n≈25. No bubble-level feature
separates them.

**(d.2) c2 sig hypothesis** — tested in `tools/region_color_census`. Hypothesis: in
phantoms, tumor reads at the c1 locus produce slightly different sigs that become
tumor-only (c2) super-nodes co-located with the phantom c1. Result: **0 c2 super-nodes
in the ±2 kb region** of EITHER class (745 real + 366 phantom). Tumor reads at phantom
loci do NOT produce co-located c2 sigs. Hypothesis FALSE.

**(d.3) The actual signal — repeat-region density** — same census reports:

| feature | REAL median / p10 / p90 | PHANTOM median / p10 / p90 |
|---|---:|---:|
| n_sigs in ±2 kb region | 126 / 113 / 10,806 | **12,080** / 171 / 174,926 |
| n_super_in_region | 16 / 6 / 1,415 | **573** / 15 / 11,133 |
| n_c3_super | 13 / 4 / 1,409 | **571** / 12 / 11,111 |
| n_c1_super | 2 / 2 / 7 | 5 / 2 / 183 |
| n_c3_node_sigs (chain occurrences) | 122 / 108 / 10,789 | **9,699** / 165 / 172,068 |

**Phantoms sit in repeat-saturated regions** — ~40× more c3 super-nodes and ~80× more
c3 chain occurrences in the surrounding ±2 kb than reals.

**(d.4) ROC threshold sweep on n_c3_super (drop if > T)**:

| T | REAL kept | PHANTOM kept |
|---:|---:|---:|
| 5 | 184 (25 %) | 1 (0 %) |
| 10 | 337 (46 %) | 21 (6 %) |
| 50 | 603 (82 %) | 125 (38 %) |
| 200 | 640 (87 %) | 139 (43 %) |

At `n_c3_super ≤ 10`: 337 real, 21 phantom — looks great ROC-wise (94 % bubble
precision). **But enabling those 337 bucket=2 D-tag emissions in vcf.cpp only added
+2 TPs and +356 FPs to the geq50 truvari run**. Reason: "real" here means *within
500 bp of any truth variant*, not *near a truth DEL of similar size*. Most of those
337 "real-near-truth" bucket=2 bubbles sit beside INS/DUP/INV truth and emit DEL
records that match no DEL truth → FPs. The filter discriminates phantoms from
"near-any-truth" cleanly, but doesn't add geq50 DEL TPs.

**Mechanism (best inference, NOT yet directly verified)** — in repeat-saturated
regions there are many closely-related LCP sigs at the locus. LCP boundary
determinism on slightly different read offsets / soft-error context selects one of
those alternatives. Normal reads happen to pool onto sig X (high n). Tumor reads
pool onto sigs Y, Z (different selections), leaving X at t=0. X then becomes a
c1 (normal-only) super-node despite the bytes existing in both samples. This
isn't sequencing error per se — both samples have the bytes; they just hit
different LCP determinism arms. Verifying requires replaying LCP on the raw
fastq reads at one locus and comparing per-strand sig selection between samples —
not done this session.

**Files used in the descent:**
- `tools/phantom_descent.{cpp,py}` — per-bubble TSV with (n,t) for every chain_id
- `tools/region_color_census.cpp` — per-bubble super-node colour count in ±2 kb
- `tools/bucket2_origin.cpp` — posmap hit rate of c1 nodes (disproved the
  "sequencing error" hypothesis: real 100 %, phantom 99.9 %)

### (e1) DUP via "coverage candidate + graph cycle confirmation"

`tools/dup_cov_graph` implements the directive: find cov-uplift clusters → induced-
subgraph cycle confirmation. Results across two param regimes:

| params | clusters | confirmed | truth covered |
|---|---:|---:|---:|
| ratio_bg 1.5 ratio_tn 1.3 cluster_kb 2 min_size 3 | 592 | 8 | 0 / 32 |
| ratio_bg 1.3 ratio_tn 1.2 cluster_kb 10 min_size 2 | 646 | 19 | 0 / 32 |

The induced-subgraph cycle check is too strict. Real-DUP super-nodes are connected
sequentially within the cluster (s1→s2→s3...) but their back-edges go OUT of the
cluster (through hub nodes in the giant SCC) and back in. The induced subgraph has
forward edges only → no cycle within it → rejected. Loosening to "any cluster member
has a cycle through ≤K hops touching another cluster member" recovers the dup_loop_
detect behavior (11/32 covered at 4.5 % precision, F1 still negative).

### (e) Cycle-based DUP detection (loose vs. strict)

`tools/dup_loop_detect` finds the shortest cycle through each cov-uplift c3 super-node.

| params | candidates | near-truth | size-matched | distinct truth covered |
|---|---:|---:|---:|---:|
| ratio_bg 2.0, hops 15, **all-uplift** | 295 | 78 | 13 | 3 / 32 |
| ratio_bg 1.7, hops 15, all-uplift | 28 | 3 | 0 | 0 / 32 |
| ratio_bg 1.5, hops 30, no uplift gate | 7,201 | 511 | 325 | 11 / 32 |
| ratio_bg 1.5, hops 30, **60 %-uplift** | 33 | 3 | 0 | 0 / 32 |

The loose run covers 11 of 32 truth DUPs at size-match level but the precision is
325 / 7,201 = 4.5 %. Adding to VCF: net F1 −0.44 (recall +6 %, FP +680). Not deployed.

### (f) Back-edge DUP detection

`tools/dup_backedge_detect` looks for direct edges `u → v` where both are c3, anchored
on the same chrom, with `u.pos > v.pos + min_svlen` (a "backward jump" in forward
coordinates). Result: **1 candidate, 0 near any truth DUP.** Real DUP cycles don't
manifest as a single back-edge — they're multi-hop cycles inside the giant SCC.

### (g) `[modeb] direct=1.5M` is correct to drop

The 1.55 M `direct_only` entries are c3 super-nodes directly connected by a c3→c3 edge
with no c1/c2 branch between them. These represent pure backbone segments where the
graph has no variant signal between the two anchors. Not a variant; correctly dropped.

---

## 6¾. Raw-reads-and-pre_prune descent (the round 3 deep dive)

This section reports what was measured by going BELOW Mode A — directly into the
pre-prune sig pool (`outputs/svp/wgs/pre_prune.bin`, 887,246,413 raw `SigCount` records,
14 GB; the sample pool before the mm=10 row-level drop) and into the raw FASTQ reads.

### (h) pre_prune at 50 hand-picked loci

`tools/preprune_at_locus` streams pre_prune.bin and accumulates `(n, t)` for every sig
that anchors in each query locus via posmap. The query set is **10 phantom + 10 real
bucket=2 D-tag ≥50 bp bubbles + 10 truth DUP centers + 10 truth INV centers + 10
control positions on chr1 ≥50 kb from any truth** — all generated by
`tools/make_locus_queries.py`. The full table is in `outputs/svp/wgs_mm10/preprune_at_locus.tsv`.

Median total_t / total_n in pre_prune at each locus class (raw, NOT post-prune):

| class | locus median t/n | observation |
|---|---:|---|
| control (clean chr1) | **0.91** | tumor ≈ normal, as expected for no variant |
| truth DUP | **2.07** | tumor ≈ 2× normal — the DUP coverage signature |
| truth INV | **0.94** | tumor ≈ normal (INV doesn't change coverage, only orientation) |
| real bucket=2 D-tag | 0.77 | mixed; query window includes anchor regions |
| phantom bucket=2 D-tag | 0.59 | tumor LESS than normal at the locus |

Per-DUP detail (all 10 truth DUPs sampled): ratios were 1.20, 1.37, 1.69, 1.80, 1.98,
2.04, 2.22, 2.64, 2.80, 2.99 — every single one elevated.

### (i) The actual mechanism behind phantoms

At phantom locus #0 (`phantom_0_bid1780`, the bucket=2 D-tag bubble at chr 11 ≈ 102626666):
- n_ref_sigs = 13 reference sigs anchor in the locus
- total_n = 280, total_t = 107 — **tumor reads DO exist at the locus** (tumor coverage ≈ 38 % of normal)
- BUT: **9 of those 13 sigs have t = 0** in pre_prune (normal-only) — the c1 super-nodes
- Only 4 sigs reach mm=10 in tumor

**So the phantom mechanism is NOT "tumor has no reads at locus".** Tumor has reads;
those reads just produce DIFFERENT sigs than normal. The reference sigs that survive
in the normal pool don't appear in the tumor pool. The pre_prune verification rules
out the prior hand-wave "phantoms are pruned-away tumor counts" — the tumor counts on
the reference sigs are LITERALLY 0 even BEFORE the mm prune.

**Why does tumor select different sigs than normal at the same locus?** LCP boundary
determinism on slightly different read context: HiFi/CCS error patterns, read start
offset, RC orientation. In repeat-saturated regions where many sigs share the same
posmap target window, tiny input differences pick different cores. Normal-pool and
tumor-pool sigs at the locus are both reference-derived (per `bucket2_origin`) but
non-overlapping sets.

### (j) Genome-wide pre_prune binning + truth DUP coverage

`tools/preprune_genome_bins` walks pre_prune.bin and bins (n, t) per 1000 bp posmap
window. Output: `preprune_bins.tsv` (3,042,693 bins, 120 MB).

Sweep on bin-level t/n ratio with cluster constraint (consecutive bins above threshold):

| T (ratio) | total bins ≥ T | bins in truth DUP | distinct DUPs covered |
|---:|---:|---:|---:|
| 1.3 | 460,810 | 542 | 21 / 32 |
| 1.5 | 184,236 | 486 | 18 / 32 |
| 2.0 | 16,882  | 291 | 14 / 32 |
| 2.5 |  1,789  | 104 |  ? |
| 3.0 |    220  |  23 |  ? |

Cluster-based DUP calling (`tools/dup_preprune_call.py`):

| T | min_consec | clusters | size-matched | F1 if merged |
|---:|---:|---:|---:|---:|
| 1.7 | 3 | 8,901 | 9 / 32 | net F1 NEGATIVE (≈ 0.02) |
| 2.0 | 5 | 1,013 | 6 / 32 | net F1 NEGATIVE |
| 2.5 | 5 |    88 | 2 / 32 | net F1 NEGATIVE |

The signal is real (truth DUP regions ARE bin-elevated) but the FP base rate at the
genome-bin level is too high — repeat regions and mappability-uneven loci show similar
elevation without being DUPs. Even strict thresholds give ≥30 FPs per TP. Not deployed.

### (k) The INV in pre_prune: n_t_zero spike

For the **truth INV at chr 11:102626684-102626950**, the surrounding bins read:

| bin (chr 11) | total_n | total_t | n_t_zero |
|---|---:|---:|---:|
| 102625000-102626000 | 644 | 882 | 0 |
| **102626000-102627000** | **700** | **665** | **10** |
| 102627000-102628000 | 693 | 920 | 0 |
| 102628000-102629000 | 772 | 976 | 0 |

The INV creates a **single-bin spike of n_t_zero = 10** (10 sigs are normal-only) flanked
by zero-spike bins. This is the INV signature in pre_prune. Same locus as phantom #0 →
that "phantom" was actually the INV.

**But the n_t_zero signal alone is not deployable as an INV detector**: 164 genome-wide
"isolated" bins (n_t_zero ≥ 5 with ≤ 2 in both surrounding bins) — 0 of 164 overlap any
of the 15 truth INVs by the simple position check. The signal exists locally at every
INV but is shared with phantoms, repeat-region asymmetries, and small DELs across the
genome. Without a way to separate the per-locus shapes, n_t_zero alone can't drive
calls.

### (l) Reclassification of phantom_descent with bubble center

The original phantom_descent used **`s1_left`** as the bubble's position for truth-
proximity classification. But s1 is a unitig-compacted chain that can span 25+ kb
(verified bid 1780: s1 spans chr 11:102601284-102626666 = 25 kb). The bubble's actual
variant locus is `min(s1_right, s2_left)` — the c1 super-node's position.

Reclassifying with `locus = min(s1_right, s2_left)`:

| old | new | count |
|---|---|---:|
| phantom → phantom | 113 |
| phantom → real | **212** |
| real → phantom | 234 |
| real → real | 499 |

Of the 212 "newly-real" bucket=2 D-tag bubbles (previously misclassified as phantom),
their nearest truth SVTYPE breakdown: **53 DEL, 15 INS, 32 INV**, 611 near SNP/small.

**Translation: 32 bucket=2 D-tag bubbles sit at truth INV positions but emit DEL records
that don't match any DEL truth → they're all FPs that could be re-routed to INV.**
This is the most concrete actionable gap I found in this session.

### (m) Raw read extraction at the phantom_0 locus

`tools/extract_locus_reads` finds reads whose canonical 21-mer hit count against the
locus's 457 reference k-mers passes a threshold. At the loose threshold (min_hits=8):
**all 200 extracted reads from each sample came from chr 1**, not chr 11. The phantom
locus shares 21-mers with a chr 1 repeat family — direct read-level evidence the
phantom sits in cross-chromosome repeat content.

Strict-threshold extraction at min_hits=50 ran for 14+ minutes scanning ~14 GB of the
29 GB FASTQ without producing any matches, then was killed. Either no reads with ≥50
distinct 21-mer hits exist in the first 14 GB (chrom-sorted FASTQs would put chr 11
reads near the middle), or the cross-chrom repeat sharing is so high that 50 unique
21-mer matches just don't occur for any single read in the sampled fraction. The
matcher logic was verified at min_hits=8. Re-running at min_hits=20 with a higher
read budget is the next step (not done this session for time reasons).

### (n) Deployable INV detector via bin n_t_zero spike — measured F1 impact

`tools/inv_from_bins.py` flags bucket=2 D-tag ≥50 bp bubbles whose locus bin has
n_t_zero ≥ T and surrounding bins have n_t_zero ≤ N. Wired into `src/vcf.cpp` via
`--inv-flags`; bucket=2 D-tag bubbles with INV flag emit as `<INV>` instead of being
dropped. svlen heuristic = bin width = 1000 bp.

Measured truvari deltas (geq50):

| T | N | emit | truvari TP | FP | per-SVTYPE   | F1 |
|---:|---:|---:|---:|---:|---|---:|
| 8 | 2 | 13 | 52 (DEL 38, INS 13, INV 1) | 14 (1 DEL, 1 INS, 12 INV) | 0.584 |
| 10 | 0-2 | 4 | 51 | 4 | unchanged base, F1 = 0.604 |

The signal CAN call 1 INV correctly (chr 19:26010322 svlen 944 ≈ bin width 1000) and
incidentally bumps DEL TPs by 2 (re-routed bucket=2 bubbles), but the 12 INV FPs net
hurt F1 (0.618 → 0.584). NOT DEPLOYED. The bin-width svlen heuristic is too coarse
for the truth INVs at 200-600,000 bp range.

---

## Files added this round

C++:
- `tools/preprune_at_locus.cpp` — stream pre_prune.bin, sum (n,t) per query locus
- `tools/preprune_genome_bins.cpp` — genome-wide pre_prune binning
- `tools/extract_locus_reads.cpp` — k-mer-based FASTQ subset by locus
- `tools/relcp_at_locus.cpp` — re-run LCP on a FASTQ subset
- `tools/modea_dump_region.cpp` — modea topology dump per region
- `tools/phantom_descent.cpp` — bucket=2 D-tag per-bubble feature dump
- `tools/region_color_census.cpp` — c1/c2/c3 super-node count per bubble region
- `tools/bucket2_origin.cpp` — c1 path posmap hit rate
- `tools/truth_dup_locate.cpp` — per-truth-DUP / INV super-node enumeration
- `tools/dup_cov_graph.cpp` — DUP detector with cycle-confirmation
- `tools/dup_loop_detect.cpp` — DUP detector via shortest cycle through cov-uplift

Python:
- `tools/make_locus_queries.py` — generate phantom/real/DUP/INV/control queries.tsv
- `tools/dup_preprune_call.py` — DUP candidates from preprune bins
- `tools/inv_from_bins.py` — INV candidates from pre_prune n_t_zero spike

Outputs of interest:
- `preprune_at_locus.tsv` — 50 hand-picked loci, full pre_prune detail
- `preprune_bins.tsv` — 3.04 M genome-wide 1 kb bins
- `phantom_descent.tsv` — 1,111 bucket=2 D-tag bubble features
- `region_color_census.tsv` — 1,058 bubble region color counts
- `modea_dump_region.tsv` — 1,801 super-node features at 49 queries
- `bubble_table.tsv` — 37,654 bubble features (regenerated each anchor/vcf change)

---

## 7. Fundamental upstream defects (caps on accuracy)

These limit achievable accuracy regardless of how the VCF emitter is rewritten.

1. **Mode B explicitly does not enumerate c3-only cycles.** `modeb.cpp::enumerate_bubbles`
   walks `for each c3 s1: for each non-c3 out-neighbour: dfs_branch(...)`. Pure c3→c3
   cycles (= the topology of every tandem DUP) are never reached. Verified by re-running
   modeb: cycle_skipped=4 (only the c1/c2 DFS-revisit kind). **ALL 32 truth DUPs sit
   inside the giant SCC of the super-graph** (96.7 % of c3 nodes are in this one SCC,
   per `tools/truth_dup_locate`). Looking at small SCCs (the 84 found by Tarjan size 2-50)
   produces 0 truth coverage. Two independent attempts (`tools/dup_cycle_detect`,
   `tools/dup_loop_detect` with shortest-cycle-through-cov-uplift-node, several param
   regimes) all give either ~3 % precision (loose) or 0 distinct truth covered (strict)
   — net F1 hurt either way. **Fixing DUP needs modeb to enumerate c3-only cycles with
   their own bucket so downstream can filter by bubble-grade features rather than raw
   topology.**

2. **Truth INVs have ZERO c2 super-nodes.** Verified for all 15 truth INVs (491 c3 + 66
   c1 + 0 c2 across all of them) with `tools/truth_dup_locate --svtype INV`. Mechanism:
   the caller does dual fw+RC LCP extraction per read. Tumor's inverted bytes appear as
   X via fw extraction and as RC(X)=X' via RC extraction; normal at the same locus
   produces both X' and X via the same dual extraction. Result: both samples produce
   both forms → every sig at the INV locus is c3 (shared). c1/c2 colouring CANNOT
   detect INVs structurally. The four prior INV approaches that depended on c2 sigs
   (block-RC, position-overlap, RC-pairing, etc.) all fail for this reason; they're
   not detector bugs, the c2 sigs literally don't exist.

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
