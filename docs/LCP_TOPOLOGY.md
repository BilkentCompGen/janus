# LCP-core topology probes: bags, alphabet symmetry, encoding granularity, locality

A second analysis pass (after `docs/CORE_ANALYSIS.md`), aimed at the *structure* of the LCP-core map rather than its uniqueness/tangle statistics. Everything here is produced by **read-only experiment binaries** under `tools/` that link the compiled `build/{io,cores,MurmurHash3}.o` and call `lib/lcptools` — **no change to `src/` or `lib/lcptools/`**. Each tool sets the lcptools global alphabet before parsing, so it probes the *same* parser svp uses.

The unifying object is the **bag**: for a genome parsed into LCP cores, group cores by their lcptools `label` and collect the distinct core *sequences* under each label —

```
bag(label) = { distinct genomic substrings of cores carrying that label }
```

Bag size = how many genuinely different sequences one label conflates. A label with bag size 1 is a faithful sequence id (it may still occur many times — a repeat); a label with bag size ≥2 is a *collision*.

## How LCP cores work (the three facts every experiment leans on)

Read from `lib/lcptools/{lps.c,core.c,encoding.c}`:

1. **Level-1 parsing is purely order-based** (`parse1`, `lps.c:182`). It scans left→right and cuts a core at a local minimum (`a>b<c`), a local maximum (`a<b>c` plus the neighbour guards `prev≤a` and `c≥next`), or a run (`b==c`), comparing only the integer `alphabet[base]` values with `<`/`>`/`==`. So the level-1 partition depends only on the **total order** of the 4 symbols, never on the specific integers.
2. **The label is a lossy digest.** At level 1 (`core.c:126`) `label = (len-2, first, 2nd-last, last)` packed at `alphabet_bit_size` bits; at deeper levels (`core.c:219`) `label = MurmurHash3_32(child_labels[first, 2nd-last, last] + count, seed 42)`. **`ulabel` is `uint32`.** Deeper-level core *comparison* (used during deepening) uses the full `bit_rep`, not the label.
3. **Records are parsed independently**, separated by N-bumpers; `extract_cores` (`src/cores.cpp:20`) keeps only `{start,end,label}`. A core's sequence is `seq.substr(start,end-start)`.

## Tools (all in `tools/`, built via `make lcp-exp`)

| tool | experiment | one-line method |
|---|---|---|
| `lcp_bags` | 1 bags · 2 rotation · 4 ref-vs-mixed · I Zipf | bag builder; `--order` permutes the alphabet; `--dump` writes multi-seq bags |
| `lcp_triplet_bags` | 3 Huffman/triplet | rewrites genome into frequency-ranked k-mer symbols, runs LCP, bags |
| `lcp_collision_anatomy` | B fiber structure · G k-mer | separates repeat / structural / 32-bit-birthday collisions by shared-affix + Jaccard |
| `lcp_perturb` | E locality | flips/ins/del one base, measures cores changed + bp radius |
| `lcp_rc_invariance` | A strand symmetry | forward parse vs reverse-complement parse, mirror-matched |
| `lcp_crosslevel` | C hierarchy | L_low cores per L_high span + refinement determinism |
| `lcp_repeat_probe` | F repeats | synthetic homopolymer/tandem response curves |
| `lcp_consistency` | J local-consistency | local-window re-parse vs whole-record parse, interior cores |

**Scope.** Whole genome, both `simulation/ref.fasta` (3.21 Gbp, 455 rec) and `simulation/truth/mixed/out.fa` (3.22 Gbp; SVs in chr19 only), run separately. Default level 5 (svp production). Tools were validated on hand-checkable synthetic inputs first, then run whole-genome. **Measured cost:** one whole-genome L5 bag run = **94 s wall, 10.5 GiB RSS** (single-threaded); anatomy = 172 s. Far cheaper than the prior flank-sweep tool because each does one extraction + one sort.

---

## Exp 1 — Bags: how lossy is the label? (`lcp_bags`)

**Method.** Parse the genome at level L (default 5). For every pure-ACGT core (N-block/barrier cores are skipped and counted), emit `(label, murmur64(core_seq))`. Sort by `(label, seqhash)`; one linear pass yields, per label, the **bag size** (distinct seqhashes) and the **multiplicity** (total cores). Histogram both; dump every bag of size ≥2 (byte-verified distinct members) to a TSV.

**Result — reproduces `CORE_ANALYSIS.md` exactly** (validates the tool at whole-genome scale; differences are the ~2.4k skipped N-cores):

| level | ACGT cores | distinct labels | distinct seqs | label loss (seqs−labels) |
|---:|---:|---:|---:|---:|
| 4 | 106,165,111 | 93,181,736 | 97,592,479 | 4,410,743 |
| 5 |  45,413,056 | 42,479,424 | 43,072,014 |   592,590 |
| 6 |  19,444,226 | 18,541,782 | 18,698,717 |   156,935 |

(`CORE_ANALYSIS` reports label loss 4,411,146 / 592,851 / 157,112 — matched.)

**Bag-size and multiplicity distributions (ref, L5):**

```
BAG-SIZE (distinct seqs / label)  =1:42,064,593  =2:370,757  3-4:31,491  5-8:7,815  9-16:2,782  17-64:1,696  65-256:272  257+:18   max=683
LABEL-FREQ (cores / label, Zipf)  =1:41,092,352  =2:1,053,967 3-4:213,219 5-8:92,604 9-16:18,989 17-64:6,308 65-256:1,318 257+:667 max=2,552
```

**Reading.** 99.0 % of labels are bag-size 1 (faithful). But the tail is real: one label conflates up to **683 distinct sequences**, and the multiplicity Zipf tail reaches **2,552 cores on a single label** — heavy-tailed, repeat-driven. The label is an excellent-but-not-perfect sequence id; §Exp B dissects the 0.98 % that collide.

## Exp 4 / I — ref vs mixed, and the multiplicity spectrum

**Method.** Same tool, both genomes; compare the distributions. (Mixed differs from ref only in chr19 SVs, so whole-genome spectra are expected near-identical — the per-genome bag is the object of interest, not a diff.)

| genome (L5) | ACGT cores | distinct labels | collisions (bag≥2) | collision rate |
|---|---:|---:|---:|---:|
| ref   | 45,413,056 | 42,479,424 | 414,831 | 0.977 % |
| mixed | 45,521,029 | 42,495,588 | 414,253 | 0.975 % |

Mixed carries **+107,973 cores** — the chr19 duplicated/inserted content — yet the collision *rate* is identical to 3 decimals: the SV perturbation is real but chr19-local, drowned at genome scale (consistent with `CORE_ANALYSIS` §B.10). The core surplus is a **constant +0.238 % at every level** (+251,522 / +107,973 / +46,188 at L4/L5/L6), i.e. the injected SV sequence produces cores at the same density as the rest of the genome. The multiplicity histogram is the genome's LCP "vocabulary" spectrum: a Zipf-like power-law body (singletons dominate, then a smooth 1/freq tail) — the same heavy tail the `LABEL-FREQ` row shows.

---

## Exp F — Homopolymer & tandem response (`lcp_repeat_probe`, synthetic)

**Method.** Feed generated single-record buffers (homopolymers; perfectly periodic `(AC)ⁿ`,`(ACG)ⁿ`,`(ACGTTGC)ⁿ`; a fixed pseudo-random 300 bp unit tandem-duplicated K×) to `extract_cores` at levels 1–6; report cores, distinct labels, max label frequency, max core length.

**Result — three regimes:**

| input | L1 cores | L1 labels | L1 max-freq | L5 cores |
|---|---:|---:|---:|---:|
| homopolymer A×{100…100 000} | **0** | 0 | 0 | **0** |
| (AC)×100 000 | 99,999 | **1** | 99,999 | **0** |
| (ACG)×100 000 | 99,999 | **1** | 99,999 | **0** |
| (ACGTTGC)×10 000 | 19,999 | **2** | 10,000 | **0** |
| random-300bp × 64 | 8,830 | 59 | 448 | **254 (4 labels, freq 64)** |

**Reading.** (a) **Homopolymers are invisible to LCP** at every level — a run of equal symbols has no local extremum, so `parse1` emits nothing. (b) **Perfectly periodic tandems collapse to zero at deep levels**: at L1 every period is one core sharing a single label, and deepening (DCT compression of byte-identical cores) annihilates them by L5. (c) **A quasi-random tandem dup** of unit U×K yields, at L5, a *constant* ~4 distinct labels each repeated **exactly K times** — this is the tandem-DUP uniqueness wall (`CORE_ANALYSIS` §B.10) reconstructed from first principles: the interior cores of every copy are byte-identical, so no signature or context can separate them.

---

## Exp B / G — Anatomy of label collisions (`lcp_collision_anatomy`)

**Method.** At level 5, take every label with ≥2 distinct core sequences. The deep label is `hash(first, 2nd-last, last child label + count)`, so it ignores the *middle* of a core. Predicted two causes of collision: **structural lossiness** (two cores share boundary children → share a long bp prefix *and* suffix, differ only in the middle) and **32-bit murmur birthday collisions** (`ulabel` is uint32; at 45 M cores the birthday bound predicts ≈210–236 k accidental pairs — sequences that share *nothing*). Separate them by the shared prefix+suffix (bp) and the 8-mer Jaccard between the two most-distant members; dump labelled examples.

**Result (ref, L5):** of 414,831 collision labels —

| class | count | share | shared affix | 8-mer Jaccard |
|---|---:|---:|---|---|
| **structural** (affix ≥16 bp) | 201,828 | 48.7 % | 65–256 bp | ~0.88 |
| **birthday-like** (affix <16 bp) | 213,003 | 51.3 % | 0–2 bp | 0.000 |

```
shared-affix bp  =0:152,365  1-4:58,229  5-16:2,481 | 17-64:5,386  65-128:88,178  129-256:108,090  257+:102
8-mer Jaccard    =0:173,955  <0.1:38,372 | 0.1-0.3:3,398  0.3-0.6:16,466  >=0.6:182,640
```

**Reading.** Both signatures are sharply **bimodal** with almost nothing between the modes — collisions are cleanly one type or the other. Roughly **half of all label collisions are pure 32-bit hash accidents** (≈213 k, matching the birthday prediction), conflating totally unrelated sequences; the other half are structural, conflating near-identical sequences that differ only in their interior. This is the exact fiber structure of the label map, and it is *why* svp keys its signature on a 64-bit murmur of the window rather than the 32-bit label (`CORE_ANALYSIS` conclusions #3–4). Dumped examples (`outputs/lcp_topology/anatomy_ref_L5.tsv`): STRUCT rows share 120–150 bp affixes at Jaccard ~0.88; BIRTHDAY rows share ≤2 bp at Jaccard 0.000.

**Mixed is nearly identical** — 414,253 collisions, 48.5 % structural / 51.5 % birthday-like, same bimodal affix/Jaccard split — so the fibre structure is a genome-intrinsic property of the LCP+uint32-label map, not an artifact of the injected SVs.

---

## Exp 2 — Alphabet rotation (`lcp_bags --order`)

**Method.** `parse1` is order-only, so permuting which integer each base maps to is a clean action of the symmetric group S₄ on the partition. Re-run the bag builder under named orders (id `A<C<G<T`, rev, three cyclic shifts) at L5, and the **full 24-permutation orbit** at L6 (cheapest level), on ref and mixed.

**Result (ref, L5 core counts):**

| order (A,C,G,T codes) | meaning | ACGT cores | Δ vs id |
|---|---|---:|---:|
| id (0,1,2,3) | A<C<G<T | 45,413,056 | — |
| rev (3,2,1,0) | T<G<C<A | 45,258,184 | −0.34 % |
| cyc1 (1,2,3,0) | T<A<C<G | 45,545,094 | +0.29 % |
| cyc2 (2,3,0,1) | G<T<A<C | 45,325,273 | −0.19 % |
| cyc3 (3,0,1,2) | C<G<T<A | 45,912,227 | +1.10 % |

**Reading.** Every ordering yields a *different* partition (core counts span 45.26 M–45.91 M, a 1.4 % range), so the LCP partition is genuinely sensitive to symbol order — not merely relabelled. Crucially, **order-reversal is not a symmetry**: id (45.413 M) ≠ rev (45.258 M). This falsifies the natural "reversing the order just swaps minima↔maxima, same boundaries" guess — because `parse1`'s LMAX rule carries extra neighbour guards (`prev≤a`, `c≥next`) that LMIN does not, so the parse is **chirally asymmetric** in the value order.

**Full 24-permutation orbit (ref, L6):** all **24 orderings give 24 *distinct* core counts**, spanning 19,306,215 (`T<C<A<G`) to 19,693,010 (`G<C<A<T`) — a spread of 386,795 cores (**2.00 %**). Reversal pairs differ throughout (e.g. `0123`↔`3210`: 19,444,226 vs 19,376,336). So the symmetric group S₄ acts essentially **freely** on the partition — no two of the 24 orders share it, decisively refuting the reversal-collapse hypothesis. (Mixed orbit: same shape, +0.24 % cores; full table in `docs/analysis/lcp_bags_rotation.csv`.) Practically, this means the alphabet order is a real, ~2 %-magnitude knob on the core set — and a union over orderings would enlarge the anchor pool, though exploiting that is out of scope here.

---

## Exp 3 — Huffman / triplet first layer (`lcp_triplet_bags`)

**Method.** Replace the per-base 2-bit alphabet with a per-**k-mer** symbol whose integer code is its **frequency rank** (0 = most frequent) — the byte-aligned analog of a Huffman code. (A true variable-length Huffman *bitstream* has no byte-aligned symbols, so `parse1`'s per-symbol order is undefined on it — that variant is deliberately dropped.) The genome is rewritten frame-0 into single-byte symbols drawn from `{0x00–0x1F, 0x40–0x5F}` (all bit5=0, so the label's `&0xDF` mask stays an identity → labels remain correct), with a sentinel byte for non-ACGT; `alphabet_bit_size` is set to 2k. k=2 (dinucleotide, 16 symbols, 4 bit) is the de-risked case; k=3 (triplet, 64 symbols, 6 bit) is the headline. A diagnostic counts level-1 cores whose `(len-2)` exceeds the label's distance field (14 bits at k=3).

**Result (whole genome, L5):**

| encoding | genome | symbol buf | symbol-cores | distinct labels | collisions | rate | L1 dist-overflow |
|---|---|---|---:|---:|---:|---:|---:|
| per-base (Exp 1) | ref | 3.21 Gbp | 45,413,056 | 42,479,424 | 414,831 | 0.977 % | n/a |
| dinuc k=2 (4 bit) | ref | 2.00× shorter | 22,563,556 | 21,755,327 | 112,495 | **0.517 %** | 4 |
| triplet k=3 (6 bit) | ref | 3.00× shorter | 14,750,882 | 14,417,784 | 54,657 | **0.379 %** | 243 |
| dinuc k=2 | mixed | 2.00× | 22,615,712 | 21,798,510 | 112,608 | 0.517 % | 4 |
| triplet k=3 | mixed | 3.00× | 14,784,738 | 14,445,155 | 54,789 | 0.379 % | 242 |

**Reading.** A coarser, frequency-ranked first layer produces **fewer, longer cores** (45 M → 22.6 M → 14.8 M) and — notably — a **monotonically lower collision rate** (0.977 % → 0.517 % → 0.379 %): each k-mer symbol carries more bits (`alphabet_bit_size` 2→4→6), so the level-1 label is richer and conflates less, and the longer cores have more distinctive content. The level-1 distance-field overflow predicted for k=3 is real but tiny — **243 cores** out of 14.75 M (the centromeric/N-block satellite runs whose `(len-2)` exceeds the 14-bit field); bagging at L5 keys on the murmur label, not the raw level-1 formula, so these are isolated, not silently merged. Mixed tracks ref to 3 decimals (chr19-local SVs). Net: the Huffman-like layer is a viable coarsening that *improves* label fidelity at the cost of sparser genome coverage.

---

## Exp E — Point-mutation perturbation radius (`lcp_perturb`)

**Method.** This measures LCP's *locality* — the property that makes svp's anchoring variant-robust. For thousands of sampled genome positions, take a ±8 kb window, flip/insert/delete the centre base, re-parse both versions, and compute the longest common prefix and suffix of the two core-hash sequences. Cores outside that changed middle block are unperturbed; report the count of changed cores and the bp radius from the edit to the block edge, per level and mutation type. (Samples whose effect reaches the window edge are discarded.)

**Result (whole-genome ref, 4,744–4,774 usable samples each, 0 edge-truncated):**

| level | mutation | mean cores changed | mean radius | max radius |
|---:|---|---:|---:|---:|
| 4 | SNP | 2.27 | **57.2 bp** | 170 bp |
| 5 | SNP | 2.29 | **134.7 bp** | 748 bp |
| 6 | SNP | 2.31 | **315.3 bp** | 779 bp |
| 5 | INS | 2.27 | 133.7 bp | 346 bp |
| 5 | DEL | 2.29 | 134.5 bp | 748 bp |

**Reading.** The number of perturbed cores is **essentially constant — ~2.3 — across all levels and all mutation types**; a point edit never changed more than ~4 cores in any of ~24,000 trials. What grows with level is the *radius* (57 → 135 → 315 bp from L4→L5→L6), tracking the median core length (63/150/352 bp in `CORE_ANALYSIS` §A.2) — deeper cores see wider context, so the same handful of perturbed cores span more bp. Mutation type barely matters (SNP/INS/DEL all ~134 bp at L5). **This is the locality guarantee quantified**: anchors more than a few hundred bp from a variant are untouched, so svp's chain walks through them — exactly why LCP anchoring is variant-robust.

---

## Exp A — Reverse-complement invariance (`lcp_rc_invariance`)

**Method.** Parse the genome forward; reverse-complement each record and parse that forward; a forward core `[s,e)` in record `[rs,re)` mirrors to `[rs+re-e, rs+re-s)`. Report the fraction of forward cores with an exactly-mirrored RC core, under both the id ordering and the complement (reversed) ordering.

**Result (whole genome):**

| genome | level | match (id order) | match (complement order) |
|---|---:|---:|---:|
| ref | 4 | 0.41 % | 0.49 % |
| ref | 5 | 0.07 % | 0.08 % |
| ref | 6 | 0.01 % | 0.01 % |
| mixed | 5 | 0.07 % | 0.08 % |

**Reading:** forward LCP cores have *no* strand symmetry — the match is at chance and *falls* with level (0.41 %→0.07 %→0.01 %, as longer cores make an exact mirror-position coincidence rarer), and the complement ordering does not rescue it (0.49 % vs 0.41 % at L4 — negligible). Parsing the RC under complement-order reduces algebraically to parsing the plain *reverse* under id-order, and `parse1`'s directional LMIN/LMAX asymmetry (same cause as Exp 2) means it is not reversal-invariant. This is exactly why lcptools ships a separate backward `parse2`/`init_lps2` for the reverse-complement strand; svp, parsing only forward, sees the two strands as unrelated.

---

## Exp C — Cross-level refinement (`lcp_crosslevel`)

**Method.** For a (low, high) level pair on one genome, anchor each low core by its midpoint and, for every high core span, count the low-core midpoints inside it and hash the multiset of their labels. Report mean low-cores per high span and the fraction of high *labels* whose covered low-label signature is identical across every occurrence (refinement determinism). Cores overlap heavily, so midpoint-anchoring is used rather than strict containment.

**Result (whole genome):**

| genome | pair | mean children/parent | coverage | refinement-deterministic |
|---|---|---:|---:|---:|
| ref | L4→L5 | 5.40 | 100 % | 42,054,394 / 42,481,516 = **98.99 %** |
| ref | L5→L6 | 5.41 | 100 % | 18,431,078 / 18,543,934 = **99.39 %** |
| mixed | L5→L6 | 5.41 | 100 % | 18,437,495 / 18,550,351 = **99.39 %** |

**Reading.** Each higher-level core spans ~5.4 lower-level cores (consistent across levels), and **~99 % of high-level labels always cover the same low-level label signature** wherever they occur — the hierarchy is a near-deterministic refinement. The ~1 % that vary are precisely the colliding/repeat labels (cf. Exp B's ~1 % collision rate and Exp 1's 99.0 % faithful labels): a label that conflates distinct sequences necessarily covers different children in different places. So the level hierarchy is as consistent as the labels themselves are unique.

---

## Exp J — Local-consistency / windowing invariance (`lcp_consistency`)

**Method.** Is a core computed inside a local window identical to the core computed from the whole record? For each record ≥3 Mbp, parse it whole (`init_lps_offset`+`lps_deepen`, what svp uses) into a map `{(start,end)→label}`, then re-parse it in overlapping 20 Mbp windows (1 Mbp trust margin ≫ core length). For every core in a window's trusted interior — ≥1 Mbp from an *artificial* cut — check it matches the whole-record core exactly (position + label). (lcptools ships `init_lps4` for chunked split-merge, but it segfaults on the >100 Mbp chromosomes; rather than modify the library we do the windowing in-tool, testing the same locality property more directly.)

**Result (whole genome, L5):** ref **44,230,361 / 44,230,361 = 100.0000 %** identical across 32 records; mixed **44,338,334 / 44,338,334 = 100.0000 %** (2:50 each). The 32 records ≥3 Mbp hold ~97 % of all cores; the 423 small contigs fit in one window and are trivially consistent. **Reading.** A core is *fully determined by a bounded local neighbourhood* — windowing the genome anywhere (≥1 Mbp from a real core) reproduces the exact same cores, positions, and labels. svp's per-record parse therefore equals what a chunked/read-mode parse produces in any interior, so future FASTQ read-mode inherits the same anchors. (Synthetic cross-checks: 50 kb-chunk × 500 kb and 20 Mbp-window × 50 Mbp both 100.0000 %.)

---

---

# Part 2 — Re-keying on the ±10 bp flank signature, and Huffman in depth

Part 1 keys bags on the raw 32-bit lcptools `label`. svp's *actual* anchor identity is the **signature** `sig(core) = murmur64(seq[start−10 … end+10], seed = label)` — a 64-bit hash of the ±10 bp window (CORE_ANALYSIS `distHL`). Part 1's Exp B showed the label collides ~414k times at L5, half of it pure uint32 birthday accident. Here every experiment is repeated with `--sig` (the flank-signature key) to measure what the production signature actually buys. **Method note:** the LCP *partition* is unchanged — `parse1` never sees the key — so partition-level metrics (core counts, the rotation orbit, RC mirror-match, windowing-invariance) are key-invariant by construction; what the signature changes is the **collision / fibre structure** built on top.

## P2.1 — Bags & collision anatomy under the signature (Exp 1, B re-keyed)

Headline contrast (ref, L5):

| key | distinct keys | collisions (≥2 distinct seqs) | max bag size |
|---|---:|---:|---:|
| raw label (uint32, core-only) | 42,479,424 | **414,831 (0.977 %)** | 683 |
| flank signature (murmur64, ±10 bp) | 43,163,867 | **0 (0.000 %)** | 1 |

**The ±10 signature has zero lossy collisions** — every distinct signature maps to exactly one core sequence (distinct-sigs = distinct-seqs = 43,163,867 exactly). Both collision classes from Exp B — the 48.7 % structural *and* the 51.3 % 32-bit-birthday — are eliminated: widening the digest from 32→64 bits and from core-only to ±10-window content makes it injective on this genome (CORE_ANALYSIS's "zero observable 64-bit murmur collisions").

This holds at **every level and on both genomes** — sig collisions are 0 for ref/mixed at L4/L5/L6 (vs the label's 4.41 M / 592 k / 157 k conflated sequences):

| level | label loss (seqs−labels) | sig collisions |
|---:|---:|---:|
| 4 | 4,410,743 | **0** |
| 5 |   592,590 | **0** |
| 6 |   156,935 | **0** |

And the collision-anatomy tool under `--sig` confirms it directly: **0 collision labels** on both ref and mixed at L5 (vs 414,831 / 414,253 under the raw label). *The signature is the fix for the label's lossiness* — the entire structural+birthday fibre structure of Exp B collapses to nothing.

## P2.2 — Rotation, cross-level, perturbation, consistency under the signature

**Rotation (Exp 2).** Core counts are **byte-identical** to the Part-1 label runs across all 24 orderings (ref+mixed) — confirming the partition is key-invariant — and **sig collisions are 0 for every ordering**. The signature is collision-free regardless of alphabet order.

**Cross-level refinement (Exp C).** Determinism *rises* under the content-faithful key, closing most of the label-collision gap:

| pair | determinism (label) | determinism (sig) |
|---|---:|---:|
| ref L4→L5 | 98.99 % | **99.80 %** |
| ref L5→L6 | 99.39 % | **99.80 %** |
| mixed L5→L6 | 99.39 % | **99.80 %** |

The residual ~0.2 % is a midpoint-binning edge effect (a low core whose midpoint sits on a high-core boundary can fall in/out across repeat copies), not key lossiness.

**Perturbation (Exp E).** The signature perturbation is **wider** than the boundary perturbation, exactly because a core's signature depends on its ±10 bp window — a SNP changes the signature of neighbours up to ~10 bp away even when their boundaries are stable:

| level/type | cores changed (label / sig) | mean radius (label / sig) |
|---|---:|---:|
| L4 SNP | 2.27 / **2.93** | 57.2 / **66.4** bp |
| L5 SNP | 2.29 / **2.56** | 134.7 / **143.5** bp |
| L6 SNP | 2.31 / **2.42** | 315.3 / **323.9** bp |
| L5 INS | 2.27 / 2.53 | 133.7 / 142.7 bp |
| L5 DEL | 2.29 / 2.55 | 134.5 / 143.3 bp |

The radius widens by ~9 bp (≈ the flank) and ~+0.3–0.6 extra cores change — still tightly bounded; svp's distant anchors remain untouched.

**Consistency (Exp J).** Still **exactly 100.0000 %** (ref 44,230,361 / mixed 44,338,334 trusted-interior cores) — the signature is windowing-invariant just like the partition.

## P2.3 — Huffman, in depth (Exp 3 expanded)

**Method.** Three angles on the "Huffman-like first layer": (a) the genuine information-theoretic properties of DNA k-mers — Shannon entropy, true Huffman avg/max code length (built from a min-heap over k-mer counts), and frequency skew, for k=2,3; (b) the **alphabet-ordering sweep** `--korder freq|lex|revfreq|rand` — since `parse1` is order-only, the choice of which integer each k-mer gets is a 64-symbol rotation, and frequency-rank (Huffman-like) is just one ordering; (c) the triplet bags re-keyed on the ±10 bp-window signature.

**Result (a) — how compressible are DNA k-mers? (ref, whole genome):**

| k | k-mers | total | Shannon entropy | max | redundancy | Huffman avg | max code | freq skew |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 16 | 1,524,657,215 | 3.9094 bit | 4 | 0.0906 bit/sym | 3.9546 bit | 5 bit | 9.7× |
| 3 | 64 | 1,016,437,626 | 5.8317 bit | 6 | 0.1683 bit/sym | 5.8683 bit | 9 bit | 17.2× |

DNA k-mers are **mildly non-uniform**: triplets carry only ~0.17 bit/symbol of redundancy (entropy 5.83 of 6 bits), but the frequency skew is real — the most common triplet is **17.2×** the rarest, and the true Huffman code runs from short codes up to **9 bits**. So a Huffman first layer is information-theoretically *possible* but buys little compression; its value here is as an *ordering*, not a code.

**Result (b) — does the k-mer ordering change the partition? (ref, L5):**

| k | freq (Huffman) | lexicographic | random | reverse-freq | spread |
|---:|---:|---:|---:|---:|---:|
| 2 | 22,563,556 | 22,538,897 | 22,682,962 | 22,790,186 | 1.1 % |
| 3 | 14,750,882 | 14,770,185 | 14,808,323 | 14,813,824 | 0.43 % |

(symbol-core counts, sentinel-filtered, matching Exp 3 above)

The ordering of the 64 (or 16) k-mer symbols is a large rotation, and — as with the base-alphabet rotation (Exp 2) — it **does** shift the partition, but only mildly (0.4–1.1 % core-count spread). **Frequency-rank (Huffman) ordering sits at the low-core end** (fewest cores at k=3; near-lowest at k=2, just above lexicographic), and reverse-frequency (anti-Huffman) gives the most. So ordering k-mers by frequency is a weak attractor toward fewer/longer cores, but the effect is small — the partition is governed far more by *k* (the symbol granularity) than by the ordering within a given *k*.

**Result (c) — triplet bags under the signature:** with the ±10 bp-window signature, the triplet first layer is also **collision-free** — 0 collisions at k=2,3 on ref+mixed; the only multiplicity is genuine repeats (max repeat-frequency 458 at k=2, 136 at k=3).

---

## Synthesis

1. **The label map's fibres are half-accidental** (Exp 1, B). At L5, 99.0 % of labels are faithful sequence ids; of the 0.98 % that collide, **48.7 % are structural** (the digest conflates sequences sharing their boundary children — shared affix 65–256 bp, Jaccard ~0.88) and **51.3 % are pure 32-bit murmur birthday accidents** (shared nothing) — matching the uint32 birthday prediction. The label space, not hash quality, is the limit, which is exactly why svp signs cores with a 64-bit window murmur, not the 32-bit label.

2. **LCP parsing is chirally asymmetric** (Exp 2, A). All 24 alphabet orderings give *distinct* partitions (2.0 % core-count spread at L6; the S₄ action is essentially free), and forward cores have *no* reverse-complement symmetry (≤0.4 % mirror-match, falling with level). Both follow from one fact: `parse1`'s LMAX rule carries neighbour guards (`prev≤a`, `c≥next`) that LMIN lacks, so the scan is not invariant under order-reversal — and RC reduces to a reversal. svp parses one strand only; the other is unrelated.

3. **A coarser, frequency-ranked first layer improves label fidelity** (Exp 3). Per-base → dinucleotide → triplet drops the collision rate 0.977 % → 0.517 % → 0.379 % (more bits/symbol → richer labels) while shrinking the core set 45 M → 15 M. The k=3 level-1 distance-field overflow is real but isolated (243 N-block cores).

4. **LCP is blind to low-complexity and perfect periodicity** (Exp F). Homopolymers yield zero cores at every level; perfectly periodic tandems collapse to nothing by L5; a quasi-random tandem dup of unit U×K yields a *constant* ~4 labels each repeated exactly K times — the tandem-DUP uniqueness wall (`CORE_ANALYSIS` §B.10), reconstructed from first principles.

5. **Locality is tight and scales with level** (Exp E). A point edit (SNP/INS/DEL alike) perturbs ~2.3 cores regardless of level; the bp radius grows 57 → 135 → 315 bp from L4→L6, tracking median core length. Anchors a few hundred bp from a variant are untouched — the basis of svp's variant-robust chaining.

6. **The level hierarchy is a near-deterministic refinement** (Exp C), and **the parse is locally consistent** (Exp J). ~99 % of high-level labels always cover the same low-level label signature (the ~1 % that vary are the colliding labels of point 1); and a core computed in a local window is identical to the whole-record core — so future read-mode (reads = windows) inherits the same cores.

7. **SVs are LCP-local** (Exp 4). Mixed differs from ref by a constant +0.238 % cores at every level (chr19 SV content at genome-average density), with identical collision rates and fibre structure — the SV signal is real but spatially confined, exactly as `CORE_ANALYSIS` found.

8. **The ±10 bp flank signature is the cure for the label's lossiness** (Part 2). Re-keying every experiment on `murmur64(seq[start−10…end+10], seed=label)` — svp's production anchor — drives label collisions from 414,831 (L5) to **0**, at every level and on both genomes; distinct-sigs = distinct-sequences exactly. Cross-level determinism rises to 99.80 %, and consistency stays 100 %. The cost is a slightly *wider* perturbation footprint (~+9 bp, ~+0.3 cores), because a signature depends on its ±10 bp window — but anchors stay tightly local. Partition-level results (rotation core counts, RC, windowing) are unchanged: the key is layered on top of the partition, not part of it.

9. **Huffman buys ordering, not compression** (Part 2). DNA k-mers are only mildly redundant — triplet entropy 5.83 of 6 bits (~0.17 bit/symbol), though with a real 17.2× frequency skew and Huffman codes up to 9 bits. As an *ordering* of the k-mer alphabet, frequency-rank (Huffman) shifts the partition only 0.4–1.1 % and sits weakly at the fewer-cores end; the symbol granularity *k* dominates far more than the ordering within a *k*. The triplet layer is also collision-free under the signature.

*Raw per-run logs in `outputs/lcp_topology/`; machine-readable CSVs in `docs/analysis/lcp_*.csv`; collision dumps in `outputs/lcp_topology/anatomy_*.tsv` and `bags_*_L5_id.tsv`. Part 2 (`--sig`) runs prefixed `sig*`/`huff_*`.*
