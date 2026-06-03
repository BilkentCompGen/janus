# Graph-inspector — first-run diagnostic findings

Compiled 2026-06-03 from the first end-to-end run of `inspect/` (Stage 1 census
+ Stage 2 inspect + Stage 3 viz) on `outputs/svp/wgs/` (112,977 bubbles from
30× whole-genome HiFi normal/tumor; truth = `simulation/truth/wgs/out.vcf`).

This memo answers the seven diagnostic questions from the plan's verification
section and lists concrete next steps for the caller.

## TL;DR

1. **The graph is fine** — 99.62% sig-recall holds; 49.55% D-bubbles + 50.45%
   I-bubbles + 0 C-bubbles, exactly matching `docs/CLAUDE.md`.
2. **The FP problem is uniform, not localized.** FP rate is 95–98% across
   EVERY region class (HLA, centromere, segdup, simple repeat, gap, AND clean
   no-flag regions). The repeats aren't the culprit; the noise is everywhere.
3. **Singletons dominate.** 63.82% of bubbles have ≤2 LCP nodes. Among FP
   I-bubbles, 52% are singletons and 49% are 3-5 nodes. Almost no FP I-bubbles
   are >5 nodes.
4. **Mode B v1 only handles single-N + single-T paths after Mode A
   compaction.** Multi-path interior (>1 N or >1 T super-node) is explicitly
   skipped. This blocks the most interesting cases:
   - DUP and INV bubbles both decompose into multi-path interior (because the
     duplicated/inverted content forms parallel paths).
   - sv17947 DUP 740: 0 N + 2 T → skipped.
   - sv17978 INV 267: 6 sub-bubbles, 5 of 6 skipped (multi-path).

## Answers to the 7 diagnostic queries

### 1. What fraction of FP bubbles are in HLA / centromere / segdup regions?

From `bubble_summary.md`:

| region | bubbles | FP rate |
|---|---|---|
| HLA | 106 | 96.23% |
| centromere | 2,533 | 96.68% |
| segdup | 5,792 | 95.55% |
| simple repeat | 6,852 | 90.67% |
| gap | 3,162 | 98.32% |
| **clean (no flags)** | **68,048** | **97.33%** |

**Surprising:** Clean regions have a HIGHER FP rate (97.33%) than simple
repeats (90.67%). Repetitive regions actually have a slightly better FP rate
because simple-repeat-overlapping bubbles tend to be anchorable larger
structures, not single-node noise.

### 2. Size distribution of TP vs FP bubbles?

Whole-population from census:

| n_nodes | total | TPs | FPs |
|---|---|---|---|
| 1 | 35,293 | very few | most |
| 2 | 36,809 | few | most |
| 3-5 | 40,588 | mixed | most |
| 6-10 | 107 | mixed | mixed |
| 11-50 | 106 | enriched | few |
| 51-500 | 48 | enriched | few |
| 501+ | 26 | mostly | few |

FP I-bubbles specifically: 20,687 singletons + 19,340 small (3-5) + 15 mid +
11 large. **The FPs cluster strongly in small sizes.** The few large bubbles
are mostly real SV-region tangles.

### 3. Are FPs concentrated on specific chroms?

No striking concentration. From top-20 by FP rate, every main chromosome
sits in the 96-98% FP band:

```
chr1: 98.41%  chr2: 98.17%  chrX: 98.27%
chr10: 97.93%  chr14: 97.74%  chr18: 97.47%
chr3: 97.38%  chr11: 97.65%  chr7: 97.31%
…
```

The alt-contig haplotypes (X_KI270880v1_alt, 17_GL000258v2_alt, etc.) are 100%
FP but only contribute 12-16 bubbles each.

### 4. What % of FP I-bubbles are size 1-2 (RC-pair singletons)?

52% (20,687 / 40,053). The companion 49% are 3-5 nodes — also small. The
"≥1 sig is reverse-complement-pair of another sig" hypothesis from
`docs/CLAUDE.md` predicts this directly: a single LCP-core appears as TWO
graph nodes (forward sig + RC sig), each color-2 (tumor-only) because the
normal coverage at that exact strand orientation drops below `min_mult`.
This creates "ghost" 2-node I-bubbles that have no biological reality.

**Implication:** strand-canonicalization of murmur sigs (already documented
as upstream work in `docs/PERFORMANCE.md`) would roughly halve the FP I-bubble
count, lifting bubble-precision from 63% toward 78% (matching sig-precision).

### 5. For random TP bubbles per svtype, does modeB correctly tag svtype+svlen?

Tested on the 4 named truth SVs from the plan:

| truth SV | bubbles matched | mode B classification |
|---|---|---|
| **sv17602** (DEL 12 @ chr1:12.6M) | 2 | bubble 11295: SV=D dN=449; bubble 57107: SV=I dT=49 |
| **sv17268** (INS 49 @ chr1:3.5M) | 2 | bubble 3004: SV=D dN=349; bubble 40616: SKIP (multi-path) |
| **sv17947** (DUP 740 @ chr1:117M) | 1 | bubble 49703: SKIP (0-N + 2-T multi-path) |
| **sv17978** (INV 267 @ chr11:102M) | 6 | 5 of 6 SKIP; bubble 102211: SV=C dN=316 dT=626 |

None of the four cleanly emit the truth svlen. Reasons:

- **sv17602**: tandem-repeat copy-number contraction — actual 12-bp DEL is
  hidden inside a 3-copy → 2-copy event represented as 449-vs-49 bp in the
  graph. modeB sees the full copy-number difference, not the 12 bp.
- **sv17268**: bubble 40616 has 2 T-path super-nodes after Mode A — the INS
  is split into multiple unitigs at internal branch points.
- **sv17947**: DUP has zero N-path and two T-path super-nodes (the
  duplicated copy + its junction adapter). v1 modeB can't handle this.
- **sv17978**: the INV decomposes into 6 sub-bubbles per
  `project_caller_v1_fix_attempts` memory. The one classified bubble has
  asymmetric N vs T lengths because the RC trim didn't quite match.

**Conclusion:** Mode B v1 handles only the cleanest case (one N-path + one
T-path super-node after Mode A). The interesting SVs (DUP/INV/tandem-repeat)
all need multi-path support to call correctly. The dN/dT numbers from v1
are diagnostic (they tell you the size of the structural difference) but
they overestimate the actual SV size when the bubble represents copy-number
variation.

### 6. For random FP I-bubbles, what does modeB call?

Not run as a batch yet (loading time ~30 s per invocation makes 20 separate
calls expensive). Predicted from the size distribution: 52% of FP I-bubbles
are singletons (1 LCP node). After Mode B's anchor-prefix-suffix absorption,
those collapse to tiny diffs (1-5 bp) — likely `SV:Z:S` (SNP) or `SV:Z:EQ`
(noise). This would be the strand-pair singleton signature.

Action item: add a batch mode to `inspect` (load data once, process a list
of bubble IDs) and run on 100 random FP I-bubbles to confirm.

### 7. INV vs DUP topology — distinguishing signature?

Yes, even at the topology level (without modeB classification):

- **DUP** (sv17947 in `inspect/sv_sv17947/`):
  - 2 SHARED anchors, 0 NORMAL super-nodes, 2 TUMOR super-nodes
  - The tumor side has the extra copy; normal side has no "filler" path
  - Both T-super-nodes connect to both anchors (parallel topology)
  - Visual signature: anchors on one side, ONLY tumor on the other

- **INV** (sv17978 in `inspect/sv_sv17978/`):
  - 6 sub-bubbles spread across the inversion span
  - Mix of D and I types within the 6
  - `mixed_strand=true` (anchors disagree on strand) is the canonical INV flag
  - Visual signature: scattered sub-bubbles + strand-mismatch flag

- **DEL** (sv17602 in `inspect/sv_sv17602/`):
  - 2 SHARED anchors, 1 NORMAL super-node (the deleted content), 1 TUMOR
    super-node (the alternate at the deletion junction)
  - Clean bipartite topology

- **INS** (sv17268 in `inspect/sv_sv17268/`):
  - 2 SHARED anchors, mix of N + T super-nodes
  - The inserted bases live in the T super-node(s)

These topology signatures suggest a much simpler INV/DUP detector than
modeB: classify a bubble by counting (n_anchors, n_N_supernodes,
n_T_supernodes) after Mode A compaction.

| signature | n_anchors | n_N | n_T | type |
|---|---|---|---|---|
| simple DEL/INS | 2 | 1 | 1 | call modeB |
| pure DEL | 2 | 1 | 0 | call as DEL of n_N.LN |
| pure INS | 2 | 0 | 1 | call as INS of n_T.LN |
| DUP candidate | 2 | 0 | 2+ | call as DUP of sum(T.LN) - (anchor overlap) |
| INV candidate | 2+ | mixed | mixed | check strand_mismatch; if true, RC + retry |
| BND/translocation | 2+ different chroms | any | any | emit BND |

This is the v2 design for Mode B.

## Mode A's "no_overlap" warnings — root cause

Mode A's exact suffix/prefix overlap algorithm fails in 4 of the 4 SV
inspections (16 warnings total). Hypothesis: the suffix/prefix overlap
direction depends on read orientation, but my walk picks an arbitrary chain
order. If the chain order is A → B but the read content has B → A, then
`A.suffix vs B.prefix` returns 0 and the merge fails.

Fix: try both directions (`A.suffix vs B.prefix` AND `B.suffix vs A.prefix`)
and accept whichever returns positive overlap. If both return positive,
the larger one wins. If both return 0, fall through to the no-merge path
(per the LCP-feedback rule: don't fake a splice).

This is a v2 fix for `inspect.cpp`. The current v1 still produces correct
modeB CLASSIFICATION (the trim sees the same L-prefix and R-suffix even
without merging), but `bandage-ng` rendering of modeA.gfa shows the chain
broken into more nodes than it should be.

## Next steps for the caller (separate plan)

The diagnostic confirms what `project_caller_v1_fix_attempts.md` already
suspected: the recall loss isn't in bubble topology, it's in:

1. **DUP detection** needs the 0-N + 2-T super-node signature in Mode A,
   not anchor-overlap geometry (which undersizes per the earlier failed
   v1 attempt).
2. **INV detection** needs `mixed_strand=true` flag + RC of the T-allele
   in Mode B — already partially implemented; the asymmetric dN/dT
   suggests the RC isn't perfect.
3. **Tandem repeats** (sv17602) need copy-number counting, not bubble
   trim. Beyond Mode B's scope; defer to a future "repeat-aware caller"
   pass.
4. **Multi-path Mode B**: walk all N and T paths in the bubble, build
   candidate alleles per path, apply trim, accept the one with the
   smallest diff size.

A new plan will be written before any caller iteration begins. The
inspector remains the truth source for "what does the graph actually
look like here?"
