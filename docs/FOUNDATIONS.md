# Foundations — LCP cores, signatures, and why they work

First-principles facts the `svp` pipeline rests on: what an LCP core is, why the graph can be built
from reads alone, and which limits are structural. Measured whole-genome (human, AMD 9950X) unless
noted; collapsed from the original LCP_TOPOLOGY / CORE_ANALYSIS / BIOLCP probes. Measurements taken
under a retired design are kept and labelled **[historical]** — they are why the design is what it is.

## What an LCP core is
A locally-consistent parse (lcptools, `biolcp` branch) recursively partitions a sequence into
**cores** — variable-length substrings whose boundaries depend only on a local window, so the
same locus is meant to yield the same core whether it is seen whole-record or inside a read
(windowing-invariance: measured 100% at the time — see *Open questions* for the caveat that
measurement has since been shown to have). svp runs at **level 5** (median core ~150 bp; L4 ~63 bp,
L6 ~352 bp).

## Three load-bearing properties
1. **Locality.** A point edit (SNP/INS/DEL) perturbs only ~2.3 cores over a ~135 bp radius at
   L5. Anchors a few hundred bp from a variant are untouched → variant-robust. SV signal is real
   but spatially confined (mixed vs ref differ by a constant +0.24% of cores). This is what makes a
   *coordinate-free* graph work: a somatic allele shows up as a handful of tumor-only cores wedged
   between untouched shared cores, and those shared neighbours are the anchors the VCF stage locates.
2. **Chiral asymmetry.** The parse is *not* reverse-complement symmetric: `parse1`'s LMAX neighbour
   guards have no LMIN mirror. svp therefore extracts **both** strands explicitly (fwd + RC) — each
   locus appears as two independent nodes `W` and `RC(W)`, with no per-emission canonicalization.
   (Mirror-match rate: ≤0.4% in the original probe; `tools/lcp_symmetry` / `tools/lcp_union` restate
   it as 0.027% of intervals mirrored. Both say the same thing: RC symmetry is essentially absent.)
3. **Adjacent cores overlap.** Consecutive L5 cores overlap ~85–105 bp on average (35 bp at L4,
   199 bp at L6), so the *core step* (bp advanced per core) is far below the core length. Three
   consequences that the current v3 design depends on:
   - a graph edge between adjacent cores is a **genuine sequence overlap**, not a zero-overlap bridge;
   - a path of *n* cores spans ≈ *n × step* bp, which is exactly how `tools/anchor_v3` converts the
     caller's **core counts** into bp (it measures `total_bp / total_cores` on the reference scan
     rather than assuming a constant);
   - the repeat-core blacklist merge (`src/blacklist.cpp`) can absorb a repeat run `X` into its
     flanks as two *overlapping* cores `(A+X)` and `(X+B)` — they still share real bytes, so
     `A→B` remains a legitimate overlap edge.

   **[historical]** In v2 this same overlap was what let `refine`/`call` splice node *bytes* by exact
   KMP suffix/prefix match (never a fixed-width shortcut); a zero-overlap seam was FATAL. v3 stores no
   sequence at all, so there is no splice left to protect — the overlap now only has to exist, not be
   recomputed.

## The signature
A node's identity is `core_sig` = a **64-bit hash of the core's own bytes — no flank** — computed once
in `src/extract.cpp` and stored as the node's *only* payload in the v3 graph (`u64 sig[K]`). The
current mixer is **MurmurHash64A** (seed 0). 64-bit collisions are unobserved across whole-genome
configs (distinct sigs == distinct sequences), so `sig ↔ node id` is a bijection and every downstream
tool re-derives node identity by re-extracting cores and hashing them the same way.

**The mixer is part of the on-disk contract.** Earlier revisions hashed the same bytes with XXH3-64;
some comments still say so. Any sig-keyed artifact — `graph.bin` / `refined.bin`, `*.altsigs.bin`, the
repeat blacklist — is only valid against the hash that produced it. Change the hash and everything
must be rebuilt together, or the blacklist silently matches nothing and anchoring silently fails.

> **[historical]** The identity was originally `core ± 10 bp flank`, murmur-hashed, to make it
> position-discriminating. Measured then: dropping the flank produced ~58k zero-overlap splice seams
> and cut ≥50 bp recall 0.495 → 0.385, and using lcptools' 32-bit boundary `label` instead was far
> worse (0.045, 193k seams). Those failures were failures of the *byte-splice*, which no longer
> exists; the core-only identity survives because adjacent cores overlap natively (property 3).

## Ceilings (structural, not fixable by hashing)
- **Genome repetition caps uniqueness:** 8.07 / 5.15 / 3.83% of cores at L4/L5/L6 have a
  byte-identical duplicate somewhere (flank=0). Deeper level shrinks it; the production answer is the
  repeat blacklist (`src/blacklist.cpp`), which merges recurrent cores into their flanks instead of
  leaving them as ambiguous hubs. `tools/uniquify_levels` measures how many stacked merge levels it
  would take to make *every* core of an assembly unique.
- **LCP is blind to low-complexity / perfect periodicity:** homopolymers yield zero cores;
  a tandem dup of unit U × K collapses to ~4 labels each repeated K times — the **tandem-DUP
  wall** (~77% uniqueness on mixed-SV inputs, uncrossable even with context rescue R=8).
- **[historical]** At the flank-era production setting (L5, flank=10) the graph was near tree-like:
  ~5.5 cycles per 1000 nodes, the cycles being genome repeats. Not re-measured for the flank-free v3
  graph — `tools/census_v3` is the tool for that (degrees, SCCs, components, unitig N50).

## Consequence for somatic calling (dual fwd+RC)
Pooling both strands lets `min_mult` go high and prune ~90% of error sigs, but it **doubles**
node count and makes inversions invisible as somatic: the inverted body becomes colour-3
(shared with normal's RC pass); only the breakpoint-junction cores are tumor-only. Folding
strand pairs is therefore a real lever — but it cannot be done by canonicalizing each emission
(`min(bytes, revcomp(bytes))`), because the fwd and RC parses do not cut at mirrored positions
(property 2): the two strands' cores are different substrings, not reverse complements of each other.

## Open questions (active probes in `tools/`)
Not yet settled; the tools exist, the conclusions do not belong in this doc until measured.
- `lcp_symmetry` — quantifies how un-mirrored the fwd/RC parses are. Gates any canonical-core /
  bidirected redesign.
- `lcp_union` — a wrapper making the parse RC-symmetric *by construction* (cut set
  `C(S) = cuts(S) ∪ mirror(cuts(revcomp(S)))`), and what that costs in core count/length.
- `lcp_consistency` — **the gate**: do two *different* reads covering the same locus emit the same
  cores? Symmetry is necessary but not sufficient; a graph only exists if the parse is locally
  consistent across read framings.
- `lcp_warmup` — shows the earlier consistency measurement had a blind spot: it only compared cores a
  window *produced*, so a core the window silently *failed* to produce scored as nothing. Plain LCP can
  look 100% windowing-consistent while still losing cores near a sequence boundary. Read the
  "measured 100%" above with that caveat (svp's `--margin 10` drops read-end cores, which mitigates
  but does not by itself prove the point).
