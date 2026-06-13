# Foundations — LCP cores, signatures, and why they work

First-principles facts the `svp` pipeline rests on. Measured whole-genome (human, AMD 9950X)
unless noted; collapsed from the original LCP_TOPOLOGY / CORE_ANALYSIS / BIOLCP probes.

## What an LCP core is
A locally-consistent parse (lcptools, `biolcp` branch) recursively partitions a sequence into
**cores** — variable-length substrings whose boundaries depend only on a local window, so the
same locus yields the same core whether seen whole-record or inside a read (windowing-invariant,
measured 100%). svp runs at **level 5** (median core ~150 bp; L4 ~63 bp, L6 ~352 bp).

## Three load-bearing properties
1. **Locality.** A point edit (SNP/INS/DEL) perturbs only ~2.3 cores over a ~135 bp radius at
   L5. Anchors a few hundred bp from a variant are untouched -> variant-robust. SV signal is real
   but spatially confined (mixed vs ref differ by a constant +0.24% of cores).
2. **Chiral asymmetry.** The parse is *not* reverse-complement symmetric (<=0.4% mirror match):
   `parse1`'s LMAX neighbour guards have no LMIN mirror. svp therefore extracts **both** strands
   explicitly (fw + RC) -- each locus appears as two nodes `W` and `RC(W)`.
3. **Adjacent cores overlap.** Consecutive L5 cores overlap ~85-105 bp on average (35 bp at L4,
   199 bp at L6). This overlap lets Mode A/B splice node bytes by exact suffix/prefix (KMP)
   match -- never a fixed-width shortcut.

## The signature (the central design choice)
A node's identity is `sig = murmur64(core +/- 10 bp flank, seed 0)`, **not** the core alone and
**not** lcptools' 32-bit boundary `label`.
- **Murmur is collision-free in practice.** Zero observable 64-bit collisions across 33
  whole-genome configs (distinct sigs == distinct sequences exactly).
- **The flank is load-bearing -- it makes the sig position-discriminating.** Two loci with the
  same core but different neighbours get different sigs, so each node ~= one locus and its stored
  +/-flank window bytes are *consistent* across reads. This is what lets the byte-splice in
  Mode A/B hold across a `min_mult`-pruned coordinate gap (the +/-10 flank bridges the gap).
  Measured: dropping the flank produces 58k zero-overlap splice seams and craters recall
  (>=50bp recall 0.495 -> 0.385); the 32-bit `label` is worse (0.045, 193k seams).
- **The `label` seed adds nothing** for flank >= 0 (proven L4/L5/L6 WG); it only carries info for
  negative flank (window subset of core).

## Ceilings (structural, not fixable by hashing)
- **Genome repetition caps uniqueness:** 8.07 / 5.15 / 3.83% of cores at L4/L5/L6 have a
  byte-identical duplicate somewhere (flank=0). Deeper level / wider flank shrinks it.
- **LCP is blind to low-complexity / perfect periodicity:** homopolymers yield zero cores;
  a tandem dup of unit U x K collapses to ~4 labels each repeated K times -- the **tandem-DUP
  wall** (~77% uniqueness on mixed-SV inputs, uncrossable even with context rescue R=8).
- At production (L5, flank=10) the graph is near tree-like: ~5.5 cycles per 1000 nodes; the
  cycles are genome repeats.

## Consequence for somatic calling (dual fw+RC)
Pooling both strands lets `min_mult` go high and prune ~90% of error sigs, but it **doubles**
node count and makes inversions invisible as somatic: the inverted body becomes colour-3
(shared with normal's RC pass); only the breakpoint-junction cores are tumor-only. Folding
strand pairs is a known future step.
