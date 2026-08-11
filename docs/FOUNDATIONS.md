# Foundations — LCP cores, signatures, and the limits that come with them

The first-principles facts the pipeline rests on: what a core is, why a graph can be built from reads
alone, and which limits are structural rather than a matter of tuning. Measured whole-genome (human,
AMD 9950X) unless noted.

## What an LCP core is

A locally-consistent parse (lcptools, branch `biolcp`) recursively partitions a sequence into
**cores** — variable-length substrings whose boundaries depend only on a local window, so the same
locus tends to yield the same core whether it is seen inside one read or another. snape runs at
**level 5**: median core ~217 bp raw, ~150 homopolymer-compressed bases (L4 ~63, L6 ~352).

## Three load-bearing properties

**1. Locality.** A point edit perturbs only ~2.3 cores over a ~135 bp radius at L5. Anchors a few
hundred bp from a variant are untouched. This is what makes a coordinate-free graph work: a somatic
allele shows up as a handful of tumour-only cores wedged between untouched shared cores, and those
shared neighbours are the anchors everything downstream hangs off — the flanks of a junction contig,
and the cores `pick` matches normal reads on.

**2. Chiral asymmetry.** The parse is *not* reverse-complement symmetric — `parse1`'s LMAX neighbour
guards have no LMIN mirror. snape therefore extracts **both** strands explicitly, and each locus
appears as two independent nodes `W` and `RC(W)` with different signatures and no core in common.
This cannot be fixed by canonicalising each emission (`min(bytes, revcomp(bytes))`), because the two
strands' cores are different substrings, not reverse complements of each other. Two consequences:
the graph is about twice the size it would otherwise be, and a large inversion's **interior** is
shared rather than somatic — the inverted sequence lands on the normal reads' mirror-family nodes.
What is callable about an inversion is its two breakpoint junctions.

**3. Adjacent cores overlap.** Consecutive L5 cores overlap ~85–105 bp on average (35 bp at L4, 199
at L6), so the *core step* — bases advanced per core — is far below the core length. Hence a graph
edge is a genuine sequence overlap and not a zero-overlap bridge, and hence `junctions` can cut a
contig straight out of one read's bytes: the span between two of that read's core boundaries IS the
sequence, with no assembly and no cross-read merge.

## The signature

A node's identity is `core_sig` = **MurmurHash64A (seed 0) over the core's homopolymer-compressed
bytes**, computed once in `src/extract.cpp` and stored as the node's only payload. 64-bit collisions
are unobserved across whole-genome configs, so `sig ↔ node id` is a bijection and any stage can
re-derive node identity by re-parsing sequence and hashing it the same way.

**The mixer and the compression are both part of the on-disk contract.** Change either and every
sig-keyed artifact — `graph.bin`, `drop.sigs` — must be rebuilt together, or they silently match
nothing.

Compression is measured, not assumed: on hcc1954's own normal reads 83.4% of read errors are indels
and those are 2.3× enriched in homopolymer runs ≥ 3, while substitution errors show no enrichment at
all (`scripts/hifi_error_profile.py`, numbers in `RESULTS.md`). Compressing runs before hashing
neutralises about two thirds of what could otherwise change a locus's identity. Core start/end stay
in raw coordinates, so nothing downstream has to think in compressed space.

## Structural ceilings

- **Genome repetition caps uniqueness.** 8.07 / 5.15 / 3.83% of cores at L4/L5/L6 have a
  byte-identical duplicate somewhere. Repeats are ordinary nodes here — there is no reference list to
  mask them with — so what separates a repeat from an event is coverage and, at the calling end, the
  tandem-repeat fraction of the contig's own bases.
- **LCP is blind to perfect periodicity.** A homopolymer yields no cores at all, and a tandem
  duplication of unit U × K collapses to ~4 labels each repeated K times. Events whose entire
  signature is periodic are not visible to this graph, whatever the thresholds.
- **A core is ~217 bp, so a core signature cannot see a 50 bp change.** Two cores match or they do
  not. That is why `bubbles` opens the arms up with minimizers rather than comparing signatures.
