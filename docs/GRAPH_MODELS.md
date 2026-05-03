# Graph Models for Genomics: De Bruijn Graphs and Variation Graphs

A reference document on the two graph families used to represent sequence data and pangenomes. Pure conceptual material — what the structures are, where they came from, what variants exist, and how they relate. No project-specific content.

---

## 1. De Bruijn Graphs (DBGs)

### Origin

Pure combinatorics, not biology. The mathematical de Bruijn graph (Nicolaas Govert de Bruijn, 1946) is a graph where the vertices are all length-n strings over an alphabet of m symbols, and there is a directed edge from `u` to `v` whenever you can produce `v` by dropping `u`'s first symbol and appending one new symbol. Eulerian circuits on it generate de Bruijn sequences (every length-n string appears exactly once).

Brought into bioinformatics by **Idury & Waterman (1995)** and made central to assembly by **Pevzner, Tang & Waterman (2001)** in their EULER assembler. The key insight: chop reads into k-mers, and the assembly problem becomes finding an Eulerian path in a DBG instead of a Hamiltonian path in an overlap graph (NP-hard in general → polynomial for the Eulerian formulation).

### Formal definition (bioinformatics flavor)

Two equivalent definitions are in circulation:

- **Node-centric:** nodes are k-mers, an edge `u→v` exists iff `suffix(u, k−1) = prefix(v, k−1)`.
- **Edge-centric:** nodes are (k−1)-mers, edges are k-mers (each k-mer is an edge from its (k−1)-prefix to its (k−1)-suffix).

Both give the same graph up to a node↔edge swap. The k-mer is the atomic unit; k is a fixed external parameter.

### Properties

- Node label = a fixed-length string. All nodes have the same label length.
- Two identical k-mers in the input collapse to a single node — the graph deduplicates by label.
- Edges are derived from sequence content alone — purely a function of the input set and k.
- **Deterministic:** same reads + same k → byte-identical graph from any tool.
- **Repeats ≥ k bp create cycles** automatically. This is unavoidable structural noise — the graph cannot distinguish a true tandem repeat from incidental k-mer reuse.
- No coordinate system, no privileged sequence, no haplotype memory. Every read that contributes a k-mer becomes anonymous once chopped.
- Eulerian-path formulation makes assembly polynomial-time (in the noise-free idealization; the practical NP-hardness comes from errors, repeats, and read constraints).

### The DBG family — variants you will encounter

| Variant | What it adds | Representative tools |
|---|---|---|
| **Classic DBG** | The basic node = k-mer construction | EULER, original Velvet |
| **Compacted DBG (cDBG)** | Collapses every maximal non-branching path (a "unitig") into one node. Same graph topology but far fewer nodes — orders of magnitude in practice | BCALM2, Bifrost, Cuttlefish 2, GGCAT |
| **Bidirected DBG** | Treats a k-mer and its reverse complement as the same node, with sided edges. Necessary for double-stranded DNA | Velvet (effectively), SPAdes, ABySS |
| **Colored DBG** | Each node/edge carries a bit-vector indicating which input samples (colors) contained that k-mer. Lets one graph represent many genomes | Cortex (Iqbal et al. 2012), Mantis, Bifrost, Rainbowfish, Bloom Filter Trie |
| **Weighted DBG** | Each edge carries a count = how many times the k-mer was observed. Used for error correction and coverage-based filtering | Most assemblers internally |
| **Succinct DBG (BOSS)** | Bowe, Onodera, Sadakane, Shibuya 2012 — represents the graph in `4m + o(m)` bits via an extension of the FM-index. Squeezed a human-genome DBG from ~300 GB to ~2.5 GB | TwoPaCo, BOSS-derived libraries |
| **Dynamic DBG** | Supports add/delete of k-mers without rebuilding | DynamicBOSS, Bifrost |
| **Variable-order DBG (VODBG)** | Chari & Mäkinen — lets you query the graph as if it were built at any k′ ≤ k_max, dynamically. Trades fixed-k brittleness for index complexity | Variable-Order BOSS |
| **Probabilistic DBG** | Backed by a Bloom filter — saves memory at the cost of false-positive edges. Designed for petabyte-scale metagenomics | khmer (PNAS 2012) |
| **String graph (related but distinct)** | Read-overlap graph after transitive-edge reduction (Myers 2005). Not a DBG but the structural cousin in the OLC paradigm | SGA, Falcon (long-read OLC) |
| **A-Bruijn graph** | "Adjacency" generalization — nodes are equivalence classes from a multiple alignment instead of fixed k-mers. The conceptual stepping stone between DBG and variation graph | Pevzner et al., used in some assemblers internally |

### Where DBGs are actually used

- **Short-read de novo assembly:** Velvet, ABySS, SOAPdenovo, ALLPATHS-LG, SPAdes — short-read assemblers are essentially universally DBG-based because reads are too short for OLC to work cheaply.
- **Variant calling without a reference / across many samples:** Cortex, DiscoSnp, McCortex (colored DBG → bubble calling).
- **RNA-seq assembly and splice graph reconstruction:** Trinity, Oases, SOAPdenovo-Trans — splice graphs come out as connected components of an RNA DBG.
- **Metagenomics:** MEGAHIT, metaSPAdes, IDBA-UD — colored or weighted DBGs scale to mixed-species samples.
- **k-mer indexing of huge sequence collections:** Mantis, Bifrost, Themisto, GGCAT — search "is this k-mer in any of these 100k bacterial genomes."
- **Pangenomics (DBG flavor):** Bifrost and mdbg build colored compacted DBGs across many genomes as a pangenome model.

---

## 2. Variation Graphs (VGs)

### Origin

A different lineage. Roots in **A-Bruijn / partial-order alignment graphs (Lee et al. 2002)** and **Cactus graphs (Paten et al. 2011)**, formalized for population-scale genomics by **Erik Garrison's vg toolkit (Garrison et al., Nat Biotech 2018)**. The motivation was the *reference bias* problem: one linear reference genome systematically misrepresents alleles that aren't in it. A graph that embeds known variation natively avoids this.

### Formal definition (Garrison 2018)

A variation graph is a triple

```
G = (N, E, P)
```

where

- **N** = a set of nodes; each node `n` carries a label `seq(n)` over `{A, C, G, T}` — a DNA segment of *arbitrary* length.
- **E** = a set of edges in a **bidirected** graph: each node has a "left" side and a "right" side, and edges connect node-side pairs. Four edge orientations exist (LL, LR, RL, RR); reverse-strand traversal is encoded by entering through the "wrong" side.
- **P** = a set of named **paths**, each path being an ordered list of `(node_id, orientation)` pairs. Paths embed concrete sequences (haplotypes, genomes, contigs, reads) into the graph: concatenate the labels along the path, taking the reverse complement on reverse-orientation visits, and you recover the original linear sequence.

A **valid traversal** must enter and exit each node through opposite sides (you read it forward or you read it backward — you don't enter and leave the same side).

### Key conceptual features

- Node labels are variable length. There is no "k."
- The graph is **bidirected**, so DNA's double-stranded nature is intrinsic. A node read forward gives one sequence; the same node read backward gives its reverse complement — no second node needed.
- **Paths are first-class objects**, stored alongside the graph. A path can revisit a node any number of times. Tandem expansions are encoded by repetition in the path list.
- Cycles are allowed but rare in well-formed VGs — they appear only for genuine repetitive structure (tandem repeats, true cycles in the population). Most regions are locally a DAG with the reference path threading through.
- **Tool-dependent topology**: unlike a DBG, the same set of haplotypes can yield different VGs from different builders. The path sequences match across tools; the graph shape doesn't.
- A variation graph subsumes a "sequence graph" (which is just `(N, E)` without paths). Adding the path set is what makes it a *variation* graph.

### The VG family — how different builders interpret the model

| Builder | How it builds | What it captures | Notes |
|---|---|---|---|
| **vg construct** | Reference + VCF → graph by splicing alts onto the ref backbone | All variants in the VCF | Original Garrison toolkit; reference path is privileged |
| **Minigraph** | Iterative reference-anchored augmentation by `minigraph -cxggs` | **SVs only** (≥ ~50 bp by default); SNPs/small indels are *not* added to topology | Compact, fast; lossy at the small-variant level |
| **Minigraph-Cactus** | Minigraph builds the SV skeleton, then Cactus realigns assemblies onto it to add base-level detail | SVs + small variants | Two-stage; produces both an SV graph and a base-level graph |
| **PGGB (Pangenome Graph Builder)** | All-pairs alignment via wfmash + seqwish → smoothxg → odgi | All variation, no reference bias | Symmetric, slow but unbiased; tandem dups become cycles instead of inserted alts |
| **Cactus (multi-genome alignment)** | Hierarchical alignment via Progressive Cactus | Originally for evolutionary multi-species alignment, retargeted for pangenomes | Predates vg; produces HAL graphs that can be converted to GFA |
| **Bifrost / mdbg** | Colored compacted DBG written out as GFA | Technically pangenome graphs but built on DBG semantics — fixed k under the hood | Sit on the boundary between the two families |
| **Seqwish** | Pairwise alignments → induced sequence graph | Building block used inside PGGB | Lossless w.r.t. its input alignments |

### Storage / serialization

The shared interchange format is **GFA (Graphical Fragment Assembly)**:

- `S` — segment lines (= nodes), carry the DNA sequence
- `L` — link lines (= edges), carry orientation and overlap length (overlap is usually 0 in pangenome use)
- `P` — path lines (per-individual walks; the original GFA1 path representation)
- `W` — walk lines (added in GFA 1.1, more compact for many haplotypes; suggested by Heng Li for pangenomes)

Indexing structures sit on top of GFA: **GBWT** (Sirén) for haplotype-aware indexing, **xg** (vg's static index), **odgi**'s succinct in-memory format. **Sequence-to-graph alignment** is hard in general (NP-hard for arbitrary VGs because the number of paths can be exponential in the number of variants), so production aligners — `vg giraffe`, `GraphAligner`, `minichain`, `vg map` — rely on haplotype indexes (GBWT) or restrict to specific graph classes.

### Where VGs are actually used

- **Reference-bias-free read mapping:** vg giraffe, GraphAligner, minigraph, PanAligner.
- **SV genotyping in populations:** vg call, Pangenie, BayesTyper — using a VG as the reference and reporting which embedded paths match the sample.
- **Human Pangenome Reference Consortium (HPRC):** the production human pangenome is published as Minigraph-Cactus and PGGB graphs.
- **Crop / livestock pangenomes:** plant and bovine super-pangenomes use pggb / minigraph-cactus.
- **Visualization and comparative genomics:** Bandage, sequenceTubeMap, odgi viz — operate directly on GFA.

---

## 3. The relationship between them

Despite being separate lineages, they are not unrelated. The recent paper *"From de Bruijn Graphs to Variation Graphs — Relationships Between Pangenome Models"* (Beller, Mäkinen et al., SPIRE 2023; arXiv 2503.14520) formalizes a constructive map from DBG to VG:

1. **Split** every node (k-mer) into a path of single-character nodes.
2. **Merge** nodes whose incident edge structure is identical, collapsing equivalent positions.
3. **Collapse** isolated edges and clean up.

The result is a variation graph that represents the same set of strings the DBG did, but with the redundancy removed and a coordinate system introduced. This is the formal sense in which a VG is "a DBG with the redundancy quotiented out plus path information added."

The reverse direction is also defined: any VG can be re-broken into a DBG of any chosen k by re-chopping its node labels and union-ing k-mers, but this throws away path identity.

The unifying abstraction proposed in some recent work is the **"sequence graph"** — `(N, E)` with sequence-labeled nodes and adjacency edges — under which both DBGs and VGs are special cases. A DBG is a sequence graph with all-equal-length labels and edges defined by overlap; a VG is a sequence graph with arbitrary-length labels, edges defined by adjacency, and an additional path set on top.

---

## 4. Side-by-side conceptual comparison

| Dimension | De Bruijn Graph | Variation Graph |
|---|---|---|
| **Era / origin** | de Bruijn 1946 (math); Idury–Waterman 1995, Pevzner 2001 (bio) | Cactus / POA roots ~2002–2011; Garrison 2018 (vg toolkit) |
| **Atomic unit** | k-mer (fixed length) | DNA segment (any length) |
| **Node label uniqueness** | Each label appears exactly once (collapse on label) | Labels may repeat across nodes (identity is the node ID) |
| **Edge definition** | Sequence-derived: (k−1) overlap | Topological: an adjacency observed in some embedded sequence |
| **Strand handling** | Bidirected via canonical k-mers, or two nodes per k-mer | Bidirected by definition; each node has L/R sides |
| **k parameter** | Mandatory, fixes resolution, brittle to choice | Absent (or only theoretical) |
| **Reference frame** | None — symmetric | Reference is a named path (usually privileged) |
| **Haplotype storage** | Not native — bolted on via colors | First-class — paths are part of `G` |
| **Cycles** | Forced by every repeat ≥ k | Only for true repetitive / tandem structure |
| **Determinism** | Same input + k → identical graph | Tool-dependent topology |
| **Build complexity** | Linear in input; fast (hours for human pangenome) | Alignment-driven; slow (days for human pangenome) |
| **Sequence-to-graph alignment** | Linear-time on the graph; well studied | NP-hard in general; needs haplotype index (GBWT) |
| **Representational redundancy** | High — each base appears in up to k overlapping nodes | Low — each base lives in one node |
| **Coordinate system** | None native — must build a separate index | Implicit via the reference path |
| **Visualization** | Hard — millions of nodes, dense cycles | Easier — ref path acts as backbone |
| **Standard format** | FASTA of unitigs + edge list, or BOSS/binary | GFA (S/L/P/W lines) |
| **Canonical assembly use** | De novo assembly from short reads | Reference-bias-free read mapping & SV genotyping |
| **Canonical pangenome use** | Bifrost, mdbg (k-mer-based pangenomes) | vg, pggb, minigraph(-cactus) (assembly-based pangenomes) |
| **Equivalent abstraction** | A "labeled multigraph with overlap edges" | A "bidirected sequence graph with embedded paths" |
| **Conceptual lineage** | Combinatorics → Eulerian-path assembly → k-mer indexing | Multiple sequence alignment → POA / Cactus → bidirected sequence graph + paths |

---

## 5. The ecosystem at a glance

```
                ┌────────────── Sequence graph G = (N, E) ──────────────┐
                │                                                       │
                │   nodes carry sequence labels, edges are adjacency    │
                │                                                       │
                ▼                                                       ▼
   ┌───────────────────────────────┐               ┌───────────────────────────────┐
   │   De Bruijn Graph             │               │   Variation Graph             │
   │   - all labels length k       │               │   - labels arbitrary length   │
   │   - edges = (k−1) overlap     │               │   - edges = observed adj.     │
   │   - parameter k mandatory     │               │   - bidirected, paths in G    │
   └───────────────────────────────┘               └───────────────────────────────┘
              │                                                       │
              ├─ Compacted DBG (BCALM2, Bifrost, Cuttlefish)           ├─ vg / vg construct
              ├─ Colored DBG (Cortex, Mantis, Bifrost)                 ├─ Minigraph (SV-only)
              ├─ Succinct DBG / BOSS (Bowe et al.)                     ├─ Minigraph-Cactus
              ├─ Variable-order DBG                                    ├─ PGGB (all-pairs)
              ├─ Probabilistic DBG (khmer)                             ├─ Cactus
              └─ Dynamic DBG (DynamicBOSS, Bifrost)                    └─ Seqwish (alignment→graph)

   Assembly:        SPAdes, Velvet, ABySS, MEGAHIT, Trinity            Read mapping:    vg giraffe, GraphAligner
   Indexing:        Mantis, Themisto, GGCAT                            SV genotyping:   vg call, Pangenie
   Variant calling: Cortex, DiscoSnp, McCortex                         Pangenomes:      HPRC, plant/livestock
```

DBGs are the workhorse of *building* sequence representations from raw data — assembly, k-mer indexing, multi-sample variant calling. VGs are the workhorse of *using* a known set of sequences as a richer, less biased reference — population-scale read mapping, pangenome-based SV genotyping, comparative genomics.

The two families are increasingly converging at the pangenome layer, where modern toolchains (Bifrost+vg, minigraph→Cactus→GBWT) cross the boundary deliberately: build with DBG-style speed, query with VG-style coordinate semantics.

---

## Sources

- [De Bruijn graph (Wikipedia)](https://en.wikipedia.org/wiki/De_Bruijn_graph)
- [How to apply de Bruijn graphs to genome assembly — Compeau, Pevzner & Tesler, Nat Biotech 2011](https://www.marcottelab.org/users/BCH339N_2016/DeBruijnPrimer.pdf)
- [Why are de Bruijn graphs useful for genome assembly? (PMC5531759)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5531759/)
- [De novo assembly and genotyping of variants using colored de Bruijn graphs — Iqbal et al., Nat Genet 2012 (Cortex)](https://www.nature.com/articles/ng.1028)
- [Succinct de Bruijn graphs — Bowe, Onodera, Sadakane, Shibuya (WABI 2012)](https://link.springer.com/chapter/10.1007/978-3-642-33122-0_18)
- [Succinct colored de Bruijn graphs (Bioinformatics 2017)](https://academic.oup.com/bioinformatics/article/33/20/3181/2995815)
- [Variable-Order de Bruijn Graphs — Boucher et al., arXiv 1411.2718](https://arxiv.org/abs/1411.2718)
- [Cuttlefish 2 — scalable compacted DBG construction (Genome Biology 2022)](https://link.springer.com/article/10.1186/s13059-022-02743-6)
- [GGCAT — fast compacted/colored DBGs (PMC10538363)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10538363/)
- [Comparison of OLC and DBG assembly classes (Briefings in Functional Genomics 2012)](https://academic.oup.com/bfg/article/11/1/25/191455)
- [Variation graph toolkit improves read mapping — Garrison et al., Nat Biotech 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6126949/)
- [vg toolkit — Variation graphs (wiki)](https://github.com/vgteam/vg/wiki/Variation-graphs)
- [Haplotype-aware graph indexes — Sirén et al. (PMC7223266)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7223266/)
- [Pangenome graphs — Eizenga et al., Annu Rev Genomics Hum Genet (PMC8006571)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8006571/)
- [ODGI: understanding pangenome graphs (Bioinformatics 2022)](https://academic.oup.com/bioinformatics/article/38/13/3319/6585331)
- [Minigraph — design and construction of reference pangenome graphs (Genome Biology 2020)](https://link.springer.com/article/10.1186/s13059-020-02168-z)
- [Minigraph-Cactus pangenome construction (Nat Biotech 2023)](https://www.nature.com/articles/s41587-023-01793-w)
- [Comparing methods for constructing and representing human pangenome graphs (Genome Biology 2023, PMC10691155)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10691155/)
- [From de Bruijn Graphs to Variation Graphs — Relationships Between Pangenome Models (arXiv 2503.14520)](https://arxiv.org/html/2503.14520)
- [On the Hardness of Sequence Alignment on De Bruijn Graphs (PubMed 36450127)](https://pubmed.ncbi.nlm.nih.gov/36450127/)
- [GFA format specification](http://gfa-spec.github.io/GFA-spec/GFA1.html)
- [Methods and Developments in Graphical Pangenomics (PMC8384392)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8384392/)
