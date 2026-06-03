# Genome variation graph visualization — literature survey

Compiled 2026-06-03 for the `graph-inspector` branch. The goal: design a
WG-scale viz tool for the svp 2-colour LCP-core graph (85M nodes, 86M
edges, 113k bubbles). This doc records what's in the field and why we
picked our specific stack.

## Universal finding

**Nobody renders the whole human-WG graph at once.** Every mature tool
operates on extracted subgraphs (region, path neighbourhood, or
specific bubble). The svp inspector follows the same pattern:
`bubble_census` (WG) → `inspect` (per-subgraph) → `inspect_viz`
(per-subgraph PNG).

## Tools surveyed

### Bandage / Bandage-NG
- Repos: <https://github.com/rrwick/Bandage>, <https://github.com/asl/BandageNG>
- Paper: Wick et al., *Bioinformatics* 2015
- Input: GFA1 (plus FASTG, Graphviz); reads `KC:i` / `RC:i` for depth
- Layout: OGDF FMMM (Fast Multipole Multilevel), force-directed 2D
- Scale: practical for assembly subgraphs ≤ ~30k nodes; the full
  human-WG graph is impractical
- Color: per-node `CL:Z:<#RGB>` (custom tag, manual setting)
- Output: Qt GUI; PNG/SVG export
- BandageNG (2020+ fork) adds Hi-C links + multigraph support
- Status: BandageNG actively maintained; original Bandage dormant
- Verdict: **good as the optional interactive viewer for extracted
  subgraphs**; we generate Mode A GFAs that load cleanly

### ODGI (Optimized Dynamic Genome/Graph Implementation)
- Repo: <https://github.com/pangenome/odgi>
- Papers: Guarracino et al., *Bioinformatics* 2022; PG-SGD layout —
  Heumos et al., *Bioinformatics* 2024
- Commands of interest:
  - `odgi build -g x.gfa -o x.og` — convert GFA1 → ODGI's native `.og`
  - `odgi extract` — subgraph by region/path
  - `odgi sort` — 1D path-guided SGD ordering
  - `odgi layout` — 2D path-guided SGD layout
  - `odgi viz -i x.og -o x.png` — 1D banded PNG (path = row, position = column)
  - `odgi draw` — render 2D layout to PNG/SVG
- Layout: **Path-Guided Stochastic Gradient Descent (PG-SGD)** — samples
  pairwise genomic distances along paths, minimizes via SGD; O(m log m)
- Scale: tested on human pangenome (HPRC); gigabase-scale graphs
- Color: per-path coloring via CSV
- Status: actively maintained (HPRC core infrastructure)
- Verdict: **the canonical WG-scale renderer**; we call it as an
  optional subprocess when present

### gfaestus
- Repo: <https://github.com/chfi/gfaestus>
- Tech: Rust + Vulkan, GPU-accelerated
- Input: laid-out GFA1 (typically from `odgi layout`)
- Scale: large laid-out pangenome subgraphs interactively
- Status: last release 2022; codebase active but low-priority maintenance
- Verdict: nice-to-have if user has Vulkan; not a critical dependency

### Sequence Tube Map
- Repo: <https://github.com/vgteam/sequenceTubeMap>
- Paper: Beyer et al., *Bioinformatics* 2019
- Tech: D3.js web app
- Input: vg's JSON (xg + path data)
- Style: London Underground-inspired path bundles
- Scale: best for small regions (< few thousand nodes)
- Verdict: too vg-tied to repurpose; not used here

### gfaviz
- Repo: <https://github.com/ggonnella/gfaviz>
- Paper: Gonnella et al., *Bioinformatics* 2019
- Input: GFA1, GFA2
- Generic graph viz; no genome-specific features
- Verdict: not specialized enough for our needs

### MoMI-G (Modular Multi-scale Integrated Genome Graph browser)
- Repo: <https://github.com/MoMI-G/MoMI-G>
- Paper: Yokoyama & Suzuki, *BMC Bioinformatics* 2019
- Tech: web UI, custom input format
- Strength: hierarchical SV-centric inspection (Interval Card Deck)
- Verdict: interesting model but custom format = high integration cost

### PGR-TK (Pangenome Research Tool Kit)
- Repo: <https://github.com/Sema4-Research/pgr-tk>
- Paper: Chin et al., *Nature Methods* 2023
- Specialty: **principal bundle decomposition** for highly polymorphic
  loci (HLA, LPA, tandem-repeat hotspots) — re-linearizes the tangle
  into per-bundle SVGs
- Verdict: useful as a side tool for HLA / repeat regions; not yet
  wired into the inspector (could be a future "complex region" mode)

### vg viz
- Part of vg toolkit; SVG output of small graph regions
- Verdict: too small-scale; vg-tied; not used

### PanGraphViewer
- Paper: bioRxiv 2023
- Tech: Python + FCOSE layout
- Verdict: layout slow on large graphs; not used

### Assembly Graph Browser (AGB)
- Repo: <https://github.com/almiheenko/AGB>
- Paper: Mikheenko et al., *Bioinformatics* 2019
- Strength: automatic subgraph partitioning for large assembly graphs
- Verdict: assembly-focused; not pangenome-specialized

### GenomeRibbon
- Repo: <https://github.com/MariaNattestad/Ribbon>
- Paper: Nattestad & Schatz, *Bioinformatics* 2021
- Orthogonal: alignment-ribbon plots, not graph viz
- Useful for read-evidence display alongside graph (out of scope here)

## Algorithmic foundations

### Unitig compaction (Mode A)
- **BCALM2**: Drezen et al. — minimizer-based disk partitioning + unitig
  walk. Standard de Bruijn unitig definition: maximal 1-in-1-out chain
  excluding endpoints.
- **BiFrost**: Holley & Melsted — in-memory colored compacted de Bruijn
  graph; 15× faster than BCALM2 on colored graphs. Colors split unitigs
  at sample-set boundaries.
- For svp: maximal 1-in/1-out chain with **color-preserving** split (we
  do not merge across N → T transitions).

### Allele-difference (Mode B)
- **Cortex** (Iqbal et al. 2012, *Nature Genetics*): the canonical
  colored-DBG variant caller. The classify step trims shared prefix +
  shared suffix from the two allele paths; residual is the variant.
- For svp's Mode B: refinement on Mode A's output — the trimmed
  prefix/suffix is **absorbed into the L/R anchor super-nodes**, not
  emitted as separate "shared" segments.

### Superbubbles
- Onodera et al. 2013 (quadratic) / Brankovic et al. 2016 (O(n+m))
- Definition: maximal subgraph with single entrance + single exit (more
  general than simple bubble cycles)
- svp's `decompose()` uses simple BFS through non-shared nodes, NOT
  superbubble decomposition (per `feedback_bubble_decomposition` memory —
  superbubble + union-find over-merged adjacent loci, dropping recall
  from 99.62% to 93-96%)

### Path-Guided SGD (PG-SGD) layout
- Heumos et al. 2024
- Uses genome paths as anchors; samples pairwise genomic distances along
  paths; SGD minimizes a stress function
- Avoids the O(n²) cost of force-directed methods at the cost of
  requiring embedded paths
- Backing algorithm for ODGI's 1D + 2D layouts

### Force-directed (FMMM)
- Fruchterman-Reingold (1991) → FMMM (Hachul & Jünger 2004 — Fast
  Multipole Multilevel Method)
- Used by Bandage via OGDF; quality vs. iterations tradeoff (3 / 12 / 30
  / 60 / 120 iterations)
- Practical scale: ~ few thousand to tens of thousands of nodes

## GFA1 + custom tags convention

| tag | usage |
|---|---|
| `LN:i:<n>` | segment length in bases (standard) |
| `KC:i:<n>` | k-mer count (standard, Bandage uses it for depth color) |
| `RC:i:<n>` | read count (standard, Bandage uses it for line width) |
| `LB:Z:<label>` | custom display label (Bandage tooltip) |
| `CL:Z:<NORMAL|TUMOR|SHARED>` | our color label (custom) |
| `CC:Z:<#RGB>` | explicit RGB hex (custom, Bandage-compatible) |
| `PS:Z:<chrom>:<start>-<end>:<strand>` | posmap-position (custom) |
| `BI:i:<n>` | bubble id (custom) |
| `NI:i:<n>` | original svp node id (custom) |
| `NN:i:<n>` | how many LCP nodes folded into this super-node (Mode A/B) |
| `SV:Z:<EQ|D|I|S|C>` | bubble-level allele-diff tag (Mode B) |
| `INV_RC` (tag presence) | flag: T-allele was reverse-complemented before trim |

## Decision rationale for svp's inspector

1. **Render with matplotlib** as the primary path — zero install
   dependency, full control over annotation overlays (truth label, HLA
   badge, dead-end red border), small subgraphs (≤ few hundred unitigs
   after Mode A) fit easily in 2D matplotlib.
2. **GFA1 with custom tags** as the data interchange — both Bandage-NG
   and ODGI accept it; future plugins can layer on extra tools.
3. **Optional Bandage-NG / ODGI subprocesses** — silently invoked if the
   binaries are on PATH; their PNGs join the per-subgraph HTML index.
4. **Custom unitig compaction** (Mode A) — color-preserving, exact
   suffix/prefix overlap merge (NOT fixed `2*flank`; see
   `feedback_lcp_core_overlap`).
5. **Custom allele-diff** (Mode B) — refinement on Mode A; shared
   prefix/suffix absorbed into anchor super-nodes; tag set EQ/D/I/S/C
   per user spec (no MNP/CPLX separation).
6. **WG diagnostic via bubble_census** — full enumeration before any
   visualization; writes census.bin + census.tsv + bubble_summary.md
   so we can choose interesting bubbles to render.

## References

- Wick, R. R. et al. (2015) "Bandage: interactive visualization of de
  novo genome assemblies." *Bioinformatics* 31(20):3350–3352.
- Guarracino, A. et al. (2022) "ODGI: understanding pangenome graphs."
  *Bioinformatics* 38(13):3319–3326.
- Heumos, S. et al. (2024) "Pangenome graph layout by Path-Guided
  Stochastic Gradient Descent." *Bioinformatics* 40(1):btad891.
- Iqbal, Z. et al. (2012) "De novo assembly and genotyping of variants
  using colored de Bruijn graphs." *Nature Genetics* 44(2):226–232.
- Brankovic, L. et al. (2016) "Linear-time superbubble identification
  algorithm for genome assembly." *Theoretical Computer Science* 609:
  374–383.
- Chin, C.-S. et al. (2023) "Multiscale analysis of pangenomes enables
  improved representation of genomic diversity for repetitive and
  clinically relevant genes." *Nature Methods* 20:1213–1221.
- Beyer, W. et al. (2019) "Sequence tube maps: making graph genomes
  intuitive to commuters." *Bioinformatics* 35(24):5318–5320.
- Gonnella, G., Niehus, N., Kurtz, S. (2019) "GfaViz: flexible and
  interactive visualization of GFA sequence graphs." *Bioinformatics*
  35(16):2853–2855.
- Yokoyama, T., Suzuki, Y. (2019) "MoMI-G: modular multi-scale
  integrated genome graph browser." *BMC Bioinformatics* 20(1):548.
- Mikheenko, A., Kolmogorov, M. (2019) "Assembly Graph Browser:
  interactive visualization of assembly graphs." *Bioinformatics*
  35(18):3476–3478.
- Garrison, E. et al. (2024) "Building pangenome graphs." *Nature
  Methods* 21:2008–2012.
- Holley, G., Melsted, P. (2020) "Bifrost: highly parallel construction
  and indexing of colored and compacted de Bruijn graphs." *Genome
  Biology* 21:249.
