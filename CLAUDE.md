# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

SVP (Structural Variant Pipeline) detects somatic structural variants (DEL, INS, DUP, INV, SNP) by building a Locally Consistent Parsing (LCP) graph from reference, normal, and tumor genomes, then finding divergent paths (bubbles) between them.

## Build & Run

```bash
make                # build the `svp` binary (g++, C++20, -O3)
make svp_clean      # build the `svp_clean` preprocessor binary (tools/svp_clean.cpp)
make help           # print available make targets
make clean          # remove build artifacts

# run
./svp --ref ref.fa --normal normal.fa --tumor tumor.fa --depth 5 --chr 19 --output out

# tests (each is a standalone C++ file in tests/)
make test findpath          # build & run a single test
make test cycledup
make test refanchor
make test svcall
make test node
make test hash

# diagnostic/trace tools (in tests/diagnostics/, not regression tests)
make diag trace_bfs
make diag dup_diag

# VCF comparison
make compare base.vcf comp.vcf ref.fa    # runs truvari via scripts/truvari.sh
python scripts/compare_vcf.py sim.vcf detected.vcf   # position-tolerant VCF diff
python scripts/deep_analysis.py                       # position-level cross-reference with type checking

# GFA visualisation: build/compress/strip via Docker odgi
./scripts/odgi.sh <graph_name>   # reads outputs/gfa/<graph_name>_graph.gfa
```

The `make` machinery treats trailing words as build targets, so `make test foo` builds and runs `tests/foo.cpp` (and likewise `make diag bar` for `tests/diagnostics/bar.cpp`). Test source files are standalone — they include the same `src/` modules and link the same libs as `svp` itself.

**Dependencies:** g++20, zlib (`-lz`), htslib (`-lhts`), parasail (`-lparasail`), lcptools (`-llcptools` — bundled in `lib/lcptools/`).

## Algorithm & Complexity

See [docs/ALGORITHM.md](docs/ALGORITHM.md) for a detailed breakdown of every pipeline phase, the algorithms used (LCP parsing, KMP merging, Needleman-Wunsch alignment, BFS cycle detection), and their time/space complexity. Overall: **O(N + B·L² + S·K)**.

## Background reading

[docs/GRAPH_MODELS.md](docs/GRAPH_MODELS.md) — conceptual reference on de Bruijn graphs and variation graphs: formal definitions, the variant families (compacted / colored / succinct / bidirected DBGs; vg, minigraph, minigraph-cactus, pggb VGs), the GFA format, the formal DBG↔VG relationship, and the wider pangenomics ecosystem. No project-specific content; pure background for choosing or migrating between graph models.

## Architecture

**Data flow:** FASTA inputs → LCP cores → colored graph → bubble detection → anchor to reference → SV calling → VCF output + GFA graph.

### Core modules (src/)

- **graph.h/cpp** — `LCPGraph` and `Node`. Nodes carry color bits: `0x1` = normal, `0x2` = tumor, `0x3` = shared. Key methods: `addString()`, `detectBubbles()`, `extractBubbles()`, `findColoredPath()` (BFS with KMP-based sequence merging), `writeGFA()`.
- **ref.h/cpp** — `RefIndex`. Maps LCP core hashes to reference positions. `anchorBubble()` walks backward/forward through predecessors to find a unique reference anchor for each bubble.
- **sv.h/cpp** — `SVCall` and SV calling logic. `callSV()` checks for INV (reverse-complement identity ≥0.90) first, then `detectSV()` does quick classification; `callMulti()` uses parasail Needleman-Wunsch alignment + CIGAR parsing for complex variants; `extractAndCallSVs()` orchestrates the full pipeline including flipped DUP detection. `detectLargeINVs()` is a separate core-index pipeline for inversions too large for bubble-based detection — classifies tumor cores against forward-ref and revcomp-ref hash sets, finds runs of rcref-only cores (≥5 cores, ≥500bp).
- **pipeline.h/cpp** — `processGenome()` (FASTA→LCP→graph), `buildRefHashes()` (reference core hash set + RefIndex population), `buildRcrefHashes()` (reverse-complement reference hash set for large INV detection).
- **vcf.h/cpp** — `writeVCF()`. VCF output formatting and INS→DUP Smith-Waterman reclassification.
- **diag.h/cpp** — `DiagLog`, `DiagBubble`, `DiagCall`. Optional per-bubble / per-call TSV trace, gated on `--debug`. Writes anchor method, skip reasons, emitted call types — used for post-hoc analysis when accuracy regresses.
- **utils.h/cpp** — `Logger` (debug macro), `FastaStreamer` (gzip-aware FASTA iterator built on `zstr.hpp`), `ArgParser`, `hash_string()`, `reverse_complement()`.
- **zstr.hpp** — bundled header-only gzip stream wrapper around zlib; consumed by `FastaStreamer` so `.fa` and `.fa.gz` work transparently.

### Key concepts

- **Bubbles**: Divergent paths in the graph where normal and tumor sequences differ. Detected by shared→divergent→shared node color transitions.
- **Anchoring**: Mapping bubble positions back to reference genome coordinates via LCP core hashes. Walks through graph predecessors/successors to find uniquely-mapped cores.
- **INV detection (bubble-based)**: Bubbles where `tumor_seq ≈ reverse_complement(normal_seq)` (SW score ≥0.90, min 10bp) are classified as inversions. Uses Smith-Waterman local alignment to handle LCP core boundary shifts.
- **INV detection (core-index)**: For large inversions that don't form clean bubbles. Classifies each tumor LCP core as forward-ref-only (cls=1), rcref-only (cls=2), both (cls=3), or novel (cls=0). Runs of rcref-only cores (gap-tolerant, MAX_GAP=3, ≥5 rcref cores, ≥500bp span) are called as INVs. Anchors to reference via flanking ref-only cores. Runs before bubble extraction; overlapping bubble calls are suppressed.
- **INS→DUP conversion**: Insertions with ≥70% Smith-Waterman similarity to adjacent reference are reclassified as duplications.
- **Flipped DUP detection**: Medium/large DUPs produce bubble starts near the DUP end and ends near the DUP start (positions flipped). Detected by the `T≥1 S≥1 N=0` successor/predecessor color pattern on anchored bubble boundaries.

### main.cpp

Thin orchestrator: parses CLI args (`--chr` defaults to 19), calls `buildRefHashes()` / `buildRcrefHashes()` / `processGenome()`, runs `detectLargeINVs()` + `extractAndCallSVs()`, then `writeVCF()` and `writeGFA()`.

### tools/svp_clean.cpp

Standalone two-stage preprocessor for **real-data SV calling** (built via `make svp_clean`). The LCP graph is sensitive to single-base noise: a real tumor genome has thousands of SNPs / small indels that fragment cores and inflate the bubble count. `svp_clean` strips them out before `svp` runs.

Algorithm: divide ref and tumor into fixed-size windows (FNV-1a hashed), match identical windows as anchors, then for each gap between anchors, replace the tumor gap with the ref gap when their sizes differ by ≤ `max_indel` bp (i.e. SNP / small-indel). Larger gaps are preserved unchanged so `svp` can still see them as SVs. Fixed windows beat LCP cores here because LCP boundaries shift around any nearby SNP — fixed windows have stable boundaries.

### Experimental sibling: `svp2/`

`svp2/` holds an alternative **6-stage modular pipeline** (s1_lcp → s2_graph → s3_bubbles → s4_simple.py → s5_align.py → s6_finalize.py) developed on the `experimental/modular-pipeline` branch. The full implementation, extractors (`chain_rt`, `chain_3w`, `bfs`, `bulge`, `superbubble`), aligners (`wfa`, `wfa_rc`, `edlib`, `parasail`, `mappy`), grid-search infrastructure, and per-sim F1 numbers live there. **`svp2/STATUS.md` is the authoritative entry point** — read it before touching that directory. On `main` the directory is largely a stub; the populated tree lives on the experimental branch.

## Simulations (insilicoSV)

Ground truth simulations live in `simulation/` and are generated by [insilicoSV](lib/insilicoSV/) (`pip install .` from `lib/insilicoSV/`, or use the project venv: `source .venv/bin/activate`).

```bash
# Generate a simulation (run from simulation/ directory)
cd simulation
insilicosv -c <config>.yaml
# Outputs: sim.fa, sim.vcf, sim.novel_insertions.fa, sim.divergence.fa
# Rename to match naming convention, e.g.:
mv sim.fa <type>_simulation.fa && mv sim.vcf <type>_simulation.vcf
```

**Existing configs** (all in `simulation/`):
- `normal_tumor_deletions.yaml` → `deletion_simulation.{fa,vcf}`
- `normal_tumor_insertions.yaml` → `insertion_simulation.{fa,vcf}`
- `normal_tumor_duplications.yaml` → `duplication_simulation.{fa,vcf}`
- `normal_tumor_snps.yaml` → `snp_simulation.{fa,vcf}`
- `normal_tumor_inversions.yaml` → `inversion_simulation.{fa,vcf}`
- `normal_tumor_mixed.yaml` → `mixed.{fa,vcf}` (DEL + INS + DUP + INV + SNP)

**Creating a new simulation type:**
1. Create a YAML config in `simulation/` following existing configs as templates.
2. Key fields: `reference` (path to `ref.fasta`), `overlap_regions` (`chr19.bed`), `variant_sets` (type, number, length_ranges), `haploid: True`, `min_intersv_dist`.
3. Run `insilicosv -c <config>.yaml` from the `simulation/` directory.
4. Rename outputs from `sim.*` to `<type>_simulation.*`.
5. The crash at `output_stats()` is a known bug — outputs (fa, vcf) are already written before it hits.

**Reference**: `simulation/ref.fasta` (human chr19, ~58MB). All simulations use this as the base genome. The ref also serves as the "normal" genome in benchmarks since variants are tumor-only.

## Benchmarking

### SVP benchmarking tools

- **`scripts/compare_vcf.py`** — Position-tolerant (50bp) VCF diff. Reports TP/FP/FN by category (small <50bp, del >=50bp, ins >=50bp, inv >=50bp). Does not understand DUP type.
- **`scripts/deep_analysis.py`** — Position-level (200bp) cross-reference with strict type matching. Reports correct-type, wrong-type, ghost FPs, missed FNs with size breakdowns. The primary accuracy analysis tool.
- **`scripts/truvari.sh`** — [Truvari](https://github.com/ACEnglish/truvari) benchmark for SV-level comparison (>=50bp). Reports TP/FP/FN with precision/recall/F1. Requires `bcftools`, `bgzip`, `tabix`, `truvari` (install via `pip install truvari` or use project venv).

```bash
# Truvari benchmark (SV-level, requires sorted/indexed VCFs)
./scripts/truvari.sh <base.vcf> <comp.vcf> <ref.fa> [output_dir]
# or via make:
make compare base.vcf comp.vcf ref.fa

# Position-tolerant diff (all variant sizes)
python scripts/compare_vcf.py <simulated.vcf> <detected.vcf>

# Deep cross-reference analysis (strict type matching, all runs)
python scripts/deep_analysis.py
```

### Real-data benchmarks

`benchmarks/real/` holds real-genome benchmarks (COLO829 tumor/normal cell line is the primary somatic SV truth set). `download.sh` fetches truth VCFs and demo data; `prepare.sh` indexes the reference and assembles HiFi BAMs into FASTA via hifiasm; `run_svp.sh` runs SVP on the prepared sample. See `benchmarks/real/README.md` for the full flow and disk/RAM requirements (full COLO829 needs ~120 GB disk and ~64 GB RAM for assembly). Sample FASTAs are gitignored.

Rival-caller comparisons (Delly, Manta, etc.) referenced from `benchmarks/real/README.md` live in a sibling `benchmarks/rivals/` tree on branches where it was wired up — that directory is not present on `main`.

Historical benchmark results are in [docs/BENCHMARKS.md](docs/BENCHMARKS.md). A presentation summary lives in [docs/PRESENTATION.md](docs/PRESENTATION.md) (`SVP_Presentation.pptx` is the slide deck).

## Output

- VCF files in `outputs/vcf/`
- GFA graph files in `outputs/gfa/`
- Truvari benchmark results in `outputs/truvari/`

## Directory layout

```
main.cpp              # thin orchestrator for the `svp` binary
src/                  # core modules (graph, ref, sv, pipeline, vcf, diag, utils, zstr.hpp)
tools/                # standalone tools (svp_clean.cpp — real-data SNP/indel scrubber)
tests/                # unit/integration tests (make test <name>)
tests/diagnostics/    # investigative trace/debug tools (make diag <name>)
scripts/              # production helpers (truvari.sh, compare_vcf.py, deep_analysis.py, odgi.sh)
docs/                 # ALGORITHM.md, BENCHMARKS.md, PRESENTATION.md, slide deck
simulation/           # insilicoSV configs + generated truth data (chr19)
benchmarks/real/      # real-data benchmarks (COLO829 etc., gitignored FASTAs)
svp2/                 # experimental 6-stage modular pipeline — see svp2/STATUS.md (lives on `experimental/modular-pipeline`)
lib/                  # bundled dependencies (lcptools, parasail, insilicoSV)
outputs/              # gitignored: vcf/, gfa/, og/, truvari/
```
