# CLAUDE.md

_Project-specific guidance. Global behavioural rules + cross-project map: `/home/claude/CLAUDE.md`._

---

`snape` — a somatic SV caller that makes the somatic decision **before** alignment. Normal and tumour
HiFi reads go into one 2-colour LCP-core graph; a core the tumour carries and the normal never does
is somatic by construction, with no reference and no assembly involved. Those cores are cut out of
the reads that carry them as junction contigs, and only the contigs (plus the ~1.5% of normal reads
that answer the germline question) are ever aligned. C++17, one binary, one Makefile, no CI.

    ./snape run <normal.svb|fq> <tumour.svb|fq> <ref.fa> <dir> --threads 16

## Hard rules

- **Reads only.** The two sample inputs are read sets (`.svb` or uncompressed FASTQ, a file or a
  directory of shards). An assembly is a ground-truth oracle, never a caller input.
- **Whole genome.** Scope comes from the data, never hard-coded.
- **16 cores.** Always `--threads 16` + `taskset -c 0-15`. The box shows 32 logical CPUs; 16 are
  ours, and oversubscribing makes every wall-time incomparable to the record.
- **`/home/datasets` is read-only** (hg38 lives at `references/genomes/hg38/hg38.fa`). Real read
  data and truth sets are under `hifi/`, untracked; scratch goes to `/home/yuio/tmp`.
- Never `rm`. Move to `/home/yuio/trash/<what>-<date>/`.
- A new version of a file is not a new file: edit in place, no `_v2` / `_new` / `_fixed`.

## What is where

    src/            the whole pipeline, one binary (see docs/PIPELINE.md for the stage map)
    docs/           PIPELINE (stages) · RESULTS (every measurement) · FOUNDATIONS (LCP facts) · CODE_STYLE
    scripts/        measurement only: the simulation read generator, the junction-level sim scorer,
                    the HiFi error profiler. Nothing here is part of calling.
    paper/          main.tex, the method paper (LNCS). Build it with the texlive image:
                    docker run --rm --user $(id -u):$(id -g) -v "$PWD":/work -w /work davp-latex \
                        -c "pdflatex -interaction=nonstopmode main.tex"   (twice, for references)
    hifi/           get.sh fetches the CASTLE pairs; the reads/truth it writes are untracked
    lib/            third-party, built in place (lcptools `biolcp` is the only build dependency)
    simulation/     insilicoSV configs + generated truth/reads (untracked)
    outputs/        run artifacts (untracked)

## The graph, in one box

```
u64 K | u64 sig[K] | u32 ncov[K] | u32 tcov[K] | u64 E | CEdge{u32 src,dst,ncov,tcov}[E]
```

- `sig` = MurmurHash64A over the core's **homopolymer-compressed** bytes. The mixer and the
  compression are both part of the on-disk contract — change either and every sig-keyed artifact
  (`graph.bin`, `drop.sigs`) has to be rebuilt together or they silently match nothing.
- **Colour is derived, never stored:** `ncov==0 && tcov>0` somatic, `tcov==0 && ncov>0` normal-only,
  both shared.
- Node id = index into the sig-sorted survivor array. Edges carry per-colour traversal counts.
- Both strands of every read go in — the parse is not RC-symmetric, so each locus has a mirror twin
  with a different signature and no core in common. A large inversion's interior is therefore
  *shared*, not somatic; what is callable is its two breakpoint junctions.

## Scoring — do not mix the two instruments

- **Real data**: Minda, the Severus paper's own tool at the paper's settings —
  `minda truthset --base hifi/<pair>/truth_clean.vcf --vcfs <calls>.vcf --min_size 50 --tolerance 500`.
  Severus is scored the same way from its **published** output, so the comparison is symmetric. Use
  the `minda_4_2` ensemble; the looser ones exist only for hcc1954.
- **Simulation**: `scripts/severus_bp_eval.py`, junction-level against insilicoSV's own
  `adjacencies.bed` (a deletion is one junction, an inversion two).
- Truvari is **not** used any more, in either place.

## Working on this

- Every default in the code was measured, and the measurement is written next to it or in
  `docs/RESULTS.md`. Change one and re-measure on real data AND the simulation — the two disagree
  often, and `hard` has no VNTRs to punish a repeat filter.
- Do not rebuild the binary in the middle of a sweep. Freeze a copy under `/home/yuio/tmp` and sweep
  with that.
- `graph.bin` is not byte-reproducible in its junk-core tail (map merge order); the surviving core
  set is deterministic.
- `-DLARGE_GENOME` must match `lib/lcptools`'s build or core structs read back corrupt.
- Do not commit without being asked.
