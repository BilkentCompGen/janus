# CLAUDE.md

_Project-specific guidance. Global behavioural rules + cross-project map: `/home/claude/CLAUDE.md`._

---

`snape` — a somatic SV caller that makes the somatic decision **before** alignment. One normal and N
tumour HiFi read sets go into one (N+1)-colour LCP-core graph; a core no tumour-free sample carries
is somatic by construction, with no reference and no assembly involved, and which tumour lines carry
it is read straight off the count planes. Those cores are cut out of the reads that carry them as
junction contigs, per line, and only the contigs (plus the ~1.5% of normal reads that answer the
germline question) are ever aligned — the normal once, however many tumours there are. C++17, one
binary, one Makefile, no CI.

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
"SNPG" | u32 ver=2 | u32 S | u32 resv | u64 K | u64 sig[K] | u8 cov[S][K]
                                     | u64 E | u32 esrc[E] | u32 edst[E] | u8 ecov[S][E]
```

- `sig` = MurmurHash64A over the core's **homopolymer-compressed** bytes. The mixer and the
  compression are both part of the on-disk contract — change either and every sig-keyed artifact
  (`graph.bin`, `drop.sigs`) has to be rebuilt together or they silently match nothing.
- **Colour is derived, never stored:** plane 0 is the normal, so `cov[0][i]==0` is somatic and the
  non-zero tumour planes ARE the support pattern — which lines carry the variant, with no coordinate
  involved. S <= 64, checked at the CLI boundary.
- **Counts saturate at 255.** No stage reads their magnitude: three of four uses are zero tests or a
  floor test, the fourth is behind a filter that ships disabled. Measured on hcc1954, clipping at
  255 leaves the surviving core set and the chunk set identical. One typedef, `cov_t` in `snape.h`.
- Node id = index into the sig-sorted survivor array. Edges carry per-sample traversal counts.
- Both strands of every read go in — the parse is not RC-symmetric, so each locus has a mirror twin
  with a different signature and no core in common. A large inversion's interior is therefore
  *shared*, not somatic; what is callable is its two breakpoint junctions.

## Scoring — do not mix the two instruments

- **Real data**: Minda, the Severus paper's own tool at the paper's settings —
  `minda truthset --base hifi/<pair>/truth_clean.vcf --vcfs <calls>.vcf --min_size 50 --tolerance 500`.
  Severus is scored the same way from its **published** output. Use the `minda_4_2` ensemble; the
  looser ones exist only for hcc1954. That is the same *command* on both sides but **not a symmetric
  comparison**: four of the compared callers are members of the ensemble and a Severus vote is in
  93.9–97.3% of the records. For a head-to-head, score both sides against `truth_no_<other>` from
  `/home/yuio/tmp/loo/` (built by `loo_truth.py` from the ensemble's own SUPP_VEC); snape is in no
  ensemble, so that is the only construction in which neither side helped build the event set.
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
