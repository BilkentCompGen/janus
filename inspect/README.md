# inspect — graph-inspector for the svp somatic SV caller

WG-scale diagnostic + visualization toolkit for the 2-colour LCP-core graph
produced by `svp` (branch `lcp-cortex-merged`). Lives on its own branch
`graph-inspector` so `lcp-cortex-merged` source stays clean.

Design doc: `~/.claude/plans/now-read-all-of-humble-spring.md`.

## Build

```
cd inspect && make
```

Depends on the parent `svp` build for `lib/lcptools/liblcptools.a`,
`build/io.o`, `build/cores.o`, `build/MurmurHash3.o`. The Makefile invokes
the parent build for the .o files automatically.

## Pipeline

### Stage 0 — schema probe (run once)
```
./inspect/seqs_probe ../outputs/svp/wgs
```
Confirms the byte layout of `seqs.bin`. Required reading before anything
that consumes node sequences.

### Stage 1 — WG bubble census
```
./inspect/bubble_census ../outputs/svp/wgs ../simulation/truth/wgs/out.vcf assets/
```
Walks all 112,977 bubbles. Outputs `census.{bin,tsv}` + `bubble_summary.md`
in the same dir. ~10 min.

### Stage 2 — per-subgraph extract
```
./inspect/inspect ../outputs/svp/wgs sv ../simulation/truth/wgs/out.vcf sv17978
./inspect/inspect ../outputs/svp/wgs bubble 42057
./inspect/inspect ../outputs/svp/wgs region 11:102600000-102700000
```
Each emits `inspect/<key>/{raw,modeA,modeB}.gfa + meta.json`.

### Stage 3 — render
```
python3 inspect/src/inspect_viz.py ../outputs/svp/wgs/inspect/sv_sv17978
```
PNGs (topology + linear) + `<key>/index.html`. Aggregate report at
`outputs/svp/wgs/inspect/index.html`.

## Two compaction modes

- **Mode A** = maximal 1-in/1-out chain compaction. Adjacent LCP-core
  sequence windows are spliced by **exact suffix/prefix byte match**
  (variable overlap; the LCP cores themselves overlap by ~85 bp on average
  at L5, so the merge windows overlap by 100+ bp — NEVER `2*flank`).
- **Mode B** = refinement on Mode A. Within each bubble, the shared
  prefix + suffix of (N-allele, T-allele) is absorbed into the L/R anchor
  super-nodes. The N-path / T-path super-nodes shrink to the actual diff;
  empty diffs collapse the edge entirely. Bubble tag set: EQ / D / I / S / C.

See `~/.claude/projects/-home-yuio-svp/memory/feedback_lcp_core_overlap.md`
and `feedback_modeB_on_modeA.md` for the rationale.

## Annotations

`assets/` holds small UCSC hg38 BED files (HLA, centromere, segdup, simple
repeat). These are used to flag in-region bubbles in the census and on
the rendered PNGs. The tools work without them — flags just default to
false.

This is an "hg38 cheat" — the simulation's normal is hg38 ref, so we
abuse it for design-time inspection. Real tumor/normal pairs without ref
identity wouldn't have these directly.
