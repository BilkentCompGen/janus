# svp

Reference-free somatic structural-variant caller. Builds a 2-colour LCP-core graph directly from
normal + tumor HiFi reads (no assembly, no alignment, no reference), error-corrects it, and calls
somatic variants off the coordinate-free graph.

Three resumable, checkpointed stages:

```
svp graph  <normal.fq.gz> <tumor.fq.gz> <dir>    reads       -> graph.bin
svp refine <dir>                                 graph.bin   -> refined.bin
svp call   <dir>                                 refined.bin -> calls.tsv + calls.altsigs.bin
svp run    <normal.fq.gz> <tumor.fq.gz> <dir>    all three, skipping completed stages
```

See `docs/PIPELINE.md` for the stage algorithms and `docs/FOUNDATIONS.md` for the LCP-core
foundations. `CLAUDE.md` maps the project.
