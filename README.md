# snape

A somatic structural-variant caller that makes the somatic decision **before** any alignment.

Normal and tumour HiFi reads are parsed into LCP cores and counted into one 2-colour graph. A core
the tumour reads carry and the normal reads never do is somatic by construction — no reference, no
assembly, no aligner has been involved yet. Those cores are then cut out of the reads that carry
them, as short junction contigs with germline anchors on either side, and only those contigs go to
an aligner. On hcc1954 that is 1.87 Gbp of contigs standing in for 196 Gbp of tumour reads.

    make install-deps-lcptools      # build lib/lcptools (biolcp) in place
    make                            # -> ./snape
    ./snape run <normal.svb|fq> <tumour.svb|fq> <ref.fa> <dir> --threads 16

Eight stages, each skipped when its output is already in `<dir>`, so a kill and a restart resume:

    graph -> bubbles -> junctions -> minimap2 -> call (no normal)
          -> pick -> minimap2 -> call  =  <dir>/calls.vcf

Needs `minimap2` and `samtools` on PATH, and system htslib to build. Reads may be uncompressed FASTQ
or the 2-bit `.svb` format `snape pack` writes; every stage takes either.

## Where it stands

Against the Severus paper's own CASTLE truth sets, scored with the paper's own tool (Minda,
`--min_size 50 --tolerance 500`), Severus scored from its published output:

| pair | snape P / R / F1 | Severus P / R / F1 |
|---|---|---|
| hcc1954 | 0.895 / 0.847 / **0.870** | 0.921 / 0.826 / 0.871 |
| h2009 | 0.863 / 0.841 / 0.852 | 0.913 / 0.894 / 0.903 |
| hcc1937 | 0.854 / 0.683 / 0.759 | 0.859 / 0.795 / 0.826 |
| h1437 | 0.728 / 0.774 / 0.750 | 0.782 / 0.842 / 0.810 |

Every row is one run of the whole pipeline at its shipped defaults; nothing is tuned per pair. On
hcc1954 that is 0.001 of F1 short of Severus, ahead on recall and behind on precision. Every number, how it was measured, and
what was tried and thrown away: **`docs/RESULTS.md`**.

## Docs

- `docs/PIPELINE.md` — the eight stages and what each one does
- `docs/RESULTS.md` — the measured record, real data and simulation
- `docs/FOUNDATIONS.md` — LCP cores, signatures, and the structural limits
- `docs/CODE_STYLE.md` — house rules for the code
