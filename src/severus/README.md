# Severus (vendored)

The parts of [Severus](https://github.com/KolmogorovLab/Severus) v1.7 that svp's contig arm calls
into, copied here so the calling pipeline is self-contained under `src/`. Upstream lives in
`lib/Severus`; nothing here is modified.

    severus.py           entry point
    severus/*.py         the 8 modules main.py imports -- bam_processing, breakpoint_finder,
                         build_graph, resolve_vntr, vcf_output, main, __init__, __version__
    vntrs/               the two human tandem-repeat beds we pass to --vntr-bed:
                           human_GRCh38_no_alt_analysis_set.trf.bed   (hg38)
                           human_hs37d5.bed                           (GRCh37, COLO829)

Deliberately NOT copied: `pon/` (the 1000G panels of normals -- hg38 and chm13 only, never used;
COLO829 is GRCh37 and has no panel) and the mouse/chm13 VNTR beds. They remain in `lib/Severus`
if a run ever needs them.

Invoked as: `python src/severus/severus.py --target-bam ... --out-dir ...`
