# Long-read somatic SV callers: competitors for snape (as of 2026-10-02)

Scope: callers of somatic structural variants from long reads (PacBio HiFi and/or ONT), tumor/normal pairs first,
plus single-sample mosaic callers that can be made somatic by subtracting the normal. Versions were checked against
GitHub releases, bioconda and quay.io/biocontainers. "Runnable" means it can run on our HiFi BAMs aligned to
`ref/hg38_parmask.fa`, within the bench limits (cores 16-31, <= 16 threads, <= ~40 GB, nothing installed on the host).
The wrappers and their status on the micro pair are in [README.md](README.md).

## Severus paper (Keskus et al., Nat Biotechnol 2025) and the callers it benchmarked

Keskus A. et al. "Severus detects somatic structural variation and complex rearrangements in cancer genomes using
long-read sequencing." Nat Biotechnol 44:247-257 (online 2025-04-04), doi:10.1038/s41587-025-02618-8
(preprint bioRxiv 10.1101/2024.03.22.586313). HiFi was aligned with minimap2 2.24 (`-k 17 -y --eqx`).

Competitors in that paper (from the PMC full text):
- long-read somatic callers: **nanomonsv** 0.6.0, **SAVANA** 1.0.3;
- long-read germline callers made somatic: **Sniffles2** 2.2 (multisample mode; `--mosaic` was tried and lowered
  F1), **cuteSV**, **SVDSS** 2.0, **DeBreak** 1.2 (each run on tumor and normal and subtracted with Jasmine 1.1.5,
  `--dup_to_ins --max_dist 500`), **SVision-pro** 2.0 (its own somatic mode);
- short-read callers (Illumina): Manta 1.6.0, SvABA 1.2.0, GRIDSS2/GRIPSS 2.13.2/2.2.

## Table

| Tool | Version (date) | Paper / year | Input | T/N | License | bioconda / container | Runnable here |
|---|---|---|---|---|---|---|---|
| **Severus** | 1.7 (2026-03-12) | Keskus et al., Nat Biotechnol 2025, 10.1038/s41587-025-02618-8; github.com/KolmogorovLab/Severus | GRCh38 BAM (minimap2/pbmm2); phased VCF + haplotagged BAMs recommended; VNTR bed | T/N (tumor-only with a PoN) | BSD-3 | bioconda 1.7; quay.io/biocontainers/severus:1.7 | **yes** (`severus.sh` unphased, `severus_phased.sh` with Clair3 + LongPhase) |
| **SAVANA** | 1.3.8 (2026-06-30) | Elrick et al., Nat Methods 22:1436, 2025, 10.1038/s41592-025-02708-0; github.com/cortes-ciriano-lab/savana | GRCh38 BAM; haplotagging optional; SNP VCF only for CNAs | T/N (`savana to`: tumor-only) | Apache-2.0 | bioconda 1.3.8; quay.io/biocontainers/savana:1.3.8 | **yes** (`savana.sh`, `--pb`) |
| **nanomonsv** | 0.9.0 (2026-03-25) | Shiraishi et al., NAR 2023 (10.1093/nar/gkad526); github.com/friend1ws/nanomonsv | GRCh38 BAM/CRAM; racon or mafft | T/N (control strongly recommended) | GPL-3.0 | bioconda 0.9.0; quay.io/biocontainers/nanomonsv:0.9.0 | **yes** (`nanomonsv.sh`, `--qv25` HiFi preset) |
| **Sniffles2** | 2.8.0 here; 2.8.1 (2026-09-10) latest | Smolka et al., Nat Biotechnol 2024; mosaic mode improved in 2.7.0 (2025-10-22); github.com/fritzsedlazeck/Sniffles | any long-read BAM | no T/N mode: multisample (SUPP_VEC) or `--mosaic` + subtraction | MIT | bioconda 2.8.1; quay 2.8.0 | **yes** (`sniffles2.sh` CASTLE-style multisample; `sniffles2_mosaic.sh`) |
| **SVision-pro** | 2.5 (2025-09-04) | Wang et al., Nat Biotechnol 43:181, 2025 (online 2024-03), 10.1038/s41587-024-02190-7; github.com/songbowang125/SVision-pro | BAM + bundled LiteUnet model (.pth), CPU | T/N (`--detect_mode somatic`) | GPL-3.0 (+ README: non-commercial academic use) | bioconda 2.5; quay.io/biocontainers/svision-pro:2.5 | **yes** (`svision-pro.sh`; the biocontainer's pysam 0.9.1 crashes on breakend output, our image adds pysam 0.21) |
| **Delly** (`lr`) | 2.6.0 (2026-08-17) | Rausch et al., Bioinformatics 2012; long-read mode not separately published; github.com/dellytools/delly | BAM, `-y pb` / `-y ont`, exclude list | T/N (`delly filter -f somatic`) | BSD-3 | bioconda 2.6.0; quay.io/biocontainers/delly:2.6.0 | **yes** (`delly.sh`) |
| **minisv** | r29, no releases (last commit 2025-10-27, "WIP") | Qin, Heinz, Li, Cancer Res Commun 6(9):2282, 2026, 10.1158/2767-9764.CRC-25-0769; bioRxiv 10.1101/2025.10.28.685154; github.com/lh3/minisv | raw reads -> PAF/GAF with `ds:Z`: GRCh38 + CHM13 (minimap2 >= 2.28 `-cx map-hifi -s50 --ds`), CHM13 pangenome (minigraph 0.21 `-cxlr`), ideally tumor -> hifiasm assembly of the normal | T/N, tumor-only, mosaic | MIT | not on bioconda; own image (images/minisv) | **partly**: without the normal assembly (hifiasm on a 30-60x normal needs well over 40 GB); `minisv.sh` |
| **SomaSV** | 0.0.1 (no tags; commit 92f1560) | bioRxiv 10.64898/2026.02.16.706063 (2026-02); github.com/eioyuou/SomaSV | tumor + normal long-read BAMs; hybrid mode adds a normal short-read BAM + gnomAD-SV PoN | T/N | MIT (README) | pip only; own image (images/somasv) | **yes**, long-read-only mode (`somasv.sh`); hybrid mode needs short reads we do not use |
| **SVScope** | no releases (commit 0ad3686, 2026-03-27) | Tu & Xie, Genome Biol 27:186, 2026, 10.1186/s13059-026-04076-0; bioRxiv 10.1101/2025.02.11.636543; github.com/Goatofmountain/SVscope | minimap2 BAMs, hg38 RepeatMasker windows (bundled), POA (spoa/abPOA/mafft), shipped random forest | T/N | Apache-2.0 | conda env file; own image (images/svscope) | **wrapped, unverified** (`svscope.sh`): crashes at its last stage on both small test pairs (no inner-alignment candidate); its random forest was trained on HCC1395/H2009/H1437/HCC1954/HCC1937 labels: in-sample on our pairs |
| **longcallD** | 0.0.11 (2026-04-06) | bioRxiv 10.64898/2026.03.20.713111 (2026); github.com/yangao07/longcallD | single BAM; `-s` mosaic mode tags SOMATIC; SVs = INS/DEL only | single-sample mosaic; no T/N | MIT | bioconda; quay.io/biocontainers/longcalld:0.0.11 | **yes** as mosaic + subtraction (`longcalld.sh`), INS/DEL only |
| **colorSV** | 0.1.0 (2024-07-09) | Le & Li, Genomics Proteomics Bioinformatics qzaf082, 2025, 10.1093/gpbjnl/qzaf082; github.com/mktle/colorSV | no BAM: hifiasm co-assembly of pooled tumor+normal reads (raw unitig graph), read names must carry the sample, CHM13 | T/N, long-range events only (translocations, ~Mb) | none stated | not on bioconda | **no**: hifiasm co-assembly of ~60-120x pooled HiFi needs >> 40 GB and many hours; reports only large rearrangements |
| **SomaGauss-SV** | code not found | bioRxiv 10.1101/2025.06.12.659252 (2025-06) | ONT | T/N | ? | none | **no**: no code; ONT only |
| **EchoSV** (donor-specific assembly) | code not found | bioRxiv 10.1101/2025.10.28.685155; 10.64898/2026.02.20.707061 | diploid donor assembly from hybrid data | T/N | ? | none | **no**: no code; needs a diploid assembly |
| cuteSV | 2.1.4 (2026-07-30) | Jiang et al., Genome Biol 2020 | BAM | germline only (README points to SomaSV for somatic) | MIT | bioconda | not wrapped (germline; Severus paper used tumor-minus-normal subtraction) |
| SVDSS | 2.1.1 | Denti et al., Nat Methods 2023 | BAM | germline only | MIT | bioconda | not wrapped (germline) |
| DeBreak | 1.3 | Chen et al., Genome Biol 2023 | BAM | germline only | MIT | bioconda | not wrapped (germline) |
| Sawfish | 2.2.1 (2025-11-04) | Saunders et al., Bioinformatics 2025 | HiFi BAM | germline / joint only | PacBio license | bioconda | not wrapped (no somatic mode) |
| pbsv | 2.11.0 | PacBio | HiFi BAM | germline | PacBio license | bioconda | not wrapped (germline) |

Not somatic SV callers (checked): ClairS / ClairS-TO 0.5.1 and DeepSomatic 1.10 (SNV/indel only), HiFiCNV and Wakhan
(copy number), LongSom (single-cell long-read RNA), LRSomatic (Nextflow pipeline that runs Severus), EPI2ME
wf-somatic-variation (ONT, runs Severus).

## Notes

- **Learned components.** SAVANA ships random-forest classifiers (default ONT-trained; with `--pb` it uses PacBio
  thresholds, which is what we run). SVision-pro uses a shipped LiteUnet neural network. SVScope uses a shipped random
  forest trained on labelled somatic SVs of five of our benchmark cell lines (HCC1395, H2009, H1437, HCC1954,
  HCC1937): any score of SVScope on those pairs is in-sample for it and must be reported as such.
- **Other recent comparisons** (who benchmarked whom): SAVANA's paper: Sniffles2 2.2.0, cuteSV 2.1.0, SVIM 1.4.2,
  nanomonsv 0.5.0, Severus 1.0, SVision-pro 2.1, GRIDSS2. SVScope: nanomonsv 0.7.0, SAVANA 1.2.1, SVision-pro 1.8,
  Severus 1.6, Sniffles2 2.0.7 (+ SvABA, GRIPSS, Manta) on HCC1954, H2009, HCC1937, H1437, HCC1395, COLO829, HG008.
  minisv: Severus 1.6, SAVANA 1.3.6, Sniffles2 2.6.3, nanomonsv 0.8.0. colorSV: Severus, Sniffles2, nanomonsv, SAVANA.
  HG008 benchmark (Front Genet 2026, 10.3389/fgene.2026.1732039): Sniffles2 2.6.1, nanomonsv 0.8.0, SAVANA 1.3.5,
  Severus 1.5.
- **Minimum set every reviewer will expect**: Severus, SAVANA, nanomonsv, Sniffles2, SVision-pro (all in the Severus
  paper and most later ones); then the 2025-2026 newcomers that have code: minisv, SomaSV, SVScope; Delly lr and
  longcallD are extra.

## Cannot be run (and why)

- **colorSV**: input is a hifiasm co-assembly of tumor + normal reads, not BAMs; for 60-120x pooled HiFi this needs
  far more than 40 GB and many hours, and the method reports only long-range events (translocations / ~Mb SVs).
- **minisv, full recipe**: the recommended tumor/normal recipe aligns the tumor reads to a hifiasm assembly of the
  normal (same memory problem). We run it with GRCh38 + CHM13 + the CHM13 pangenome only (the README's documented
  fallback "when the normal assembly is not available").
- **SomaSV, hybrid mode**: needs a normal short-read BAM (and the gnomAD-SV PoN); we run its long-read-only mode.
- **SomaGauss-SV, EchoSV**: no public code (EchoSV also needs a diploid donor assembly; SomaGauss-SV is ONT-only).
- **Germline-only callers** (cuteSV, SVDSS, DeBreak, Sawfish, pbsv): no somatic mode; could only enter through
  tumor-minus-normal subtraction as in the Severus paper. Not wrapped.
