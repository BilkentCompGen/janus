<!-- Real public somatic-tumor PacBio HiFi datasets + benchmark SV truth sets.
     Compiled by background research agent, 2026-06-08. For testing svp on REAL data
     as an alternative to the in-silico Sim-it benchmark. -->

# Real public somatic tumor PacBio HiFi datasets (+ SV truth)

> The famous COLO829/HCC1395 PacBio data on ENA/SRA (PRJEB27698 / SRP162370) is **Sequel
> CLR, not HiFi**. The genuine HiFi came later from PacBio's public **Revio 2023Q2** release
> — that's what the Severus/SAVANA/nanomonsv papers used. All raw HiFi ships as PacBio
> unaligned `hifi_reads.bam` (convert with `samtools fastq` / `bam2fastq`), not `.fastq.gz`.

## Recommended pick: COLO829 (melanoma) — turnkey for somatic SV

The field's canonical somatic-SV reference. Hand-curated, orthogonally validated (PCR /
capture / Bionano) truth; the exact set Severus/SAVANA/nanomonsv report F1 against
(Severus F1≈0.79 on PacBio HiFi). Open Revio HiFi for tumor + matched normal, GRCh38 truth on
Zenodo, no application.

- **HiFi reads (open, Revio, GRCh38):** `https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829/`
  - Tumor `COLO829/`: `m84039_230312_025934_s1.hifi_reads.bc2026.bam` (66 GB) + `m84039_230328_000836_s3.hifi_reads.bc2025.bam` (53 GB) — ~37×
  - Normal `COLO829-BL/`: `m84039_230327_230708_s1.hifi_reads.bc2007.bam` (61 GB) + `m84039_230327_233730_s2.hifi_reads.bc2007.bam` (64 GB) — ~40×
- **Truth VCF (open, 68 somatic SVs: 38 DEL/7 DUP/7 INV/3 INS/13 TRA):** Valle-Inclán 2022, *Cell Genomics*. Zenodo `10.5281/zenodo.7515830`:
  - GRCh38 (use this): `https://zenodo.org/records/7515830/files/truthset_somaticSVs_COLO829_hg38lifted.vcf?download=1`
  - Refreshed 2025 T2T/GRCh38 benchmark (54 stable SVs): Zenodo `10.5281/zenodo.14367730`
- Contigs un-prefixed (`1`,`2`,…) — `bcftools annotate --rename-chrs` if your ref is `chr`-prefixed.

```bash
mkdir -p colo829/{tumor,normal,truth} && cd colo829
B=https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829/COLO829
wget -c "$B/m84039_230312_025934_s1.hifi_reads.bc2026.bam" -P tumor/
wget -c "$B/m84039_230328_000836_s3.hifi_reads.bc2025.bam" -P tumor/
N=https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829/COLO829-BL
wget -c "$N/m84039_230327_230708_s1.hifi_reads.bc2007.bam" -P normal/
wget -c "$N/m84039_230327_233730_s2.hifi_reads.bc2007.bam" -P normal/
for f in tumor/*.bam;  do samtools fastq -@8 "$f"; done | gzip > tumor/COLO829_tumor_hifi.fastq.gz
for f in normal/*.bam; do samtools fastq -@8 "$f"; done | gzip > normal/COLO829BL_normal_hifi.fastq.gz
wget -c "https://zenodo.org/records/7515830/files/truthset_somaticSVs_COLO829_hg38lifted.vcf?download=1" \
     -O truth/truthset_somaticSVs_COLO829_hg38lifted.vcf
```

## Alternatives

- **HG008 (GIAB, pancreatic)** — best breadth + actively-curated: open Revio HiFi (tumor ~116×, normal ~35×) at `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Revio_20240125/`; NIST draft SV benchmark v0.3 (149 SVs + high-conf BED) on GRCh38-GIABv3 / GRCh37 / CHM13. Pre-aligned BAMs provided. SNV/indel somatic truth still in development.
- **HCC1395 (SEQC2, breast)** — largest truth (1,788 SVs: 717 DEL/230 DUP/551 INS/133 INV/146 TRA). Open Revio HiFi at `https://downloads.pacbcloud.com/public/revio/2023Q2/HCC1395/`; consensus SV set lives in the *Genome Biology* 2022 paper supplement (DOI 10.1186/s13059-022-02816-6), not a clean FTP VCF.
- **CASTLE panel (Severus)** — 6 cancer pairs, mostly dbGaP controlled (PRJNA1086849); skip unless you need extra lines.

## Sources
PacBio Revio COLO829/HCC1395: https://downloads.pacbcloud.com/public/revio/2023Q2/ ·
COLO829 truth Zenodo 7515830 / 14367730 · Valle-Inclán 2022 (PMC9903816) ·
GIAB HG008: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/ ·
HG008 SV benchmark v0.3 · SEQC2 HCC1395 (PMC9746098) · Severus (PMC12483193).
