# Real-Data Benchmarks

Benchmark SVP on real genomes with validated ground truth.

## Quick start

```bash
# 1. Download truth sets + demo data + reference (~1GB)
./download.sh quick

# 2. Index reference
./prepare.sh index-ref

# 3. Run SVP on demo data (fast test)
./run_svp.sh colo829-demo
```

## Full COLO829 benchmark

```bash
# 1. Download HiFi BAMs (~120GB — need disk space)
./download.sh colo829-hifi

# 2. Assemble into FASTA for SVP + align BAMs for rival callers
./prepare.sh assemble-colo829    # needs hifiasm, ~64GB RAM, hours
./prepare.sh prep-rival-bams     # needs minimap2, hours

# 3. Run SVP
./run_svp.sh colo829

# 4. Run rival callers (see ../rivals/)
cd ../rivals && ./run_all.sh
```

## Datasets

### COLO829 (melanoma cell line) — somatic SV benchmark

The gold standard for somatic SV caller benchmarking. Matched tumor (COLO829) and normal (COLO829BL) cell lines with a validated truth set of 68 somatic SVs.

| Resource | Source | Size |
|----------|--------|------|
| Truth set VCF (2022) | [Zenodo 4716169](https://zenodo.org/records/4716169) | ~1MB |
| Truth set VCF (2024 update) | [Zenodo 10819636](https://zenodo.org/records/10819636) | ~1MB |
| PacBio Revio HiFi BAMs (60X) | [PacBio downloads](https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829/) | ~120GB |
| colorSV demo FASTAs (subset) | [GitHub release](https://github.com/mktle/colorSV/releases/tag/v0.1.0) | ~50MB |

**References:**
- Valle-Inclan et al. 2022, *Cell Genomics* — [doi:10.1016/j.xgen.2022.100167](https://doi.org/10.1016/j.xgen.2022.100167)
- Closing the gaps, *Genome Research* 2025 — updated truth set with GRCh38/T2T coordinates

**Note:** No pre-assembled whole-genome FASTA exists for COLO829. SVP needs FASTA input, so the HiFi BAMs must be assembled with hifiasm first (`prepare.sh assemble-colo829`). Rival callers use the aligned BAMs directly.

### HG002 (GIAB Ashkenazi son) — germline SV benchmark

Well-characterized germline sample with a Tier 1 SV truth set (12,745 INS+DEL >= 50bp). T2T diploid assembly available as direct FASTA input for SVP.

| Resource | Source | Size |
|----------|--------|------|
| Tier 1 SV truth set (v0.6) | [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/) | ~5MB |
| T2T diploid assembly (hg002v1.1) | [Human Pangenomics S3](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz) | ~2GB gz |

**Caveat:** GIAB v0.6 truth set uses GRCh37 coordinates. A liftover to GRCh38 is needed for accurate comparison. This is a germline (not somatic) benchmark — SVP uses GRCh38 as "normal" and HG002 as "tumor" to detect germline variants.

### GRCh38 reference genome

| Resource | Source | Size |
|----------|--------|------|
| No-alt analysis set | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/) | ~850MB gz |

## Directory structure

```
real/
  download.sh       # download datasets (modular, by name)
  prepare.sh        # index, assemble, extract chromosomes
  run_svp.sh        # run SVP and compare against truth
  samples/          # FASTAs, BAMs, assemblies (gitignored)
  ground_truth/     # truth set VCFs
  results/          # SVP output VCFs + timing
```

## Available download.sh commands

```bash
./download.sh --list            # show all datasets
./download.sh quick             # truth sets + demo + GRCh38 (~1GB)
./download.sh truth             # just truth sets
./download.sh colo829-hifi      # full HiFi BAMs (~120GB)
./download.sh all               # everything
```
