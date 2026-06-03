# inspect/assets — hg38 annotation BED files

UCSC hg38 tables, pulled via the curl commands below. Used to flag
in-region bubbles in the bubble census and overlay region badges on
rendered PNGs.

Loader normalizes chrom names: UCSC uses `chr1`/`chrX`/`chrM`; svp's
contigs use `1`/`X`/`M`. The loader strips the `chr` prefix on read.

Provenance (pulled 2026-06-03, run from this directory):

```bash
# HLA / MHC region (single record; reference Genome Reference Consortium)
printf "6\t28477797\t33448354\tMHC\n" > hg38_hla.bed

# Centromeres
curl -fsSL "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz" \
  | gunzip | awk 'BEGIN{OFS="\t"} {gsub(/^chr/,"",$2); print $2, $3, $4, $5}' \
  > hg38_centromere.bed

# Heterochromatin / telomere / centromere gaps
curl -fsSL "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz" \
  | gunzip | awk 'BEGIN{OFS="\t"} {gsub(/^chr/,"",$2); print $2, $3, $4, $8}' \
  > hg38_gap.bed

# Segmental duplications (filter to >=1 kb)
curl -fsSL "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz" \
  | gunzip | awk 'BEGIN{OFS="\t"} ($4-$3)>=1000 {gsub(/^chr/,"",$2); print $2, $3, $4, "segdup"}' \
  > hg38_segdup.bed

# Simple repeats (filter to >=100 bp)
curl -fsSL "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz" \
  | gunzip | awk 'BEGIN{OFS="\t"} ($4-$3)>=100 {gsub(/^chr/,"",$2); print $2, $3, $4, $5}' \
  > hg38_simpleRepeat.bed
```

(`hgdownload.cse.ucsc.edu` and HTTPS variants failed with timeouts in the
2026-06-03 pull; HTTP via `hgdownload.soe.ucsc.edu` worked.)

Actual pulled sizes (2026-06-03):

| file | size | records |
|---|---|---|
| hg38_hla.bed | 24 B | 1 |
| hg38_centromere.bed | 3.4 KB | 109 |
| hg38_gap.bed | 26 KB | 829 |
| hg38_segdup.bed | 1.9 MB | 69,870 |
| hg38_simpleRepeat.bed | 6.0 MB | 239,921 |

The simpleRepeat and segdup tables ship with the repo (~8 MB total committed)
so the inspector tool is self-contained — no UCSC dependency at runtime.
