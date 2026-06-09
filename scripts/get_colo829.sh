#!/usr/bin/env bash
# get_colo829.sh — fetch the COLO829 real somatic HiFi benchmark (PacBio Revio 2023Q2).
#   tumor  = COLO829     (2 hifi_reads BAMs, ~128 GB) -> colo829/tumor/COLO829_tumor_hifi.fastq.gz
#   normal = COLO829-BL  (2 hifi_reads BAMs, ~134 GB) -> colo829/normal/COLO829BL_normal_hifi.fastq.gz
#   truth  = Valle-Inclan 2022 GRCh38-lifted 68-SV VCF (Zenodo 7515830)
#
# Resumable + idempotent:
#   - wget -c resumes partial BAMs; a BAM is re-fetched only if its byte size != server Content-Length.
#   - per-sample FASTQ conversion is skipped if the .fastq.gz already exists.
#   - after a sample converts OK, its source BAMs are deleted to cap peak disk (~400 GB vs ~520 GB).
# Usage: scripts/get_colo829.sh [threads=8]    (run in screen/background; ~262 GB download)
set -euo pipefail
ROOT=/home/yuio/svp
DIR="$ROOT/colo829"
THREADS="${1:-8}"
B=https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829/COLO829
N=https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829/COLO829-BL
TUMOR_BAMS=(m84039_230312_025934_s1.hifi_reads.bc2026.bam m84039_230328_000836_s3.hifi_reads.bc2025.bam)
NORMAL_BAMS=(m84039_230327_230708_s1.hifi_reads.bc2007.bam m84039_230327_233730_s2.hifi_reads.bc2007.bam)
TRUTH_URL="https://zenodo.org/records/7515830/files/truthset_somaticSVs_COLO829_hg38lifted.vcf?download=1"

mkdir -p "$DIR"/{tumor,normal,truth}
log(){ echo "[colo829 $(date +%F_%T)] $*"; }

remote_size(){ curl -sIL "$1" | awk 'tolower($1)=="content-length:"{v=$2} END{gsub(/\r/,"",v); print v}'; }
local_size(){ [ -f "$1" ] && stat -c %s "$1" || echo 0; }

fetch_bam(){ # url_base filename destdir
  local url="$1/$2" dest="$3/$2"
  local want; want=$(remote_size "$url")
  local have; have=$(local_size "$dest")
  if [ -n "$want" ] && [ "$have" = "$want" ]; then log "have $2 ($have B) -> skip"; return 0; fi
  log "fetch $2 (have ${have} / want ${want:-?} B)"
  wget -c -q --show-progress -O "$dest" "$url"
  have=$(local_size "$dest")
  [ -z "$want" ] || [ "$have" = "$want" ] || { log "SIZE MISMATCH $2: $have != $want"; return 1; }
  log "done $2 ($have B)"
}

convert(){ # sample_label destdir out_fastq bam...
  local label="$1" destdir="$2" out="$3"; shift 3
  if [ -s "$out" ]; then log "$label fastq exists -> skip convert"; return 0; fi
  log "$label: BAM -> $out (samtools fastq, $THREADS threads, pigz)"
  for f in "$@"; do samtools fastq -@ "$THREADS" "$destdir/$f"; done | pigz -p "$THREADS" > "$out.tmp"
  mv -f "$out.tmp" "$out"
  log "$label fastq: $(du -h "$out" | cut -f1)"
  log "$label: removing source BAMs to reclaim disk"
  for f in "$@"; do rm -f "$destdir/$f"; done
}

log "START threads=$THREADS  dest=$DIR"
df -h /home | tail -1

# truth first (tiny)
if [ -s "$DIR/truth/truthset_somaticSVs_COLO829_hg38lifted.vcf" ]; then
  log "truth VCF exists -> skip"
else
  log "fetch truth VCF"
  wget -c -q -O "$DIR/truth/truthset_somaticSVs_COLO829_hg38lifted.vcf" "$TRUTH_URL"
fi

# tumor
for f in "${TUMOR_BAMS[@]}"; do fetch_bam "$B" "$f" "$DIR/tumor"; done
convert tumor "$DIR/tumor" "$DIR/tumor/COLO829_tumor_hifi.fastq.gz" "${TUMOR_BAMS[@]}"

# normal
for f in "${NORMAL_BAMS[@]}"; do fetch_bam "$N" "$f" "$DIR/normal"; done
convert normal "$DIR/normal" "$DIR/normal/COLO829BL_normal_hifi.fastq.gz" "${NORMAL_BAMS[@]}"

log "DONE"
log "  tumor : $DIR/tumor/COLO829_tumor_hifi.fastq.gz"
log "  normal: $DIR/normal/COLO829BL_normal_hifi.fastq.gz"
log "  truth : $DIR/truth/truthset_somaticSVs_COLO829_hg38lifted.vcf"
df -h /home | tail -1
