#!/usr/bin/env bash
# The germline check asks the NORMAL READS at the reference locus, and it was being asked at 11x.
# Measured on hcc1954 by thinning the very same BAM: 3x leaves 290 false calls, 5.5x leaves 172,
# 11x leaves 136. The curve has not flattened, so the question was being asked with too little
# evidence, not answered wrongly. This aligns the rest of the normal and merges it in.
#
# Reads 1..2M are already in the existing BAM, so only the remainder is aligned and the two are
# merged rather than starting over.
set -uo pipefail
cd /home/yuio/svp
S=${1:-hcc1954}
REF=/home/datasets/references/genomes/hg38/hg38.fa
JX=/home/yuio/tmp/jx/$S
OLD=$JX/normal2M.bam
[ "$S" = hcc1954 ] && OLD=/home/yuio/tmp/fp_rootcause/normal_mirrors/n2M.bam
EXTRA=${EXTRA:-6000000}                      # reads to add on top of the 2M already aligned
OUT=$JX/normal_deep.bam

echo "[$S] +${EXTRA} okuma hizalaniyor  $(date +%T)"
./svp dump hifi/$S/n.svb | tail -n +8000001 | head -$((EXTRA * 4)) \
  | taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y --secondary=no "$REF" - 2> "$JX/mm2_deep.log" \
  | samtools sort -@ 4 -m 3G -T "$JX/stdeep" -o "$JX/normal_extra.bam" -
samtools index -@ 8 "$JX/normal_extra.bam"
echo "[$S] birlestiriliyor  $(date +%T)"
samtools merge -@ 8 -f "$OUT" "$OLD" "$JX/normal_extra.bam"
samtools index -@ 8 "$OUT"
echo "[$S] hazir: $(samtools idxstats "$OUT" | awk '{s+=$3} END{print s" hizalama"}')  $(date +%T)"
