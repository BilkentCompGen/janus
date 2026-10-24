#!/usr/bin/env bash
# Kromozom basina uretilen haplotipleri TEK bir diploid genoma birlestirir ve
# germline truth VCF'ini kurar.
#
#   assemble_germline.sh <GERMLINE_DIR>
#
# Cikti:
#   <DIR>/germline.fa        46 kontig: chr1_hapA, chr1_hapB, ... chrX_hapA, chrX_hapB
#   <DIR>/germline.fa.fai
#   <DIR>/germline.truth.vcf.gz   tum kromozomlarin import VCF'lerinin birlesimi (hg38 koordinati)
#   <DIR>/assembly_report.txt
set -euo pipefail
D="${1:?germline dir}"
D=$(cd "$D" && pwd)
CH=$(ls -d "$D"/chr* 2>/dev/null | grep -E '/chr[0-9XY]+$' | sort -V)
[ -z "$CH" ] && { echo "kromozom dizini yok" >&2; exit 1; }

OUT="$D/germline.fa"
: > "$OUT"
n=0
for c in $CH; do
  name=$(basename "$c")
  for h in hapA hapB; do
    f="$c/sim.$h.fa"
    if [ -s "$f" ]; then cat "$f" >> "$OUT"; n=$((n+1));
    else echo "EKSIK: $f" >&2; fi
  done
done
samtools faidx "$OUT"

# --- truth VCF: her kromozomun import VCF'i zaten hg38 koordinatinda ve acik GT'li ---
TV="$D/germline.truth.vcf"
{
  echo '##fileformat=VCFv4.2'
  echo "##source=phylo/build_germline.sh (HGSVC3 v1.0 PAV, HG00733)"
  awk '/^##contig=/{print}' "$D/import/$(basename "$(echo $CH | awk '{print $1}')").vcf" 2>/dev/null || true
  for c in $CH; do
    awk '/^##contig=/{print}' "$D/import/$(basename $c).vcf" 2>/dev/null
  done | sort -u
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">'
  echo '##INFO=<ID=END,Number=1,Type=Integer,Description="END">'
  echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">'
  echo '##INFO=<ID=INSSEQ,Number=1,Type=String,Description="INSSEQ">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00733\n'
  for c in $CH; do
    grep -v '^#' "$D/import/$(basename $c).vcf" 2>/dev/null || true
  done
} > "$TV"
bgzip -f "$TV" && tabix -f -p vcf "$TV.gz"

{
  echo "germline assembly raporu — $(date '+%Y-%m-%d %H:%M:%S')"
  echo
  echo "kontig sayisi        : $(wc -l < "$OUT.fai")   (beklenen: 2 x kromozom sayisi)"
  echo "toplam uzunluk       : $(awk '{s+=$2} END{printf "%.3f Gb", s/1e9}' "$OUT.fai")"
  echo "truth VCF kayit      : $(bcftools index -n "$TV.gz")"
  echo
  echo "kontig                uzunluk"
  awk '{printf "  %-18s %13d\n", $1, $2}' "$OUT.fai"
  echo
  echo "haplotip uzunluk farklari (hapA - hapB, net indel etkisi):"
  awk '{L[$1]=$2} END{for(k in L){ if(k ~ /_hapA$/){b=k; sub(/_hapA$/,"_hapB",b);
       if(b in L) printf "  %-10s %+8d\n", substr(k,1,length(k)-5), L[k]-L[b]}}}' "$OUT.fai" | sort -V
} > "$D/assembly_report.txt"
cat "$D/assembly_report.txt"
