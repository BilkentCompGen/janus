#!/usr/bin/env bash
# FAZ C — gercekci diploid germline genomu kur.
#
#   build_germline.sh <SAMPLE> <SRC_VCF_DIR> <OUT_DIR> [<n_parallel>]
#
# Girdi : HGSVC3 v1.0 T2T-CHM13 callset'leri (PAV 2.4.0.1, 65 fazli genom, DIZI-COZUNURLU).
#         GRCh38 surumunun birebir karsiligi; yalnizca koordinat sistemi farkli.
#         SNV + indel + SV(insdel) + INV, hepsi h1|h2 fazli.
# Cikti : <OUT_DIR>/
#           <chrom>/sim.hapA.fa, sim.hapB.fa, sim.vcf, sim.hap*.paf, adjacencies.bed
#           import/<chrom>.vcf         (insilicoSV'ye verilen kesin varyant listesi)
#           verify/<chrom>.txt         (BAGIMSIZ dizi dogrulamasi -- GECTI olmali)
#           germline.fa                (48 kontig: chr*_hapA + chr*_hapB)
#
# Tasarim gerekcesi: experiments/01 (zigosite), experiments/03 (ucdan uca dogrulama),
# docs/07_architecture.md (kontig = haplotip).
#
# Kromozom basina ayri kosu: bellek 24 kat duser, 24 cekirdege dagilir ve bir kromozomdaki
# hata digerlerini bozmaz. Germline varyantlari kromozomlar arasi gecmez, bu yuzden gecerlidir.
set -uo pipefail

SAMPLE="${1:?sample}"
SRC="${2:?hgsvc3 dir}"
OUT="${3:?out dir}"
NP="${4:-11}"
# Yollar MUTLAK olmali: asagida her kromozom icin kendi dizinine `cd` ediliyor ve goreli bir
# $OUT o noktadan sonra cozulemiyor (yeniden yonlendirme sessizce patliyor, insilicoSV hic
# calismiyor, rc=1). 2026-08-18'de tam olarak bu oldu.
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
SRC=$(cd "$SRC" && pwd)
PHYLO=$(cd "$(dirname "$0")/.." && pwd)
REF="$PHYLO/ref/chm13.primary.fa"
PY=$PHYLO/.venv/bin/python
ISV=$PHYLO/.venv/bin/insilicosv
SEED=20260818

mkdir -p "$OUT"/{import,verify,logs,chrref}
export SAMPLE SRC OUT PHYLO REF PY ISV SEED

one_chrom(){
  # ZATEN TAMAMLANMIS KROMOZOMU ATLA. Onceki surumde atlama yoktu: zincir yeniden
  # baslatildiginda 23 kromozomun hepsi bastan kuruluyordu (~24 dk bosa is). Tamamlanma
  # olcutu log'daki TAMAM satiri -- yani hem dizi yazilmis hem BAYT-BAYT dogrulama gecmis
  # olmali; yarim kalan bir dizin bu testi gecemez.
  if [ -s "$OUT/logs/$1.log" ] && grep -q "$1 TAMAM" "$OUT/logs/$1.log" 2>/dev/null \
     && [ -s "$OUT/$1/sim.hapA.fa" ] && [ -s "$OUT/$1/sim.hapB.fa" ]; then
    echo "SKIP $1 (zaten tamam)"; return 0
  fi
  c="$1"
  log="$OUT/logs/$c.log"
  : > "$log"
  echo "[$(date +%H:%M:%S)] $c basladi" >> "$log"

  # --- 1. kromozom referansi ---
  cr="$OUT/chrref/$c.fa"
  if [ ! -s "$cr" ]; then
    samtools faidx "$REF" "$c" > "$cr" 2>>"$log"
    samtools faidx "$cr" 2>>"$log"
  fi

  # --- 2. dort callset'ten bu ornegin bu kromozomdaki varyantlarini topla ---
  raw="$OUT/import/${c}.raw.vcf"
  if [ ! -s "$raw" ]; then
    # Kaynak dosyalar: HGSVC3 cok-ornekli callset'ler, ya da zaten cikarilmis tek-ornek VCF.
    srcs=$(ls "$SRC"/variants_T2T-CHM13_snv_snv_alt_*.vcf.gz \
              "$SRC"/variants_T2T-CHM13_indel_insdel_alt_*.vcf.gz \
              "$SRC"/variants_T2T-CHM13_sv_insdel_alt_*.vcf.gz 2>/dev/null)
    [ -z "$srcs" ] && srcs=$(ls "$SRC"/${SAMPLE}.autosomes.vcf.gz "$SRC"/${SAMPLE}*.vcf.gz 2>/dev/null | head -1)
    for f in $srcs; do
      [ -s "$f" ] || continue
      bcftools view -r "$c" -s "$SAMPLE" -Ou "$f" 2>>"$log" \
        | bcftools view -a -c1 -Ov - 2>>"$log"
    done | awk 'BEGIN{h=0} /^##/{if(!h) print; next} /^#CHROM/{if(!h){print; h=1} next} {print}' \
         > "$raw.tmp" 2>>"$log"
    bcftools sort -Ov -o "$raw" "$raw.tmp" 2>>"$log" && mv -f "$raw.tmp" "$raw.tmp.sorted-input"
  fi
  n_raw=$(grep -vc '^#' "$raw" 2>/dev/null || echo 0)
  echo "[$(date +%H:%M:%S)] $c ham varyant: $n_raw" >> "$log"

  # --- 3. insilicoSV import formatina cevir ---
  imp="$OUT/import/${c}.vcf"
  "$PY" "$PHYLO/scripts/vcf_to_insilicosv.py" "$raw" "$imp" \
      --drop-overlaps --seqless-ins drop --ref "$cr" \
      --report "$OUT/import/${c}.dropped.tsv" >> "$log" 2>&1
  n_imp=$(grep -vc '^#' "$imp" 2>/dev/null || echo 0)
  echo "[$(date +%H:%M:%S)] $c import kaydi: $n_imp" >> "$log"

  # --- 4. insilicoSV ---
  d="$OUT/$c"; mkdir -p "$d"
  cat > "$d/config.yaml" <<YAML
reference: "$cr"
random_seed: $SEED
allow_hap_overlap: True
output_paf: True
output_adjacencies: True
variant_sets:
    - import: "$imp"
YAML
  ( cd "$d" && /usr/bin/time -v "$ISV" -c config.yaml >> "$log" 2> "$d/time.log" )
  rc=$?
  if [ $rc -ne 0 ] || [ ! -s "$d/sim.hapA.fa" ]; then
    echo "[$(date +%H:%M:%S)] $c BASARISIZ rc=$rc" >> "$log"; tail -20 "$d/time.log" >> "$log"
    echo "FAIL $c" ; return 1
  fi
  grep -E 'Maximum resident set size|Elapsed \(wall' "$d/time.log" >> "$log"

  # --- 5. BAGIMSIZ dizi dogrulamasi (uretecin raporuna bakmadan) ---
  "$PY" "$PHYLO/scripts/verify_haplotypes.py" --ref "$cr" \
      --hapA "$d/sim.hapA.fa" --hapB "$d/sim.hapB.fa" \
      --import-vcf "$imp" --chrom "$c" > "$OUT/verify/$c.txt" 2>&1
  vrc=$?
  if [ $vrc -ne 0 ]; then
    echo "[$(date +%H:%M:%S)] $c DOGRULAMA KALDI" >> "$log"; echo "VERIFY-FAIL $c"; return 2
  fi
  echo "[$(date +%H:%M:%S)] $c TAMAM ($n_imp varyant, dogrulama GECTI)" >> "$log"
  echo "OK $c $n_raw $n_imp"
}
export -f one_chrom

CHROMS="${CHROMS:-$(echo chr{1..22} chrX)}"
echo "kromozomlar: $CHROMS"
printf '%s\n' $CHROMS | xargs -P "$NP" -I{} bash -c 'one_chrom {}'
echo "=== ozet ==="
grep -h 'TAMAM\|BASARISIZ\|KALDI' "$OUT"/logs/*.log 2>/dev/null | sort
