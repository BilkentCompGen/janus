#!/usr/bin/env bash
# FAZ E — bir genomdan gercekci HiFi okumalari uretir (Badread).
#
#   sim_reads.sh <GENOME_FA> <OUT_DIR> <DEPTH_X> [<threads>]
#
# ORNEK: sim_reads.sh genomes/tumor/T01/genome_cnv.fa reads/T01 20 30
#
# ---------------------------------------------------------------------------------------
# DERINLIK TANIMI — bu isin en klasik hatasi
# ---------------------------------------------------------------------------------------
# Badread'in `--quantity Nx` degeri KENDISINE VERILEN FASTA'nin toplam uzunluguna goredir.
# Bizim FASTA'mizda her haplotip (ve WGD sonrasi her kopya) AYRI KONTIG oldugu icin,
# `--quantity 30x` demek lokus basina 30 x (kopya sayisi) kat demektir -- yani diploidde 60x.
#
# Dogru hedef: HAPLOID REFERANSIN her pozisyonunun DEPTH kez kaplanmasi:
#     hedef_toplam_baz = DEPTH x 3,054,815,472        (CHM13v2.0 chr1-22+chrX)
# Bu yuzden asagida `--quantity` BAZ CINSINDEN veriliyor, "x" ile DEGIL.
#
# Ploidi 2.9 olan bir yaprakta bu, kopya basina DEPTH/2.9 kat demektir; CN=4 bolge dogal
# olarak CN=2 bolgenin iki kati okuma alir. Gercekci olan da budur.
# ---------------------------------------------------------------------------------------
# PARAMETRELER (spec/parameters.yaml reads:, docs/05)
#   --error_model / --qscore_model pacbio2021
#   --identity 28.3,3   -> DELESYON oranina gore eslenmis (SV/indel benchmark'i icin dogru hedef).
#                          docs/05 §3.2: --identity yalnizca TOPLAM hatayi olcekler; sub:ins:del
#                          karisimi 7-mer tablosuna gomulu, dolayisiyla tek bir deger birden fazla
#                          sinifi ayni anda eslestiremez ve HANGI sinifin eslendigi soylenmelidir.
#   --length 17000,4400 -> gercek Revio HiFi uzunluk dagilimi
#   OLCULDU (2026-08-18, 400 Mbaz kalibrasyonu): 137.3 kb/s/cekirdek, 0.869 bayt/baz,
#   ortalama uzunluk 17,253 / N50 17,948 / medyan kimlik %99.857.
# ---------------------------------------------------------------------------------------
# IS BIRIMI: KONTIG x MIKTAR PARCASI
# ---------------------------------------------------------------------------------------
# Kontig BASINA ayri kosu sart: badread referansi tumuyle bellege alir; tum genomu (8.9 Gb)
# 30 surece vermek ~270 GB RSS demek olurdu. En buyuk kontig ~250 Mb, yani 30 x 250 MB = 7.5 GB.
# ASLA kontig-ALTI bolmeyin (docs/05 §5.1): referansi kesmek parca sinirlarinda kapsamayi
# COKERTIR.
#
# Ama MIKTARI bolmek referansi bolmek DEGILDIR: ayni tam kontig uzerinde, farkli tohumla,
# hedefin bir parcasi kadar okuma uretilir. Okumalar bagimsiz cekimler oldugu icin parcalarin
# birlesimi tek bir buyuk kosuyla istatistiksel olarak ayni seydir. Bunun iki faydasi var:
#
#   1. COKME DAYANIKLILIGI. 2026-08-18'de bir reboot, o an yazilan 26 shard'in HEPSINI kesti
#      ve her biri ~1.3 Gbaz'lik tek parca oldugu icin ~1.5 saat is coop gitti. Parca 250 Mbaz
#      olunca cokmede cekirdek basina en fazla ~30 dakika kaybedilir.
#   2. YUK DENGESI. chr1 izi (~250 Mb) tek parca olsaydi 1.7 Gbaz uretmesi ~3.4 saat surer ve
#      ornegin son saatlerinde diger cekirdekler bos kalirdi. Esit parcalar bunu bitirir.
#
# ATOMIK YAZIM: cikti once .part'a yazilir, badread+gzip 0 donerse nihai ada tasinir. Nihai
# adin varligi = o parcanin eksiksiz oldugu garantisi. Onceki surum "dosya bos degilse tamam"
# sayiyordu ve yukaridaki reboot'ta 26 KESIK gzip'i sessizce kabul edecekti.
set -uo pipefail

GEN="${1:?genome fa}"
OUT="${2:?out dir}"
DEPTH="${3:?depth x}"
NT="${4:-24}"
SEEDBASE="${SEEDBASE:-20260818}"
CHUNK="${CHUNK:-250000000}"     # parca basina hedef baz (~30 dk/cekirdek)

BADREAD_PY=$(cd "$(dirname "$0")/.." && pwd)/.venv/bin/python3
HAP_REF_BP=$(python3 "$(dirname "$0")/refinfo.py" | awk '/haploid/{gsub(",","",$3); print $3}')
                               # T2T-CHM13v2.0 chr1-22 + chrX, refinfo.py'den (gomulu sayi YOK)

[ -x "$BADREAD_PY" ] || { echo "Badread venv python yok: $BADREAD_PY" >&2; exit 1; }

# TEK KOSU KILIDI (flock). 2026-08-19 02:21'de AYNI cikti dizinine IKI sim_reads birden
# baslatildi (benim baslatma hatam) ve ikisi de ayni .part dosyalarina yaziyordu -- yani
# sessizce bozuk gzip uretiliyordu. Artik ikinci kosu rc=3 ile HEMEN cikiyor.
mkdir -p "$OUT"
exec 9>"$OUT/.lock"
if ! flock -n 9; then
  echo "!! $OUT icin baska bir sim_reads kosuyor -- bu kosu baslatilmiyor (rc=3)" >&2
  exit 3
fi
[ -s "$GEN.fai" ] || samtools faidx "$GEN"
mkdir -p "$OUT/shards" "$OUT/ref"

TOTAL_TARGET=$(python3 -c "print(int($DEPTH * $HAP_REF_BP))")
FASTA_BP=$(awk '{s+=$2} END{print s}' "$GEN.fai")
echo "genom      : $GEN"
echo "FASTA bazi : $FASTA_BP"
echo "hedef baz  : $TOTAL_TARGET   (= ${DEPTH}x x haploid referans ${HAP_REF_BP})"
echo "  -> FASTA'nin $(python3 -c "print(round($TOTAL_TARGET/$FASTA_BP,3))") kati"

# --- gorev listesi: kontig x parca, uzunlukla orantili, parca <= CHUNK ---
python3 - "$GEN.fai" "$TOTAL_TARGET" "$FASTA_BP" "$CHUNK" "$SEEDBASE" > "$OUT/tasks.tsv" <<'PY'
import sys, zlib, math
fai, tot, fb, chunk, seedbase = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
for line in open(fai):
    f = line.split("\t")
    c, L = f[0], int(f[1])
    q = int(tot * L / fb)
    if q < 1000:
        continue
    k = max(1, math.ceil(q / chunk))
    base = q // k
    for j in range(k):
        qq = base + (q - base * k if j == k - 1 else 0)
        seed = seedbase + zlib.crc32(f"{c}#{j}".encode()) % 2000000000
        print(f"{c}\t{j}\t{qq}\t{seed}")
PY
NTASK=$(wc -l < "$OUT/tasks.tsv")
echo "gorev      : $NTASK parca (<= $CHUNK baz), $NT paralel"

# --- kontig FASTA'larini BIR KEZ cikar (parcalar paylasiyor) ---
if [ ! -s "$OUT/ref/.done" ]; then
  cut -f1 "$OUT/tasks.tsv" | sort -u | while read -r c; do
    [ -s "$OUT/ref/$c.fa" ] || samtools faidx "$GEN" "$c" > "$OUT/ref/$c.fa"
  done
  touch "$OUT/ref/.done"
fi

one_piece(){
  set -o pipefail
  c=$(echo "$1" | cut -f1); j=$(echo "$1" | cut -f2)
  q=$(echo "$1" | cut -f3); seed=$(echo "$1" | cut -f4)
  o="$OUT/shards/${c}.p${j}.fastq.gz"
  [ -s "$o" ] && return 0
  p="${o}.part"
  "$BADREAD_PY" -m badread simulate \
      --reference "$OUT/ref/${c}.fa" --quantity "$q" \
      --error_model pacbio2021 --qscore_model pacbio2021 \
      --identity 28.3,3 --length 17000,4400 \
      --seed "$seed" 2>> "$OUT/shards/${c}.p${j}.log" | gzip -1 > "$p"
  rc=$?
  if [ $rc -ne 0 ] || [ ! -s "$p" ]; then
    echo "  ${c}.p${j} BASARISIZ (rc=$rc)"
    return 1
  fi
  mv -f "$p" "$o"          # <-- ancak burada "tamamlanmis" sayilir
  echo "  ${c}.p${j} tamam ($(stat -c %s "$o") bayt)"
}
export -f one_piece; export OUT BADREAD_PY

DONE0=$(ls "$OUT"/shards/*.fastq.gz 2>/dev/null | wc -l)
[ "$DONE0" -gt 0 ] && echo "  ($DONE0 parca zaten tamam, atlanacak)"
cat "$OUT/tasks.tsv" | xargs -d'\n' -P "$NT" -I{} bash -c 'one_piece "$@"' _ {}

NDONE=$(ls "$OUT"/shards/*.fastq.gz 2>/dev/null | wc -l)
if [ "$NDONE" -ne "$NTASK" ]; then
  echo "!! $NDONE/$NTASK parca tamam -- EKSIK, summary yazilmiyor" >&2
  exit 2
fi

# --- olculmus ozet (hedefe guvenme, SAY) ---
python3 - "$OUT" <<'PY' | tee "$OUT/summary.txt"
import glob, gzip, os, sys
d = sys.argv[1]
n = bp = 0
lens = []
for f in sorted(glob.glob(os.path.join(d, "shards", "*.fastq.gz"))):
    with gzip.open(f, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                n += 1; L = len(line) - 1; bp += L; lens.append(L)
lens.sort()
def n50(v):
    h = sum(v)/2; s = 0
    for x in reversed(v):
        s += x
        if s >= h: return x
    return 0
print(f"parca sayisi : {len(glob.glob(os.path.join(d,'shards','*.fastq.gz')))}")
print(f"okuma sayisi : {n:,}")
print(f"toplam baz   : {bp:,}")
print(f"ortalama uz. : {bp//max(n,1):,}")
print(f"medyan uz.   : {lens[len(lens)//2] if lens else 0:,}")
print(f"N50          : {n50(lens):,}")
PY

# --- kontig FASTA'lari artik gereksiz ---
rm -rf "$OUT/ref"
echo "bitti: $OUT"
