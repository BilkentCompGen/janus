#!/usr/bin/env bash
# Faz D'yi SEVIYE SEVIYE paralel kosar.
#
#   run_tree_parallel.sh <GERMLINE_FA> <OUT_DIR> [<max_parallel>]
#
# Agac 39 dugum ve 9 derinlik; sirali kosarsak 39 x ~6 dk = ~4 saat. Ama bir dugum yalnizca
# EBEVEYNINI bekler, kardeslerini degil. Bu yuzden dugumler seviyeye gore gruplanip her seviye
# paralel kosuluyor: kritik yol ~10 adim.
#
# BELLEK SINIRI: her dugum ebeveyn genomunu Python'a yukluyor (~16 GB RSS 6 Gb'lik genom icin).
# 157 GB RAM'de guvenli ust sinir 6 -- ama WGD sonrasi genom ikiye katlanacagi icin varsayilan 4.
set -uo pipefail
GERM="${1:?germline fa}"
OUT="${2:?out dir}"
NP="${3:-4}"
PHYLO=$(cd "$(dirname "$0")/.." && pwd)
cd "$PHYLO"

bash scripts/check_patches.sh || { echo "yamalar eksik, durduruldu" >&2; exit 1; }

mkdir -p "$OUT"
python3 scripts/realize_tree.py --germline "$GERM" --out "$OUT" --plan-only > "$OUT/plan.txt"
cat "$OUT/plan.txt"

# seviye listesi: plan.json'daki parent haritasindan
python3 - "$OUT/plan.json" > "$OUT/levels.txt" <<'PY'
import json, sys
p = json.load(open(sys.argv[1]))
par = p["parent"]
depth = {}
def d(n):
    if n in depth: return depth[n]
    depth[n] = 0 if par[n] is None else d(par[n]) + 1
    return depth[n]
for n in p["nodes"]: d(n)
mx = max(depth.values())
for lv in range(mx + 1):
    print(" ".join(n for n in p["nodes"] if depth[n] == lv))
PY

# --- GOVDE: tek dugum, 2N (WGD ARTIK GOVDEDE DEGIL) ---
# WGD subklonal: agacin icindeki iki dalda (spec copy_number.wgd.events) olusur ve
# realize_tree.py o dugume gelince apply_wgd.py'yi KENDISI cagirir. Burada ozel bir
# adim yok; govde sadece sirali kosulur cunku herkes onu bekliyor.
if [ ! -s "$OUT/MRCA/genome.fa" ]; then
  echo "=== govde: MRCA (2N) ==="
  python3 scripts/realize_tree.py --germline "$GERM" --out "$OUT" --only MRCA
fi
[ -s "$OUT/MRCA/genome.fa" ] || { echo "MRCA basarisiz" >&2; exit 2; }
echo "=== govde tamam, ploidi 2N (46 kontig) ==="
# Cogu dugum 2N (6 Gb -> ~16 GB RSS); WGD kladlarindakiler 4N (12 Gb -> ~32 GB RSS).
# 157 GB RAM'de karisik seviyeler icin guvenli sinir 3.
NP_POST=${NP_POST:-3}

lv=0
while read -r line; do
  [ -z "$line" ] && continue
  lv=$((lv+1))
  echo "=== seviye $lv: $line ==="
  # govde dugumleri yukarida sirali kosuldu, seviye listesinden cikar
  line=$(echo "$line" | tr ' ' '\n' | grep -vE '^MRCA$' | tr '\n' ' ')
  [ -z "$(echo $line)" ] && { echo "  (govde seviyesi, atlandi)"; continue; }
  printf '%s\n' $line | xargs -P "$NP_POST" -I{} \
    python3 scripts/realize_tree.py --germline "$GERM" --out "$OUT" --only {}
  # seviye sonu kontrolu: hepsi uretildi mi?
  miss=""
  for n in $line; do
    [ -s "$OUT/$n/genome.fa" ] || miss="$miss $n"
  done
  if [ -n "$miss" ]; then
    echo "!! seviye $lv EKSIK:$miss -- duruyor (cocuklari kurulamaz)" >&2
    exit 2
  fi
  echo "=== seviye $lv tamam ($(date '+%H:%M:%S')) ==="
done < "$OUT/levels.txt"
echo "TUM AGAC TAMAM $(date '+%H:%M:%S')"
