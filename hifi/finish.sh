#!/usr/bin/env bash
# The two runs whose pack never happened: get.sh called tools/svb_pack, which became `svp pack` when
# the binaries were collapsed into one, so both transfers finished and then died at the pack step.
#
# Verify against ENA's own md5 before trusting a .gz that has been sitting on disk, then pack. A
# mismatch means the transfer is bad however complete it looks, so that file is re-fetched from
# scratch and packed after.
set -uo pipefail
cd /home/yuio/svp/hifi
SVP=/home/yuio/svp
API=https://www.ebi.ac.uk/ena/portal/api/filereport
TRASH=/home/yuio/trash/hifi-gz-$(date +%Y%m%d-%H%M%S); mkdir -p "$TRASH"

#   accession   sample      ENA md5
RUNS="
SRR28305160 hcc1954_n 2b30e6e279cdb54562f6d77161467ac3
SRR31537482 hs578_t   5dec906e48dfee01272469be9a267664
"

redownload() {
    local acc=$1 url
    url=$(curl -sS "$API?accession=$acc&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2 | cut -f2)
    [ -n "$url" ] || { echo "[FATAL] $acc: ENA URL vermedi"; return 1; }
    aria2c -x16 -s16 -k10M --file-allocation=none --console-log-level=warn --summary-interval=900 \
        --max-tries=0 --retry-wait=30 --continue=false -o "$acc.gz.part" "https://$url" || return 1
    mv "$acc.gz.part" "$acc.gz"
}

one() {
    local acc=$1 name=$2 want=$3
    local out=${name%_*}/${name##*_}.svb
    mkdir -p "${name%_*}"
    [ -s "$out" ] && { echo "[$name] .svb zaten var"; return 0; }

    echo "[$name] md5 dogrulaniyor $(date +%T)"
    local got; got=$(md5sum "$acc.gz" | cut -d' ' -f1)
    if [ "$got" != "$want" ]; then
        echo "[$name] MD5 TUTMADI: $got != $want -- sifirdan indiriliyor $(date +%T)"
        mv "$acc.gz" "$TRASH/$acc.bad.gz"
        redownload "$acc" || return 1
        got=$(md5sum "$acc.gz" | cut -d' ' -f1)
        [ "$got" = "$want" ] || { echo "[$name] FATAL: yeniden indirmeden sonra da md5 tutmuyor"; return 1; }
    fi
    echo "[$name] md5 OK ($got), paketleniyor $(date +%T)"

    pigz -dc "$acc.gz" | "$SVP/svp" pack "$out.tmp" 2> "${out%.svb}.pack.log" || return 1
    mv "$out.tmp" "$out"
    mv "$acc.gz" "$TRASH/"
    echo "[$name] bitti $(date +%T)  $(ls -l --block-size=1 "$out" | awk '{printf "%.2f GB",$5/1e9}')  $(tail -1 "${out%.svb}.pack.log")"
}

echo "########## 2 kosu $(date +%T) ##########"
echo "$RUNS" | while read -r acc name md5; do
    [ -z "$acc" ] && continue
    one "$acc" "$name" "$md5" &
done
wait
echo "########## bitti $(date +%T) ##########"
ls -l --block-size=1 */*.svb 2>/dev/null | awk '{printf "  %-16s %6.2f GB\n",$9,$5/1e9}'
