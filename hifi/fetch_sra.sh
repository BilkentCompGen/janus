#!/usr/bin/env bash
# Fetches the CASTLE HiFi runs from NCBI's AWS copy (ENA throttles and refuses connections) and converts each to fastq.gz.
# Per run: aria2c downloads <name>.<acc>.sra with NCBI's md5, fasterq-dump streams it through pigz into <name>.<acc>.fastq.gz,
# the spot counts are compared, and the .sra is deleted. Rerun to resume: finished runs are skipped, partial downloads resume
# from aria2c's control file, an unfinished conversion restarts. Downloads wait for disk; conversions run one at a time.
set -uo pipefail
OUT=/home/yuio/snape/hifi/raw; mkdir -p "$OUT/tmp"; cd "$OUT"
TK=/home/yuio/snape/hifi/tools/sratoolkit/bin
SDL="https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?location=s3.us-east-1&type=src&acc="
RESERVE=$((100 * 1000 * 1000 * 1000))
RUNS="
SRR28305163 hcc1954_t
SRR28305160 hcc1954_n
SRR28305172 h2009_t
SRR28305169 h2009_n
SRR28305185 hcc1937_t
SRR28305182 hcc1937_n
SRR28305179 h1437_t
SRR28305175 h1437_n
"

free() { df --output=avail -B1 "$OUT" | tail -1; }

fetch() {
    local acc=$1 name=$2 s="$2.$1.sra" g="$2.$1.fastq.gz"
    [ -s "$g" ] && { echo "[$name] present"; return 0; }
    local meta url md5 bytes
    meta=$(curl -sS --retry 5 "$SDL$acc")
    url=$(printf '%s' "$meta" | grep -o '"link": "[^"]*"' | head -1 | cut -d'"' -f4)
    md5=$(printf '%s' "$meta" | grep -o '"md5": "[^"]*"' | head -1 | cut -d'"' -f4)
    bytes=$(printf '%s' "$meta" | grep -o '"size": [0-9]*' | head -1 | grep -o '[0-9]*$')
    [ -n "$url" ] && [ -n "$md5" ] && [ -n "$bytes" ] || { echo "[FATAL] $acc: no link/md5/size from NCBI"; return 1; }
    if [ ! -s "$s" ] || [ -e "$s.aria2" ]; then
        [ -e "$s.aria2" ] || while [ "$(free)" -lt $(( bytes + RESERVE )) ]; do echo "[$name] waiting for disk $(date +%T)"; sleep 600; done
        local tries=0
        until aria2c -x16 -s16 -k20M --file-allocation=falloc --console-log-level=warn --summary-interval=600 --max-tries=0 \
                     --retry-wait=30 --auto-file-renaming=false --checksum=md5="$md5" -o "$s" "$url"; do
            tries=$((tries + 1)); echo "[$name] aria2c failed ($tries), retrying $(date +%T)"
            [ $tries -ge 20 ] && { echo "[FATAL] $name: download failed 20 times"; return 1; }
            sleep 60
        done
        echo "[$name] downloaded $(date +%T) $((bytes / 1000000000)) GB"
    fi
    exec 9> "$OUT/tmp/convert.lock"; flock 9
    while [ "$(free)" -lt $(( bytes * 13 / 10 + RESERVE / 2 )) ]; do echo "[$name] waiting for disk to convert $(date +%T)"; sleep 600; done
    echo "[$name] converting $(date +%T)"
    "$TK/fasterq-dump" -Z -e 8 -t "$OUT/tmp" "$s" 2> "$name.convert.log" | pigz -p 8 > "$g.tmp"
    local r w
    r=$(grep -o 'spots read *: *[0-9,]*' "$name.convert.log" | grep -o '[0-9,]*$' | tr -d ,)
    w=$(grep -o 'reads written *: *[0-9,]*' "$name.convert.log" | grep -o '[0-9,]*$' | tr -d ,)
    flock -u 9
    if [ -n "$r" ] && [ "$r" = "$w" ] && pigz -t "$g.tmp"; then
        mv "$g.tmp" "$g"; mkdir -p /home/yuio/trash/sra-converted; mv "$s" /home/yuio/trash/sra-converted/
        echo "[$name] done $(date +%T) $r reads $(ls -l --block-size=1 "$g" | awk '{printf "%.1f GB", $5/1e9}')"
    else
        echo "[FATAL] $name: conversion check failed (read $r, written $w); .sra kept"; return 1
    fi
}

echo "########## start $(date '+%F %T') ##########"
echo "$RUNS" | while read -r acc name; do
    [ -z "$acc" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge 5 ]; do wait -n; done
    fetch "$acc" "$name" &
    sleep 5
done
wait
echo "########## end $(date '+%F %T') ##########"
