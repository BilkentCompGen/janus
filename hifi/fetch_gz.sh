#!/usr/bin/env bash
# Downloads the original PacBio HiFi fastq.gz files of the CASTLE pairs from ENA and keeps them (no packing, no removal).
# Resumable: rerun the script; finished files are skipped and partial ones resume from aria2c's control file.
# Disk guard: a download starts only when free space exceeds its size plus RESERVE; the file is preallocated (falloc),
# so running downloads have already claimed their space. Waiting downloads poll every 10 minutes.
set -uo pipefail
OUT=/home/yuio/snape/hifi/raw; mkdir -p "$OUT"; cd "$OUT"
API=https://www.ebi.ac.uk/ena/portal/api/filereport
RESERVE=$((60 * 1000 * 1000 * 1000))
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

fetch() {
    local acc=$1 name=$2 f="$2.$1.fastq.gz"
    [ -s "$f" ] && [ ! -e "$f.part.aria2" ] && [ ! -e "$f.part" ] && { echo "[$name] present"; return 0; }
    local meta url md5 bytes
    meta=$(curl -sS --retry 5 "$API?accession=$acc&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes&format=tsv" | tail -n +2)
    url=$(printf '%s' "$meta" | cut -f2); md5=$(printf '%s' "$meta" | cut -f3); bytes=$(printf '%s' "$meta" | cut -f4)
    [ -n "$url" ] && [ -n "$md5" ] && [ -n "$bytes" ] || { echo "[FATAL] $acc: no URL/md5/size from ENA"; return 1; }
    if [ ! -e "$f.part.aria2" ]; then
        while [ $(( $(df --output=avail -B1 "$OUT" | tail -1) )) -lt $(( bytes + RESERVE )) ]; do
            echo "[$name] waiting for disk: need $((bytes / 1000000000)) GB + reserve $(date +%T)"; sleep 600
        done
    fi
    local tries=0 bad=0 rc
    while :; do
        echo "[$name] $acc downloading $(date +%T) $((bytes / 1000000000)) GB md5=$md5"
        aria2c -x8 -s8 -k10M --file-allocation=falloc --console-log-level=warn --summary-interval=900 \
               --max-tries=0 --retry-wait=30 --auto-file-renaming=false --checksum=md5="$md5" \
               -o "$f.part" "https://$url"
        rc=$?
        [ $rc -eq 0 ] && break
        tries=$((tries + 1))
        if [ $rc -eq 32 ]; then
            bad=$((bad + 1)); echo "[$name] md5 mismatch ($bad), restarting from zero $(date +%T)"
            mkdir -p /home/yuio/trash/gz-bad-md5; mv "$f.part" "$f.part.aria2" /home/yuio/trash/gz-bad-md5/
            [ $bad -ge 3 ] && { echo "[FATAL] $name: md5 failed 3 times"; return 1; }
        else
            echo "[$name] aria2c exit $rc, resuming from control file $(date +%T)"
        fi
        [ $tries -ge 20 ] && { echo "[FATAL] $name: 20 attempts failed"; return 1; }
        sleep 60
    done
    mv "$f.part" "$f"
    echo "[$name] done $(date +%T) $(ls -l --block-size=1 "$f" | awk '{printf "%.1f GB", $5/1e9}')"
}

echo "########## start $(date '+%F %T') ##########"
echo "$RUNS" | while read -r acc name; do
    [ -z "$acc" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge 8 ]; do wait -n; done
    fetch "$acc" "$name" &
    sleep 5
done
wait
echo "########## end $(date '+%F %T') ##########"
