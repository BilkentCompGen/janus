#!/usr/bin/env bash
# fetch_reads.sh: the original HiFi reads of the five benchmark pairs as fastq.gz, in /home/yuio/snape/hifi/raw.
#   CASTLE (HCC1954, NCI-H2009, HCC1937, NCI-H1437): the SRA runs from NCBI's AWS copy (ENA gave 0 MB/s on 2026-09-30,
#   NCBI's S3 74-87 MB/s per run), NCBI md5 checked by aria2c, then fastq-dump --stdout | pigz, spot counts compared.
#   HG008: the BCM Revio unaligned BAMs from GIAB (N-D 2 movies, T 2 movies; 99 + 105 GB, the read sets our svbs were
#   packed from), md5 from GIAB's checksum.md5, then samtools fastq | pigz, one fastq.gz per sample.
# Resumable: just rerun. Finished files are skipped; a partial download resumes from its .aria2 control file (never from
# the file length alone); an unfinished conversion restarts; SRA conversions run side by side, HG008's one at a time. Each source file goes to ~/trash/reads-converted once its
# fastq.gz is verified (the trash is on the same disk: empty it to free the space). Downloads wait for disk, conversions
# run one at a time. Progress: raw/fetch.log, aria2c's own lines per file in raw/tmp/<file>.dl.log; speed: status.sh.
set -uo pipefail
OUT=/home/yuio/snape/hifi/raw; mkdir -p "$OUT/tmp"; cd "$OUT"
TK=/home/yuio/snape/hifi/tools/sratoolkit/bin
SDL="https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?location=s3.us-east-1&type=src&acc="
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Revio_20240313
TRASH=/home/yuio/trash/reads-converted; mkdir -p "$TRASH"
RESERVE=$((80 * 1000 * 1000 * 1000))
PAR=6                                    # jobs at once (a job waiting to convert keeps its slot); each aria2c opens 16 connections
free() { df --output=avail -B1 "$OUT" | tail -1; }
log() { echo "$(date '+%F %T') $*"; }

get() {   # $1 name, $2 url, $3 md5, $4 bytes: aria2c into $1, resumable through $1.aria2
    local f=$1 url=$2 md5=$3 bytes=$4 tries=0 rc
    [ -s "$f" ] && [ ! -e "$f.aria2" ] && return 0
    [ -e "$f.aria2" ] || while [ "$(free)" -lt $(( bytes + RESERVE )) ]; do log "[$f] waiting for disk"; sleep 600; done
    while :; do
        log "[$f] downloading $((bytes / 1000000000)) GB"
        aria2c -x16 -s16 -k20M --file-allocation=falloc --console-log-level=warn --summary-interval=60 --show-console-readout=false --max-tries=0 \
               --retry-wait=20 --lowest-speed-limit=200K --auto-file-renaming=false --checksum=md5="$md5" -o "$f" "$url" >> "tmp/$f.dl.log" 2>&1
        rc=$?; [ $rc -eq 0 ] && { log "[$f] downloaded, md5 ok"; return 0; }
        tries=$((tries + 1))
        if [ $rc -eq 32 ]; then log "[$f] md5 mismatch, starting over"; mkdir -p "$TRASH/bad"; mv "$f" "$f.aria2" "$TRASH/bad/" 2>/dev/null
        else log "[$f] aria2c exit $rc, resuming from the control file"; fi
        [ $tries -ge 30 ] && { log "[FATAL] $f: 30 attempts"; return 1; }
        sleep 30
    done
}

convert_lock() { exec 9> "$OUT/tmp/convert.lock"; flock 9; }

sra() {   # $1 accession, $2 sample name (pair_n / pair_t)
    local acc=$1 name=$2 s="$2.$1.sra" g="$2.fastq.gz" meta url md5 bytes r w
    { [ -s "$g" ] || [ -s "/home/yuio/snape/hifi/${name%_*}/${name##*_}.fastq.gz" ]; } && { log "[$name] present"; return 0; }
    meta=$(curl -sS --retry 5 "$SDL$acc")
    url=$(printf '%s' "$meta" | grep -o '"link": "[^"]*"' | head -1 | cut -d'"' -f4)
    md5=$(printf '%s' "$meta" | grep -o '"md5": "[^"]*"' | head -1 | cut -d'"' -f4)
    bytes=$(printf '%s' "$meta" | grep -o '"size": [0-9]*' | head -1 | grep -o '[0-9]*$')
    [ -n "$url" ] && [ -n "$md5" ] && [ -n "$bytes" ] || { log "[FATAL] $acc: no link/md5/size from NCBI"; return 1; }
    get "$s" "$url" "$md5" "$bytes" || return 1
    while [ "$(free)" -lt $(( bytes * 13 / 10 + RESERVE / 2 )) ]; do log "[$name] waiting for disk to convert"; sleep 600; done
    log "[$name] converting"
    # fastq-dump streams (fasterq-dump stages the whole uncompressed fastq, ~400 GB per run, in its temp directory)
    "$TK/fastq-dump" --stdout --split-spot "$s" 2> "$name.convert.log" | pigz -p 8 > "$g.tmp"
    r=$(grep -o '^Read [0-9]* spots' "$name.convert.log" | grep -o '[0-9]*')
    w=$(grep -o '^Written [0-9]* spots' "$name.convert.log" | grep -o '[0-9]*')
    if [ -n "$r" ] && [ "$r" = "$w" ] && pigz -t "$g.tmp"; then
        mv "$g.tmp" "$g"; mv "$s" "$TRASH/"; log "[$name] done: $r reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $name: conversion check failed (read $r, written $w); $s kept"; return 1; fi
}

hg008() {   # $1 n|t, $2 GIAB folder, $3.. movie names: both movies, then one fastq.gz for the sample
    local s=$1 dir=$2 g="hg008_$1.fastq.gz" m f md5 bytes n=0 got
    shift 2
    { [ -s "$g" ] || [ -s "/home/yuio/snape/hifi/hg008/$s.fastq.gz" ]; } && { log "[hg008_$s] present"; return 0; }
    local sums; sums=$(curl -sS --retry 5 "$GIAB/checksum.md5")
    for m in "$@"; do
        f="hg008_$s.$m.hifi_reads.bam"
        md5=$(printf '%s\n' "$sums" | grep "$dir/$m.hifi_reads.bam$" | awk '{print $1}')
        bytes=$(curl -sSI --retry 5 "$GIAB/$dir/$m.hifi_reads.bam" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
        [ -n "$md5" ] && [ -n "$bytes" ] || { log "[FATAL] hg008_$s $m: no md5/size from GIAB"; return 1; }
        get "$f" "$GIAB/$dir/$m.hifi_reads.bam" "$md5" "$bytes" || return 1
    done
    convert_lock
    log "[hg008_$s] converting"
    for m in "$@"; do samtools fastq -@ 4 "hg008_$s.$m.hifi_reads.bam" 2>> "hg008_$s.convert.log"; done | pigz -p 8 > "$g.tmp"
    flock -u 9
    for m in "$@"; do n=$((n + $(samtools view -c -@ 4 "hg008_$s.$m.hifi_reads.bam"))); done
    got=$(( $(pigz -dc "$g.tmp" | wc -l) / 4 ))
    if [ "$n" -gt 0 ] && [ "$n" = "$got" ]; then
        mv "$g.tmp" "$g"; for m in "$@"; do mv "hg008_$s.$m.hifi_reads.bam" "$TRASH/"; done
        log "[hg008_$s] done: $got reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] hg008_$s: $n reads in the BAMs, $got in the fastq.gz; BAMs kept"; return 1; fi
}

log "########## start ##########"
JOBS="
hg008 n HG008-N-D_ubams m84059_240304_183205_s3 m84059_240304_203144_s4
hg008 t HG008-T_ubams m84101_240226_212100_s2 m84101_240226_232018_s3
sra SRR28305163 hcc1954_t
sra SRR28305160 hcc1954_n
sra SRR28305172 h2009_t
sra SRR28305169 h2009_n
sra SRR28305185 hcc1937_t
sra SRR28305182 hcc1937_n
sra SRR28305179 h1437_t
sra SRR28305175 h1437_n
"
while read -r fn args; do
    [ -z "$fn" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge $PAR ]; do wait -n; done
    $fn $args &
    sleep 5
done <<< "$JOBS"
wait
log "########## end ##########"
