#!/usr/bin/env bash
# fetch_more.sh <jobs file>: more tumor/normal HiFi pairs as fastq.gz, one per sample: hifi/<pair>/<n|t>.fastq.gz.
# Job lines (whitespace separated, # comments):
#   pacb  <pair> <n|t> <url> [url...]    PacBio release unaligned BAMs (downloads.pacbcloud.com), verified by S3 ETag
#   giab  <pair> <n|t> <dir-url> <movie.bam> [movie.bam...]   GIAB unaligned BAMs, verified by the folder's checksums.md5
#   sra   <pair> <n|t> <run accession>   SRA/DDBJ run from NCBI's AWS copy, verified by NCBI's md5
# Sources are downloaded into raw/<pair>/ (resumable through aria2's control files; never resumed from length alone),
# converted (samtools fastq or fastq-dump --stdout, | pigz), checked by read count, then the source goes to TRASH.
# At most DL downloads and CV conversions run at once; each waits for disk first. Rerun to resume; finished samples skip.
set -uo pipefail
H=/home/yuio/snape/hifi; RAW=$H/raw; TK=$H/tools/sratoolkit/bin; cd "$H"
SDL="https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?location=s3.us-east-1&type=src&acc="
TRASH=${TRASH:-/home/yuio/trash/reads-converted}; mkdir -p "$TRASH" "$RAW/tmp"
DL=4; CV=3; RESERVE=$((150 * 1000 * 1000 * 1000))
free() { df --output=avail -B1 "$H" | tail -1; }
log() { echo "$(date '+%F %T') $*"; }
slot() {   # $1 kind, $2 count: holds one of the kind's slots on fd 8 (dl) or 7 (cv) until released
    local k=$1 n=$2 i fd; [ $k = dl ] && fd=8 || fd=7
    while :; do for i in $(seq 0 $((n - 1))); do
        eval "exec $fd> $RAW/tmp/$k.$i.lock"; flock -n $fd && return 0
    done; sleep 20; done
}
unslot() { [ $1 = dl ] && { flock -u 8; exec 8>&-; } || { flock -u 7; exec 7>&-; }; }

etag_ok() {   # $1 file, $2 S3 etag "md5-parts": md5 of the concatenated md5s of its parts (part size a whole MiB count)
    python3 - "$1" "$2" <<'EOF'
import hashlib, math, os, sys
f, tag = sys.argv[1], sys.argv[2].strip('"'); want, n = tag.split('-'); n = int(n); size = os.path.getsize(f)
for mib in (8, 16, 32, 64, 128, 256, 512):
    ps = mib << 20
    if math.ceil(size / ps) != n: continue
    d = b''
    with open(f, 'rb') as h:
        while True:
            b = h.read(ps)
            if not b: break
            d += hashlib.md5(b).digest()
    sys.exit(0 if hashlib.md5(d).hexdigest() == want else 1)
sys.exit(2)
EOF
}

get() {   # $1 file, $2 url, $3 bytes, $4 md5 or "" : resumable aria2c download
    local f=$1 url=$2 bytes=$3 md5=$4 tries=0 ck=()
    [ -s "$f" ] && [ ! -e "$f.aria2" ] && return 0
    [ -n "$md5" ] && ck=(--checksum=md5="$md5")
    [ -e "$f.aria2" ] || while [ "$(free)" -lt $(( bytes + RESERVE )) ]; do log "[$f] waiting for disk"; sleep 600; done
    until aria2c -x16 -s16 -k20M --file-allocation=falloc --console-log-level=warn --summary-interval=60 \
          --show-console-readout=false --max-tries=0 --retry-wait=20 --auto-file-renaming=false "${ck[@]}" \
          -d "$(dirname "$f")" -o "$(basename "$f")" "$url" >> "$RAW/tmp/$(basename "$f").dl.log" 2>&1; do
        tries=$((tries + 1)); log "[$f] aria2c failed ($tries), resuming from the control file"
        [ $tries -ge 30 ] && { log "[FATAL] $f: 30 attempts"; return 1; }; sleep 30
    done
}

bams_to_fastq() {   # $1 pair, $2 sample, $3.. bam files: one fastq.gz, read counts compared, sources to trash
    local p=$1 s=$2 g=$H/$1/$2.fastq.gz n=0 got b; shift 2
    slot cv $CV
    while [ "$(free)" -lt $(( $(du -cB1 "$@" | tail -1 | cut -f1) * 12 / 10 + RESERVE / 3 )) ]; do log "[$p $s] waiting for disk to convert"; sleep 600; done
    log "[$p $s] converting"
    for b in "$@"; do samtools fastq -@ 4 "$b" 2>> "$RAW/$p/$s.convert.log"; done | pigz -p 8 > "$g.tmp"
    unslot cv
    for b in "$@"; do n=$((n + $(samtools view -c -@ 4 "$b"))); done
    got=$(( $(pigz -dc "$g.tmp" | wc -l) / 4 ))
    if [ "$n" -gt 0 ] && [ "$n" = "$got" ]; then
        mv "$g.tmp" "$g"; for b in "$@"; do mv "$b" "$TRASH/"; done; log "[$p $s] done: $got reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $p $s: $n reads in the BAMs, $got in the fastq.gz; BAMs kept"; return 1; fi
}

pacb() {   # $1 pair, $2 sample, $3.. urls
    local p=$1 s=$2 u f tag bytes files=(); shift 2
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"
    for u in "$@"; do
        f=$RAW/$p/$s.$(basename "$u"); files+=("$f")
        tag=$(curl -sSI --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="etag:"{print $2}')
        bytes=$(curl -sSI --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
        [ -e "$f.ok" ] && continue
        slot dl $DL; get "$f" "$u" "$bytes" "" || { unslot dl; return 1; }; unslot dl
        if etag_ok "$f" "$tag"; then touch "$f.ok"; log "[$f] downloaded, ETag ok"
        else log "[FATAL] $f: ETag $tag does not match; moved aside"; mkdir -p "$TRASH/bad"; mv "$f" "$TRASH/bad/"; return 1; fi
    done
    bams_to_fastq "$p" "$s" "${files[@]}" && rm -f "${files[@]/%/.ok}"
}

giab() {   # $1 pair, $2 sample, $3 folder url, $4.. bam names
    local p=$1 s=$2 dir=$3 m f md5 bytes sums files=(); shift 3
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"; sums=$(curl -sS --retry 5 "$dir/checksums.md5")
    for m in "$@"; do
        f=$RAW/$p/$s.$m; files+=("$f")
        md5=$(printf '%s\n' "$sums" | grep -F "$m" | grep -v "\.pbi" | awk '{print $1}' | head -1)
        bytes=$(curl -sSI --retry 5 "$dir/$m" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
        [ -n "$md5" ] && [ -n "$bytes" ] || { log "[FATAL] $p $s $m: no md5/size from GIAB"; return 1; }
        slot dl $DL; get "$f" "$dir/$m" "$bytes" "$md5" || { unslot dl; return 1; }; unslot dl
        log "[$f] downloaded, md5 ok"
    done
    bams_to_fastq "$p" "$s" "${files[@]}"
}

sra() {   # $1 pair, $2 sample, $3 run accession
    local p=$1 s=$2 acc=$3 f g=$H/$1/$2.fastq.gz meta url md5 bytes r w
    [ -s "$g" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"; f=$RAW/$p/$s.$acc.sra
    meta=$(curl -sS --retry 5 "$SDL$acc")
    url=$(printf '%s' "$meta" | grep -o '"link": "[^"]*"' | head -1 | cut -d'"' -f4)
    md5=$(printf '%s' "$meta" | grep -o '"md5": "[^"]*"' | head -1 | cut -d'"' -f4)
    bytes=$(printf '%s' "$meta" | grep -o '"size": [0-9]*' | head -1 | grep -o '[0-9]*$')
    [ -n "$url" ] && [ -n "$md5" ] && [ -n "$bytes" ] || { log "[FATAL] $acc: no link/md5/size from NCBI"; return 1; }
    slot dl $DL; get "$f" "$url" "$bytes" "$md5" || { unslot dl; return 1; }; unslot dl
    log "[$f] downloaded, md5 ok"
    slot cv $CV
    while [ "$(free)" -lt $(( bytes * 13 / 10 + RESERVE / 3 )) ]; do log "[$p $s] waiting for disk to convert"; sleep 600; done
    log "[$p $s] converting"
    "$TK/fastq-dump" --stdout --split-spot "$f" 2> "$RAW/$p/$s.convert.log" | pigz -p 8 > "$g.tmp"
    unslot cv
    r=$(grep -o '^Read [0-9]* spots' "$RAW/$p/$s.convert.log" | grep -o '[0-9]*')
    w=$(grep -o '^Written [0-9]* spots' "$RAW/$p/$s.convert.log" | grep -o '[0-9]*')
    if [ -n "$r" ] && [ "$r" = "$w" ] && pigz -t "$g.tmp"; then
        mv "$g.tmp" "$g"; mv "$f" "$TRASH/"; log "[$p $s] done: $r reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $p $s: conversion check failed (read $r, written $w); $f kept"; return 1; fi
}

log "########## start $1 ##########"
while read -r fn args; do
    [ -z "$fn" ] || [[ $fn == \#* ]] && continue
    $fn $args &
    sleep 3
done < "$1"
wait
log "########## end $1 ##########"
