#!/usr/bin/env bash
set -uo pipefail
H=/home/yuio/snape/hifi; RAW=$H/raw; T=$H/truth; TK=$H/tools/sratoolkit/bin; cd "$H"; mkdir -p "$RAW/tmp" "$T"
SDL="https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?location=s3.us-east-1&type=src&acc="
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic
PB=https://downloads.pacbcloud.com/public/revio
DL=4; CV=3; RESERVE=$((150 * 1000 * 1000 * 1000))

JOBS="
sra hcc1954 t SRR28305163
sra hcc1954 n SRR28305160
sra h2009 t SRR28305172
sra h2009 n SRR28305169
sra hcc1937 t SRR28305185
sra hcc1937 n SRR28305182
sra h1437 t SRR28305179
sra h1437 n SRR28305175
giab hg008 t $GIAB/HG008/Liss_lab/BCM_Revio_20240313/HG008-T_ubams $GIAB/HG008/Liss_lab/BCM_Revio_20240313/checksum.md5 m84101_240226_212100_s2.hifi_reads.bam m84101_240226_232018_s3.hifi_reads.bam
giab hg008 n $GIAB/HG008/Liss_lab/BCM_Revio_20240313/HG008-N-D_ubams $GIAB/HG008/Liss_lab/BCM_Revio_20240313/checksum.md5 m84059_240304_183205_s3.hifi_reads.bam m84059_240304_203144_s4.hifi_reads.bam
giab hg009 t $GIAB/HG009/NIST/HG009-T_bulk/20250715p16/BCM_Revio_HG009T-p16_20260313 $GIAB/HG009/NIST/HG009-T_bulk/20250715p16/BCM_Revio_HG009T-p16_20260313/checksums.md5 m84113_260203_223550_s3.demux.bc2061--bc2061.bam m84113_260204_003844_s4.demux.bc2061--bc2061.bam
giab hg009 n $GIAB/HG009/NIST/HG009-N_bulk/20250728p4/BCM_Revio_HG009N-WT-p4_20260313 $GIAB/HG009/NIST/HG009-N_bulk/20250728p4/BCM_Revio_HG009N-WT-p4_20260313/checksums.md5 m84059_260206_135527_s1.demux.bc2067--bc2067.bam m84059_260206_155840_s2.demux.bc2067--bc2067.bam
pacb hcc1395 t $PB/2023Q2/HCC1395/HCC1395/m84039_230414_235240_s2.hifi_reads.default.bam $PB/2023Q2/HCC1395/HCC1395/m84039_230415_002321_s3.hifi_reads.default.bam
pacb hcc1395 n $PB/2023Q2/HCC1395/HCC1395-BL/m84039_230412_183900_s2.hifi_reads.bc2002.bam $PB/2023Q2/HCC1395/HCC1395-BL/m84039_230415_005427_s4.hifi_reads.bc2002.bam
pacb colo829 t $PB/2023Q2/COLO829/COLO829/m84039_230312_025934_s1.hifi_reads.bc2026.bam $PB/2023Q2/COLO829/COLO829/m84039_230328_000836_s3.hifi_reads.bc2025.bam
pacb colo829 n $PB/2023Q2/COLO829/COLO829-BL/m84039_230327_230708_s1.hifi_reads.bc2007.bam $PB/2023Q2/COLO829/COLO829-BL/m84039_230327_233730_s2.hifi_reads.bc2007.bam
sra h209 t DRR1078565
sra h209 n DRR1078566
sra h2126 t DRR1078567
sra h2126 n DRR1078568
sras hs578t t SRR31537482 SRR40571718
sras hs578t n SRR31537479 SRR40571717
"

free() { df --output=avail -B1 "$H" | tail -1; }
log() { echo "$(date '+%F %T') $*"; }
slot() {
    local k=$1 n=$2 i fd; [ $k = dl ] && fd=8 || fd=7
    while :; do for i in $(seq 0 $((n - 1))); do eval "exec $fd> $RAW/tmp/$k.$i.lock"; flock -n $fd && return 0; done; sleep 20; done
}
unslot() { [ $1 = dl ] && { flock -u 8; exec 8>&-; } || { flock -u 7; exec 7>&-; }; }

etag_ok() {
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

get() {
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

bams_to_fastq() {
    local p=$1 s=$2 g=$H/$1/$2.fastq.gz n=0 got b; shift 2
    slot cv $CV
    while [ "$(free)" -lt $(( $(du -cB1 "$@" | tail -1 | cut -f1) * 12 / 10 + RESERVE / 3 )) ]; do log "[$p $s] waiting for disk to convert"; sleep 600; done
    log "[$p $s] converting"
    for b in "$@"; do samtools fastq -@ 4 "$b" 2>> "$RAW/$p/$s.convert.log"; done | pigz -p 8 > "$g.tmp"
    unslot cv
    for b in "$@"; do n=$((n + $(samtools view -c -@ 4 "$b"))); done
    got=$(( $(pigz -dc "$g.tmp" | wc -l) / 4 ))
    if [ "$n" -gt 0 ] && [ "$n" = "$got" ]; then
        mv "$g.tmp" "$g"; rm -f -- "$@" "${@/%/.ok}"; log "[$p $s] done: $got reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $p $s: $n reads in the BAMs, $got in the fastq.gz; BAMs kept"; return 1; fi
}

pacb() {
    local p=$1 s=$2 u f tag bytes files=(); shift 2
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"
    for u in "$@"; do
        f=$RAW/$p/$s.$(basename "$u"); files+=("$f"); [ -e "$f.ok" ] && continue
        tag=$(curl -sSI --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="etag:"{print $2}')
        bytes=$(curl -sSI --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
        slot dl $DL; get "$f" "$u" "$bytes" "" || { unslot dl; return 1; }; unslot dl
        etag_ok "$f" "$tag" && { touch "$f.ok"; log "[$f] downloaded, ETag ok"; } || { log "[FATAL] $f: ETag $tag does not match"; return 1; }
    done
    bams_to_fastq "$p" "$s" "${files[@]}"
}

giab() {
    local p=$1 s=$2 dir=$3 sums m f md5 bytes files=(); sums=$(curl -sS --retry 5 "$4"); shift 4
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"
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

sra() {
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
        mv "$g.tmp" "$g"; rm -f -- "$f"; log "[$p $s] done: $r reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $p $s: conversion check failed (read $r, written $w); $f kept"; return 1; fi
}

sras() {
    local p=$1 s=$2 i=0 a parts=(); shift 2
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    for a in "$@"; do i=$((i + 1)); sra "$p" "${s}_part$i" "$a" || return 1; parts+=("$H/$p/${s}_part$i.fastq.gz"); done
    cat "${parts[@]}" > "$H/$p/$s.fastq.gz.tmp" && pigz -t "$H/$p/$s.fastq.gz.tmp" && mv "$H/$p/$s.fastq.gz.tmp" "$H/$p/$s.fastq.gz" || { log "[FATAL] $p $s: joining failed"; return 1; }
    mkdir -p /home/yuio/trash/fetch-parts && mv "${parts[@]}" /home/yuio/trash/fetch-parts/ && log "[$p $s] done: ${#parts[@]} runs joined, parts in ~/trash/fetch-parts"
}
zen() {
    mkdir -p "$T/$2"
    curl -s "https://zenodo.org/api/records/$1" | python3 -c "
import json,sys
for f in json.load(sys.stdin)['files']: print(f['key'], f['checksum'].split(':')[1], f['links']['self'])" |
    while read -r k md5 url; do
        [ -s "$T/$2/$k" ] && [ "$(md5sum "$T/$2/$k" | cut -d' ' -f1)" = "$md5" ] && continue
        aria2c -x16 -s16 -k1M --file-allocation=none --console-log-level=warn --summary-interval=0 --max-tries=0 \
               --retry-wait=10 --auto-file-renaming=false --allow-overwrite=true -d "$T/$2" -o "$k" "$url" > /dev/null
        [ "$(md5sum "$T/$2/$k" | cut -d' ' -f1)" = "$md5" ] && log "[truth $2] $k md5 ok" || log "[FATAL] truth $2 $k md5 mismatch"
    done
}
web() {
    local d=$1 u; shift; mkdir -p "$T/$d"
    for u in "$@"; do [ -s "$T/$d/$(basename "$u")" ] || curl -sSfL --retry 5 -o "$T/$d/$(basename "$u")" "$u" || log "[FATAL] truth $d $u"; done
}
truth() {
    zen 5019338 colo829
    zen 21231495 prcgap
    web hcc1395 $(for n in 3 4 5 15; do echo "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-022-02816-6/MediaObjects/13059_2022_2816_MOESM${n}_ESM.xlsx"; done)
    web hg009 $(for f in HG009T-p16.severus.somatic.vcf.gz HG009T-p16.severus.somatic.vcf.gz.tbi HG009T-p16.breakpoint_clusters.tsv; do
        echo "$GIAB/HG009/NIST/HG009-T_bulk/20250715p16/analysis/GRCh38/BCM_Revio_HG009T-p16_20260313/HG009-T_vs_HG009-N-WT-p4/structural_variants/$f"; done)
    log "[truth] done"
}

{
log "########## start ${*:-all} ##########"
[ $# -eq 0 ] && truth &
while read -r fn p rest; do
    [ -z "$fn" ] && continue
    [ $# -gt 0 ] && [[ " $* " != *" $p "* ]] && continue
    $fn $p $rest &
    sleep 3
done <<< "$JOBS"
wait
log "########## end ##########"
} 2>&1 | tee -a "$H/fetch.log"
