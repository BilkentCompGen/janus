#!/usr/bin/env bash
# realign.sh: replaces the BAMs the other callers use (hifi/<pair>/{n,t}.bam). Those were aligned from FASTA dumped out of
# the svbs, so they had no base qualities; these come from the original reads (hifi/raw/<pair>_{n,t}.fastq.gz, written by
# fetch_reads.sh). Same aligner and settings as before: minimap2 -ax map-hifi -Y --secondary=no against
# ref/hg38_parmask.mmi (the same 194 sequences), samtools sort + index.
# Per pair it waits until both fastq.gz are complete, aligns normal and tumor side by side (16 threads each), checks that
# the BAM holds every read of the fastq.gz (primary records = reads counted at conversion), and only then moves the old
# BAMs to ~/trash/bams-from-svb and the new ones into place. Resumable: rerun; replaced pairs are skipped (marker
# hifi/_aln2/<pair>/replaced), an unfinished alignment restarts. Times per sample go to hifi/_aln2/<pair>/steps.txt.
set -uo pipefail
H=/home/yuio/snape/hifi; RAW=$H/raw; A=$H/_aln2; IDX=/home/yuio/snape/ref/hg38_parmask.mmi; T=16
TRASH=/home/yuio/trash/bams-from-svb; mkdir -p "$TRASH"
log() { echo "$(date '+%F %T') $*"; }
reads_in() { grep -o "\[$1\] done: [0-9]* reads" "$RAW/fetch.log" | tail -1 | grep -o '[0-9]* reads' | cut -d' ' -f1; }

one() {   # $1 pair, $2 n|t
    local p=$1 s=$2 d=$A/$1 fq=$RAW/$1_$2.fastq.gz
    [ -s "$d/$s.ok" ] && return 0
    mkdir -p "$d"
    log "[$p $s] aligning $fq"
    /usr/bin/time -f "align_$s %e s %M KB" bash -c "minimap2 -ax map-hifi -t $T -Y --secondary=no $IDX $fq 2> $d/mm2.$s.log \
        | samtools sort -@ 4 -m 2G -T $d/st.$s -o $d/$s.bam.tmp - && mv $d/$s.bam.tmp $d/$s.bam && samtools index -@ 8 $d/$s.bam" \
        2>> "$d/steps.txt" || { log "[FATAL] $p $s: alignment failed"; return 1; }
    local want got; want=$(reads_in "${p}_$s"); got=$(samtools view -c -@ 8 -F 0x900 "$d/$s.bam")
    [ -n "$want" ] && [ "$want" = "$got" ] || { log "[FATAL] $p $s: $got reads in the BAM, $want in the fastq.gz"; return 1; }
    echo "$got" > "$d/$s.ok"; log "[$p $s] aligned: $got reads"
}

log "########## start ##########"
for p in hcc1954 hg008 h2009 hcc1937 h1437; do
    [ -e "$A/$p/replaced" ] && { log "[$p] already replaced"; continue; }
    until [ -s "$RAW/${p}_n.fastq.gz" ] && [ -s "$RAW/${p}_t.fastq.gz" ]; do sleep 600; done
    one $p n & one $p t & wait
    [ -s "$A/$p/n.ok" ] && [ -s "$A/$p/t.ok" ] || { log "[FATAL] $p: not replaced"; continue; }
    st=$(date +%Y%m%d-%H%M%S)
    for s in n t; do
        for x in bam bam.bai; do [ -e "$H/$p/$s.$x" ] && mv "$H/$p/$s.$x" "$TRASH/$p.$s.$x.$st"; done
        mv "$A/$p/$s.bam" "$H/$p/$s.bam"; mv "$A/$p/$s.bam.bai" "$H/$p/$s.bam.bai"
    done
    touch "$A/$p/replaced"; log "[$p] BAMs replaced (old ones in $TRASH)"
done
log "########## end ##########"
