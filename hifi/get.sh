#!/usr/bin/env bash
# PacBio Revio HiFi tumour/normal pairs from PRJNA1086849 -- the Severus paper's own BioProject.
#
# WHY HiFi and not the ONT R10 that came before it: the COLO829 R10 graph measured a core
# REPRODUCIBILITY of 6%. A germline core with 116x of tumour over it carried a median tcov of 7,
# because an L5 core spans ~217 bp and R10's ~1% error rate rarely leaves that many bases untouched.
# Signal (7) landed on top of the read-error floor (3), so precision came out at 0.003 and no
# threshold -- min-mult or VAF ratio -- could separate them. The HiFi simulation measures 97% on the
# same statistic, and these runs sample at Q39.7-Q40.0, ~0.01% error.
#
# The archive names every PacBio file "_subreads.fastq.gz" regardless of content. These are NOT CLR
# subreads: sampled reads are 6.7-19 kb at Q35-40, which is CCS/HiFi. Revio is a HiFi-only instrument.
#
# COLO829 and HCC1395 are absent here on purpose -- neither has PacBio in this BioProject, only
# Illumina and ONT.
#
# The first version of this script streamed curl -> pigz -> svb_pack so no .gz ever hit disk. It was
# measured at 16 Mbase/s, a 30-hour run: ENA hands a single connection ~2 MB/s, and curl opens one.
# aria2c's -x16 is what got 35 MB/s per run on the COLO829 downloads, and it cannot write to a pipe
# because multi-connection transfers arrive out of order. So the .gz goes to disk after all.
#
# md5, 2026-08-30: an earlier pass of this script finished four runs at EXACTLY the ENA byte count
# and all four unpacked to garbage. Cause is the combination that was here before -- -x16 writes
# pieces out of order, --file-allocation=none makes the file's size the highest offset written, and
# --continue=true on a restart whose .aria2 control file was gone reads that size as "already have
# this much" and appends the remainder. Right length, wrong bytes, no error. Two changes: ENA's
# fastq_md5 goes to aria2c --checksum so a bad transfer fails loudly at the end of the transfer
# instead of silently at the science, and --continue is gone so resume happens only through the
# control file, which is the only mechanism that knows which pieces are real. Exit 32 is aria2c's
# checksum failure: that one restarts from zero, everything else resumes.
set -uo pipefail
cd /home/yuio/snape/hifi
SNAPE=/home/yuio/snape
API=https://www.ebi.ac.uk/ena/portal/api/filereport
TRASH=/home/yuio/trash/hifi-gz; mkdir -p "$TRASH"
LOG=$SNAPE/file-removals.txt

# accession sample expected_bases
RUNS="
SRR28305163 hcc1954_t   196448478603
SRR28305160 hcc1954_n   195080358066
SRR28305172 h2009_t     201768619118
SRR28305169 h2009_n     211158885701
SRR28305185 hcc1937_t   185423384135
SRR28305182 hcc1937_n   173700624976
SRR28305179 h1437_t     199756958598
SRR28305175 h1437_n     220468470134
"

trash_partial() {   # $1 = accession, $2 = why
    local acc=$1 why=$2 stamp f
    stamp=$(date +%Y%m%d-%H%M%S)
    for f in "$acc.gz.part" "$acc.gz.part.aria2"; do
        [ -e "$f" ] || continue
        mv "$f" "$TRASH/$stamp-$f"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  mv $f $TRASH/$stamp-$f  ($why)" >> "$LOG"
    done
}

fetch_pack() {   # $1 = accession, $2 = sample name, $3 = expected bases
    # one folder per sample pair: hcc1954/t.svb, hcc1954/n.svb, ...
    local acc=$1 name=$2 want=$3
    local out=${name%_*}/${name##*_}.svb
    mkdir -p "${name%_*}"
    [ -s "$out" ] && { echo "[$name] zaten var"; return 0; }

    local meta url md5
    meta=$(curl -sS "$API?accession=$acc&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n +2)
    url=$(printf '%s' "$meta" | cut -f2)
    md5=$(printf '%s' "$meta" | cut -f3)
    [ -n "$url" ] && [ -n "$md5" ] || { echo "[FATAL] $acc: ENA URL/md5 vermedi"; return 1; }

    local tries=0 bad=0 rc
    while :; do
        echo "[$name] $acc indiriliyor $(date +%T)  md5=$md5"
        # --max-tries=0: ENA refuses connections intermittently and a 5-try default silently ends the
        # transfer 10 GB in. No --continue: only the .aria2 control file may resume this.
        aria2c -x16 -s16 -k10M --file-allocation=none --console-log-level=warn \
               --summary-interval=900 --max-tries=0 --retry-wait=30 \
               --auto-file-renaming=false --checksum=md5="$md5" \
               -o "$acc.gz.part" "https://$url"
        rc=$?
        [ $rc -eq 0 ] && break
        tries=$((tries + 1))
        if [ $rc -eq 32 ]; then                   # aria2c: checksum validation failed
            bad=$((bad + 1))
            echo "[$name] md5 TUTMADI ($bad. kez), sifirdan $(date +%T)"
            trash_partial "$acc" "aria2c md5 dogrulamasi basarisiz, yeniden indiriliyor"
            [ $bad -ge 3 ] && { echo "[FATAL] $name: uc denemede de md5 tutmadi"; return 1; }
        else
            echo "[$name] aria2c $rc ile dustu, kontrol dosyasindan devam $(date +%T)"
        fi
        [ $tries -ge 8 ] && { echo "[FATAL] $name: 8 denemede inmedi"; return 1; }
        sleep 60
    done
    mv "$acc.gz.part" "$acc.gz"

    echo "[$name] paketleniyor $(date +%T)"
    pigz -dc "$acc.gz" | "$SNAPE/snape" pack "$out.tmp" 2> "${out%.svb}.pack.log" || return 1
    local got                                     # "[svb_pack] out.svb: N reads, M bases"
    got=$(awk '/bases/{print $(NF-1)}' "${out%.svb}.pack.log")
    [ "$got" = "$want" ] || echo "[$name] UYARI: $got baz paketlendi, ENA $want diyor"
    mv "$out.tmp" "$out"
    mv "$acc.gz" "$TRASH/"
    echo "$(date '+%Y-%m-%d %H:%M:%S')  mv $acc.gz $TRASH/  (svb paketlendi, ham .gz artik gereksiz)" >> "$LOG"
    echo "[$name] bitti $(date +%T)  $(ls -l --block-size=1 "$out" | awk '{printf "%.2f GB",$5/1e9}')  $got baz"
}

# Four at a time: measured on this link one aria2c gets 35 MB/s and five get 89 MB/s combined, so the
# cap is per-transfer at ENA's end.
echo "########## 10 kosu, 1846.6 Gbp HiFi, 864.3 GB .gz, aria2c $(date +%T) ##########"
echo "$RUNS" | while read -r acc name want; do
    [ -z "$acc" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge 4 ]; do wait -n; done
    fetch_pack "$acc" "$name" "$want" &
done
wait

echo "########## bitti $(date +%T) ##########"
ls -l --block-size=1 */*.svb 2>/dev/null | awk '{printf "  %-16s %6.2f GB\n",$9,$5/1e9}'
