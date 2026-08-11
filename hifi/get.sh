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
# because multi-connection transfers arrive out of order. So the .gz goes to disk after all. It fits:
# 864 GB of .gz plus 465 GB of .svb is 1.33 TB against 1.8 TB free -- but only because each .gz is
# trashed the moment its .svb lands, and trashing frees nothing on this filesystem, so the total is
# what has to fit, not the peak.
set -uo pipefail
cd /home/yuio/snape/hifi
SNAPE=/home/yuio/snape
API=https://www.ebi.ac.uk/ena/portal/api/filereport
TRASH=/home/yuio/trash/hifi-gz; mkdir -p "$TRASH"

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
SRR31537482 hs578_t     175649373095
SRR31537479 hs578_n      87055294667
"

fetch_pack() {   # $1 = accession, $2 = sample name, $3 = expected bases
    # one folder per sample pair: hcc1954/t.svb, hcc1954/n.svb, ...
    local acc=$1 name=$2 want=$3
    local out=${name%_*}/${name##*_}.svb
    mkdir -p "${name%_*}"
    [ -s "$out" ] && { echo "[$name] zaten var"; return 0; }
    local url
    url=$(curl -sS "$API?accession=$acc&result=read_run&fields=fastq_ftp&format=tsv" | tail -n +2 | cut -f2)
    [ -n "$url" ] || { echo "[FATAL] $acc: ENA URL vermedi"; return 1; }
    echo "[$name] $acc indiriliyor $(date +%T)"
    # --max-tries=0: ENA refuses connections intermittently and a 5-try default silently ends the
    # transfer 10 GB in. --continue picks a partial back up instead of restarting 95 GB.
    aria2c -x16 -s16 -k10M --file-allocation=none --console-log-level=warn --summary-interval=900 \
        --max-tries=0 --retry-wait=30 --continue=true -o "$acc.gz.part" "https://$url" || return 1
    mv "$acc.gz.part" "$acc.gz"
    echo "[$name] paketleniyor $(date +%T)"
    pigz -dc "$acc.gz" | "$SNAPE/snape" pack "$out.tmp" 2> "${out%.svb}.pack.log" || return 1
    local got                                     # "[svb_pack] out.svb: N reads, M bases"
    got=$(awk '/bases/{print $(NF-1)}' "${out%.svb}.pack.log")
    [ "$got" = "$want" ] || echo "[$name] UYARI: $got baz paketlendi, ENA $want diyor"
    mv "$out.tmp" "$out"
    mv "$acc.gz" "$TRASH/"
    echo "[$name] bitti $(date +%T)  $(ls -l --block-size=1 "$out" | awk '{printf "%.2f GB",$5/1e9}')  $got baz"
}

# Four at a time: measured on this link one aria2c gets 35 MB/s and five get 89 MB/s combined, so the
# cap is per-transfer at ENA's end. Four leaves cores for the simulation arms running alongside.
echo "########## 10 kosu, 1846.6 Gbp HiFi, aria2c $(date +%T) ##########"
echo "$RUNS" | while read -r acc name want; do
    [ -z "$acc" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge 4 ]; do wait -n; done
    fetch_pack "$acc" "$name" "$want" &
done
wait

echo "########## bitti $(date +%T) ##########"
ls -l --block-size=1 */*.svb 2>/dev/null | awk '{printf "  %-16s %6.2f GB\n",$9,$5/1e9}'
