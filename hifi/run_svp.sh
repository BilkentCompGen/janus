#!/usr/bin/env bash
# svp on the five real HiFi tumour/normal pairs: graph -> bubbles -> junctions -> minimap2 -> call.
#
# STRICTLY sequential. A real pair is ~2.2x the simulation's read volume and the sim graph peaked at
# ~52 GB RSS, so two of these at once would not fit in 130 GB. Every stage is skipped if its output is
# already there, so a kill and a restart pick up where this left off.
#
# Order is smallest-balanced-first: hcc1937 is the cheapest pair whose normal and tumour have
# comparable depth AND that has a CASTLE truthset, so a resource or quality problem shows up there
# before the expensive pairs run. hs578 is last -- its normal is 87 Gbp against a 176 Gbp tumour, so
# germline cores can fall under the min-mult floor in the normal and read as somatic, and it has no
# truthset to catch that with.
set -uo pipefail
cd /home/yuio/svp
SVP=/home/yuio/svp
REF=/home/datasets/references/genomes/hg38/hg38.fa      # chr-prefixed: the CASTLE truth is chr1, not 1
JXROOT=/home/yuio/tmp/jx

stage() {  # stage <label> <output-to-skip-on> <cmd...>
    local label=$1 out=$2; shift 2
    if [ -s "$out" ]; then echo "  [skip] $label: $out var"; return 0; fi
    echo "  [$label] $(date +%T)"
    /usr/bin/time -v -o "$TIMES/$label.txt" taskset -c 0-15 "$@"
    local rc=$?
    [ $rc -eq 0 ] || { echo "  [FAIL] $label rc=$rc"; return $rc; }
    awk '/Elapsed|Maximum resident/{print "        "$0}' "$TIMES/$label.txt"
}

for S in hcc1937 hcc1954 h2009 h1437 hs578; do
    N=$SVP/hifi/$S/n.svb T=$SVP/hifi/$S/t.svb
    OUT=$SVP/outputs/svp/$S JX=$JXROOT/$S TIMES=$JX/times
    mkdir -p "$OUT" "$JX" "$TIMES"
    echo "########## $S  $(date +%T) ##########"
    df -h /home | awk 'NR==2{print "  disk: "$4" bos"}'

    stage graph "$OUT/graph.bin" "$SVP/svp" graph "$N" "$T" "$OUT" --threads 16 --min-mult 3 || continue

    # bubbles and junctions re-parse the reads, and both readers want FASTQ, so .svb is streamed
    stage bubbles "$OUT/drop.sigs" "$SVP/svp" bubbles "$OUT/graph.bin" \
        <("$SVP/svp" dump "$N") <("$SVP/svp" dump "$T") 16 21 22 40 12 20 || continue

    stage junctions "$JX/ctg.fa" "$SVP/svp" junctions "$OUT/graph.bin" \
        <("$SVP/svp" dump "$T") "$JX/ctg.fa" 16 25 5 300 0 || continue

    if [ ! -s "$JX/ctg.bam" ]; then
        echo "  [align] $(date +%T)"
        /usr/bin/time -v -o "$TIMES/align.txt" taskset -c 0-15 minimap2 -ax map-hifi -t 14 -Y \
            --secondary=no -R "@RG\tID:t\tSM:t" "$REF" "$JX/ctg.fa" 2> "$JX/mm2.log" \
            | samtools sort -@ 4 -m 3G -T "$JX/st" -o "$JX/ctg.bam" - && samtools index -@ 8 "$JX/ctg.bam"
        awk '/Elapsed|Maximum resident/{print "        "$0}' "$TIMES/align.txt"
    else
        echo "  [skip] align: ctg.bam var"
    fi

    stage call "$JX/calls.vcf" "$SVP/svp" call "$JX/ctg.bam" "$JX/calls.vcf" || continue

    echo "  [ozet] $S: contig $(/usr/bin/grep -c '^>' "$JX/ctg.fa" 2>/dev/null), \
call $(/usr/bin/grep -vc '^#' "$JX/calls.vcf" 2>/dev/null), \
truth $(/usr/bin/grep -vc '^#' "$SVP/hifi/$S/truth_minda_4_2.vcf" 2>/dev/null || echo yok)  $(date +%T)"
done
echo "########## hepsi bitti $(date +%T) ##########"
