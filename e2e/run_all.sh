#!/usr/bin/env bash
# Bes CASTLE cifti, ucdan uca, her asama /usr/bin/time -v altinda.
# Sirayla ve bos kutuda: es zamanli kosarlarsa ne sure ne tepe RSS karsilastirilabilir olur.
# Asamalar main.cpp:run_pipeline'in birebir ayni cagrilari; ayri ayri kosuluyorlar ki
# her birinin kendi suresi ve tepe RSS'i olculebilsin.
# Her asama ciktisi varsa atlanir, yani kesilirse kaldigi yerden devam eder.
set -u
cd /home/yuio/snape
REF=/home/datasets/references/genomes/hg38/hg38.fa
T=16
MM=14          # minimap2: threads-2, align() ile ayni
ST=4           # samtools sort -@, align() ile ayni
SNAPE=/home/yuio/snape/snape
PAIRS="hcc1954 h2009 hcc1937 h1437 hs578"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S')  $*" | tee -a /home/yuio/snape/activity_log.md; }

stage() {   # $1 = etiket, $2 = cikti dosyasi, geri kalani komut
    local tag=$1 out=$2; shift 2
    if [ -s "$out" ]; then echo "   [$tag] atlandi ($out var)"; return 0; fi
    local tf="$D/time.$tag"
    /usr/bin/time -v "$@" > "$D/log.$tag" 2> "$tf"
    local rc=$?
    local w r
    w=$(grep -a "Elapsed (wall" "$tf" | sed 's/.*): //')
    r=$(grep -a "Maximum resident" "$tf" | awk '{printf "%.2f", $NF/1048576}')
    echo "   [$tag] rc=$rc  $w  ${r} GB"
    [ $rc -eq 0 ] || { log "FATAL $P/$tag rc=$rc"; return 1; }
}

align() {   # $1 = etiket, $2 = girdi (fa/fastq), $3 = cikti bam, $4 = sort tmp on eki
    local tag=$1 in=$2 bam=$3 tmp=$4
    if [ -s "$bam" ]; then echo "   [$tag] atlandi ($bam var)"; return 0; fi
    local tf="$D/time.$tag"
    /usr/bin/time -v sh -c "minimap2 -ax map-hifi -t $MM -Y --secondary=no '$REF' '$in' \
        | samtools sort -@ $ST -m 3G -T '$tmp' -o '$bam.tmp' -" > "$D/log.$tag" 2> "$tf"
    local rc=$?
    samtools quickcheck "$bam.tmp" || rc=1
    [ $rc -eq 0 ] || { log "FATAL $P/$tag rc=$rc"; return 1; }
    mv "$bam.tmp" "$bam"
    samtools index -@ $T "$bam"
    local w r
    w=$(grep -a "Elapsed (wall" "$tf" | sed 's/.*): //')
    r=$(grep -a "Maximum resident" "$tf" | awk '{printf "%.2f", $NF/1048576}')
    echo "   [$tag] rc=0  $w  ${r} GB"
}

log "e2e: bes cift basliyor, T=$T"
for P in $PAIRS; do
    D=e2e/$P
    mkdir -p "$D"
    printf 'name\trole\tpath\n'                                   >  "$D/samples.in.tsv"
    printf 'n\tnormal\t/home/yuio/snape/hifi/%s/n.svb\n' "$P"     >> "$D/samples.in.tsv"
    printf 't\ttumor\t/home/yuio/snape/hifi/%s/t.svb\n'  "$P"     >> "$D/samples.in.tsv"
    log "e2e $P basladi"
    echo "===== $P $(date '+%H:%M:%S')"

    stage graph      "$D/graph.bin"   $SNAPE graph --samples "$D/samples.in.tsv" "$D" --threads $T || continue
    stage bubbles    "$D/drop.sigs"   $SNAPE bubbles "$D/graph.bin" "$D/samples.tsv" $T           || continue
    stage junctions  "$D/ctg.t.fa"    $SNAPE junctions "$D/graph.bin" "$D/samples.tsv" "$D" $T    || continue
    align mm2_ctg    "$D/ctg.t.fa"    "$D/ctg.t.bam" "$D/stc.t"                                   || continue
    stage call_cand  "$D/candidates.vcf" $SNAPE call "t=$D/ctg.t.bam" - "$D/candidates.vcf" 10 0 50 50 4 || continue
    stage pick       "$D/normal_at_candidates.fastq" $SNAPE pick "$D/ctg.t.fa" "$D/candidates.vcf" \
                     "hifi/$P/n.svb" "$D/normal_at_candidates.fastq" $T                           || continue
    align mm2_norm   "$D/normal_at_candidates.fastq" "$D/normal.bam" "$D/stn"                     || continue
    stage call       "$D/calls.vcf"   $SNAPE call "t=$D/ctg.t.bam" "$D/normal.bam" "$D/calls.vcf" \
                     10 0 50 50 5 6 300 0.01 0.40 0.50 0.25 2                                     || continue

    log "e2e $P bitti, $(grep -vc '^#' $D/calls.vcf 2>/dev/null) cagri"
    echo "  -> $(grep -vc '^#' $D/calls.vcf 2>/dev/null) cagri, $(du -sh $D | cut -f1)"
done
log "e2e: bes cift bitti"
