#!/usr/bin/env bash
# v6: graph.bin'de seq YOK, index.bin var, bubbles/junctions/pick ucu de indeksten okuyor.
# Karsilastirma: e2e/hcc1954 (v5, seq grafta, uc asama okumalari tariyordu), F1 0.878.
set -u
cd /home/yuio/snape
REF=/home/datasets/references/genomes/hg38/hg38.fa
MMI=/home/yuio/snape/e2e/hg38.map-hifi.mmi
T=16; MM=14; ST=4
D=e2e/v6_hcc1954; mkdir -p $D
cp -f e2e/hcc1954/samples.in.tsv $D/

st() { local tag=$1 out=$2; shift 2
    [ -s "$out" ] && { echo "   [$tag] atlandi"; return 0; }
    /usr/bin/time -v "$@" > "$D/log.$tag" 2> "$D/time.$tag"; local rc=$?
    printf "   [%-9s] rc=%d %10s %8s GB\n" "$tag" "$rc" \
        "$(grep -a 'Elapsed (wall' $D/time.$tag | sed 's/.*): //')" \
        "$(grep -a 'Maximum resident' $D/time.$tag | awk '{printf "%.2f",$NF/1048576}')"
    return $rc
}
al() { local tag=$1 in=$2 bam=$3 tmp=$4
    [ -s "$bam" ] && { echo "   [$tag] atlandi"; return 0; }
    /usr/bin/time -v sh -c "minimap2 -ax map-hifi -t $MM -Y --secondary=no '$MMI' '$in' | samtools sort -@ $ST -m 3G -T '$tmp' -o '$bam.tmp' -" > "$D/log.$tag" 2> "$D/time.$tag"
    samtools quickcheck "$bam.tmp" || return 1
    mv "$bam.tmp" "$bam"; samtools index -@ $T "$bam"
    printf "   [%-9s] rc=0 %10s %8s GB\n" "$tag" \
        "$(grep -a 'Elapsed (wall' $D/time.$tag | sed 's/.*): //')" \
        "$(grep -a 'Maximum resident' $D/time.$tag | awk '{printf "%.2f",$NF/1048576}')"
}

echo "===== v6 hcc1954 $(date '+%H:%M:%S')"
st graph     "$D/graph.bin"    ./snape graph --samples "$D/samples.in.tsv" "$D" --threads $T || exit 1
st bubbles   "$D/drop.sigs"    ./snape bubbles "$D/graph.bin" "$D/samples.tsv" $T            || exit 1
st junctions "$D/ctg.t.fa"     ./snape junctions "$D/graph.bin" "$D/samples.tsv" "$D" $T     || exit 1
al mm2_ctg   "$D/ctg.t.fa"     "$D/ctg.t.bam" "$D/stc.t"                                     || exit 1
st call_cand "$D/candidates.vcf" ./snape call "t=$D/ctg.t.bam" - "$D/candidates.vcf" 10 0 50 50 4 || exit 1
st pick      "$D/normal_at_candidates.fastq" ./snape pick "$D/ctg.t.fa" "$D/candidates.vcf" \
             hifi/hcc1954/n.svb "$D/normal_at_candidates.fastq" $T                           || exit 1
al mm2_norm  "$D/normal_at_candidates.fastq" "$D/normal.bam" "$D/stn"                        || exit 1
st call      "$D/calls.vcf"    ./snape call "t=$D/ctg.t.bam" "$D/normal.bam" "$D/calls.vcf" \
             10 0 50 50 5 6 300 0.01 0.40 0.50 0.25 2                                        || exit 1

echo "-- skor"
/home/yuio/.venv/bin/python e2e/score.py "$D/calls.vcf" truth/hcc1954.truth.vcf
echo "   v5 referansi: cagri 1066  P 0.894  R 0.862  F1 0.878"
ls -l --block-size=1 $D/graph.bin $D/index.bin e2e/hcc1954/graph.bin 2>/dev/null | awk '{printf "   %-32s %.2f GB\n",$9,$5/1e9}'
echo "$(date '+%Y-%m-%d %H:%M:%S')  v6 hcc1954 bitti" >> /home/yuio/snape/activity_log.md
