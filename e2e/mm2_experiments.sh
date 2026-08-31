#!/usr/bin/env bash
# minimap2 boru hattinin %52.3'u. Uc yon, hepsi ayni girdide, bos kutuda.
#
# Baslangic noktasi (docs/NOTES.md, hcc1954, 16 thread):
#   contig  2.04 Gbp / 4:00  = 8.50 Mbp/s, ortalama sorgu 9,173 bp
#   normal  3.79 Gbp / 16:25 = 3.85 Mbp/s, ortalama sorgu 17,814 bp
# Baz basina 2.21 kat fark, uzunluk orani 1.94 -- maliyet uzunlukta superlineer.
set -u
cd /home/yuio/snape
REF=/home/datasets/references/genomes/hg38/hg38.fa
MMI=/home/yuio/snape/e2e/hg38.map-hifi.mmi
D=e2e/mm2exp; mkdir -p $D
CTG=e2e/hcc1954/ctg.t.fa
NFQ=e2e/hcc1954/normal_at_candidates.fastq

log() { echo "$(date '+%Y-%m-%d %H:%M:%S')  $*" | tee -a /home/yuio/snape/activity_log.md; }

run() {   # $1 etiket, geri kalani komut
    local tag=$1; shift
    [ -s "$D/$tag.time" ] && { printf "   %-24s (atlandi)\n" "$tag"; return; }
    /usr/bin/time -v "$@" > $D/$tag.out 2> $D/$tag.time
    local w r c
    w=$(grep -a "Elapsed (wall" $D/$tag.time | sed 's/.*): //')
    r=$(grep -a "Maximum resident" $D/$tag.time | awk '{printf "%.2f", $NF/1048576}')
    c=$(grep -a "Percent of CPU" $D/$tag.time | awk '{print $NF}')
    printf "   %-24s %10s  %7s GB  %6s\n" "$tag" "$w" "$r" "$c"
}

echo "===== minimap2 deneyleri $(date '+%H:%M:%S'), yuk: $(uptime | sed 's/.*average: //')"
log "mm2 deneyleri basladi"
printf "   %-24s %10s  %7s  %6s\n" "deney" "duvar" "RSS" "CPU"

echo "-- 1. FASTA indeks vs onceden kurulmus .mmi (ayni cikti mi, ne kadar kazandiriyor)"
run ctg_fa_t14  sh -c "minimap2 -ax map-hifi -t 14 -Y --secondary=no '$REF' '$CTG' > $D/ctg_fa.sam"
run ctg_mmi_t14 sh -c "minimap2 -ax map-hifi -t 14 -Y --secondary=no '$MMI' '$CTG' > $D/ctg_mmi.sam"
if [ -s $D/ctg_fa.sam ] && [ -s $D/ctg_mmi.sam ]; then
    if cmp -s <(grep -av "^@PG" $D/ctg_fa.sam) <(grep -av "^@PG" $D/ctg_mmi.sam)
        then echo "   ciktilar ayni (@PG disinda)"; else echo "   CIKTILAR FARKLI"; fi
    mv -f $D/ctg_fa.sam $D/ctg_mmi.sam /home/yuio/trash/ 2>/dev/null
fi

echo "-- 2. thread taramasi (contig'ler, .mmi)"
for t in 8 14 22 30; do
    run ctg_t$t sh -c "minimap2 -ax map-hifi -t $t -Y --secondary=no '$MMI' '$CTG' > /dev/null"
done

echo "-- 3. bayraklar (contig'ler, -t 30)"
run ctg_t30_K2G   sh -c "minimap2 -ax map-hifi -t 30 -K 2G -Y --secondary=no '$MMI' '$CTG' > /dev/null"
run ctg_t30_2io   sh -c "minimap2 -ax map-hifi -t 30 -2 -Y --secondary=no '$MMI' '$CTG' > /dev/null"
run ctg_t30_chain sh -c "minimap2 -ax map-hifi -t 30 --max-chain-skip 10 --max-chain-iter 1000 -Y --secondary=no '$MMI' '$CTG' > /dev/null"

echo "-- 4. normal okumalar: kirpilmamis vs kirpilmis (asil kaldirac)"
if [ -s "$NFQ" ]; then
    run norm_t14      sh -c "minimap2 -ax map-hifi -t 14 -Y --secondary=no '$MMI' '$NFQ' > /dev/null"
    run norm_t30      sh -c "minimap2 -ax map-hifi -t 30 -Y --secondary=no '$MMI' '$NFQ' > /dev/null"
    for tr in 1500 3000; do
        f=$D/norm_trim$tr.fastq
        if [ ! -s "$f" ]; then
            ./snape pick e2e/hcc1954/ctg.t.fa e2e/hcc1954/candidates.vcf \
                hifi/hcc1954/n.svb "$f" 16 5 $tr > $D/pick_trim$tr.out 2> $D/pick_trim$tr.time
        fi
        grep -a "\[pick\]" $D/pick_trim$tr.time | sed 's/^/   /'
        run norm_trim${tr}_t30 sh -c "minimap2 -ax map-hifi -t 30 -Y --secondary=no '$MMI' '$f' > /dev/null"
    done
else
    echo "   $NFQ yok, atlandi"
fi
log "mm2 deneyleri bitti"
