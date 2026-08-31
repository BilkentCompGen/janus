#!/usr/bin/env bash
# minimap2'nin en buyuk kalemi: secilmis normal okumalarin hizalanmasi (butcede 16:25, %42.1).
# pick okumalari BUTUN yaziyor (ortalama 17,814 bp) ama germline testi cagri konumunun
# +-500 bp'sine bakiyor. Kirpma ne kadar kazandiriyor, ve dogruluga ne yapiyor.
#
# e2e arka planda kosuyor (16 thread). Bu betik -t 8 kullaniyor ve bos 16 cekirdege siginiyor;
# butun kollar ayni yuk altinda oldugu icin ORANLAR gecerli, mutlak sureler degil.
set -u
cd /home/yuio/snape
MMI=/home/yuio/snape/e2e/hg38.map-hifi.mmi
D=e2e/mm2exp; mkdir -p $D
T=8

run() { local tag=$1; shift
    [ -s "$D/$tag.time" ] && { printf "   %-22s (atlandi)\n" "$tag"; return; }
    nice -n 5 /usr/bin/time -v "$@" > $D/$tag.out 2> $D/$tag.time
    printf "   %-22s %10s  %6s GB  %6s\n" "$tag" \
        "$(grep -a 'Elapsed (wall' $D/$tag.time | sed 's/.*): //')" \
        "$(grep -a 'Maximum resident' $D/$tag.time | awk '{printf "%.2f",$NF/1048576}')" \
        "$(grep -a 'Percent of CPU' $D/$tag.time | awk '{print $NF}')"
}

echo "===== kirpma deneyi $(date '+%H:%M:%S'), -t $T, yuk $(uptime | sed 's/.*average: //')"
printf "   %-22s %10s  %6s  %6s\n" "kol" "duvar" "RSS" "CPU"

NFQ=e2e/hcc1954/normal_at_candidates.fastq
awk 'NR%4==2{n++; b+=length($0)} END{printf "   kirpilmamis: %d okuma, %.3f Gbp, ortalama %.0f bp\n",n,b/1e9,b/n}' "$NFQ"
run norm_full sh -c "minimap2 -ax map-hifi -t $T -Y --secondary=no '$MMI' '$NFQ' | samtools sort -@ 2 -m 2G -T $D/s0 -o $D/norm_full.bam -"

for tr in 3000 1500 500; do
    f=$D/norm_trim$tr.fastq
    if [ ! -s "$f" ]; then
        nice -n 5 ./snape pick e2e/hcc1954/ctg.t.fa e2e/hcc1954/candidates.vcf \
            hifi/hcc1954/n.svb "$f" $T 5 $tr > $D/pick$tr.out 2> $D/pick$tr.time
    fi
    grep -a "\[pick\]" $D/pick$tr.time | sed 's/^/   /'
    run norm_trim$tr sh -c "minimap2 -ax map-hifi -t $T -Y --secondary=no '$MMI' '$f' | samtools sort -@ 2 -m 2G -T $D/s$tr -o $D/norm_trim$tr.bam -"
done

echo
echo "-- dogruluk: her hizalamayla call'i tekrar kosup F1"
for v in full trim3000 trim1500 trim500; do
    b=$D/norm_$v.bam
    [ -s "$b" ] || continue
    [ -s "$b.bai" ] || samtools index -@ 4 "$b"
    [ -s "$D/calls_$v.vcf" ] || nice -n 5 ./snape call "t=e2e/hcc1954/ctg.t.bam" "$b" "$D/calls_$v.vcf" \
        10 0 50 50 5 6 300 0.01 0.40 0.50 0.25 2 > $D/call_$v.out 2> $D/call_$v.time
    printf "   %-10s " "$v"
    /home/yuio/.venv/bin/python e2e/score.py "$D/calls_$v.vcf" truth/hcc1954.truth.vcf
done
