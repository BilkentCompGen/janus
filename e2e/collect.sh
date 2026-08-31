#!/usr/bin/env bash
# e2e sonuclarini tek tabloya topla: asama basina sure ve tepe RSS, cagri sayisi, skor.
set -u
cd /home/yuio/snape
PAIRS="hcc1954 h2009 hcc1937 h1437 hs578"
STAGES="graph bubbles junctions mm2_ctg call_cand pick mm2_norm call"

get() {   # $1 = dosya, $2 = wall|rss|cpu
    local f=$1
    [ -s "$f" ] || { echo "-"; return; }
    case $2 in
        wall) grep -a "Elapsed (wall" "$f" | sed 's/.*): //' | head -1;;
        rss)  grep -a "Maximum resident" "$f" | awk '{printf "%.2f", $NF/1048576}' | head -1;;
        cpu)  grep -a "Percent of CPU" "$f" | awk '{print $NF}' | head -1;;
    esac
}

secs() {  # m:ss.xx veya h:mm:ss -> saniye
    awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1 }' <<< "$1"
}

echo "===== asama basina duvar suresi (16 thread) ====="
printf "%-11s" "cift"; for s in $STAGES; do printf "%11s" "$s"; done; printf "%11s\n" "TOPLAM"
for P in $PAIRS; do
    printf "%-11s" "$P"; tot=0
    for s in $STAGES; do
        w=$(get "e2e/$P/time.$s" wall)
        printf "%11s" "$w"
        [ "$w" = "-" ] || tot=$(awk -v a="$tot" -v b="$(secs "$w")" 'BEGIN{print a+b}')
    done
    printf "%11s\n" "$(awk -v t="$tot" 'BEGIN{printf "%d:%02d", int(t/60), t%60}')"
done

echo
echo "===== asama basina tepe RSS (GB) ====="
printf "%-11s" "cift"; for s in $STAGES; do printf "%11s" "$s"; done; printf "%11s\n" "ZIRVE"
for P in $PAIRS; do
    printf "%-11s" "$P"; mx=0
    for s in $STAGES; do
        r=$(get "e2e/$P/time.$s" rss)
        printf "%11s" "$r"
        [ "$r" = "-" ] || mx=$(awk -v a="$mx" -v b="$r" 'BEGIN{print (b>a)?b:a}')
    done
    printf "%11s\n" "$mx"
done

echo
echo "===== cikti buyuklukleri ====="
printf "%-11s %10s %10s %12s %12s %10s\n" "cift" "graph.bin" "ctg.t.fa" "ctg.t.bam" "normal.bam" "cagri"
for P in $PAIRS; do
    g=$(ls -l --block-size=1 e2e/$P/graph.bin 2>/dev/null | awk '{printf "%.2f GB",$5/1e9}')
    c=$(ls -l --block-size=1 e2e/$P/ctg.t.fa 2>/dev/null | awk '{printf "%.2f GB",$5/1e9}')
    cb=$(ls -l --block-size=1 e2e/$P/ctg.t.bam 2>/dev/null | awk '{printf "%.2f GB",$5/1e9}')
    nb=$(ls -l --block-size=1 e2e/$P/normal.bam 2>/dev/null | awk '{printf "%.2f GB",$5/1e9}')
    v=$(grep -vc '^#' e2e/$P/calls.vcf 2>/dev/null || echo -)
    printf "%-11s %10s %10s %12s %12s %10s\n" "$P" "${g:--}" "${c:--}" "${cb:--}" "${nb:--}" "$v"
done

echo
echo "===== skor (kendi eslestiricimiz, Minda degil; min_size 50, tolerance 500) ====="
for P in hcc1954 h2009 hcc1937 h1437; do
    [ -s "e2e/$P/calls.vcf" ] && [ -s "truth/$P.truth.vcf" ] || { echo "$P: atlandi"; continue; }
    /home/yuio/.venv/bin/python e2e/score.py "e2e/$P/calls.vcf" "truth/$P.truth.vcf"
done
echo "hs578: CASTLE'da truth yok, skorlanmadi"
