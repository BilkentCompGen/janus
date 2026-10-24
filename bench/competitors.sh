#!/bin/bash
# Fig. 2a, 2c and Fig. 4a: the other callers on the pairs. Their inputs are the reads aligned with minimap2 (align) and,
# for Severus and SAVANA, haplotagged after Clair3 + LongPhase phasing of the normal (phase); every step is measured.
#   bench/competitors.sh align [pair ...]       minimap2 of both read sets of a pair (runs/align/<pair>/<sample>, datasets/hifi/<pair>/<sample>.bam)
#   bench/competitors.sh phase [pair ...]       Clair3 + LongPhase on the normal, haplotag both BAMs (runs/phase/<pair>)
#   bench/competitors.sh call [tool] [pair ...] a caller on the pairs (runs/competitors/<tool>/<pair>)
source "$(dirname "$0")/env.sh"
X=${1:?align, phase or call}; shift
case $X in
align)
    aligned() { minimap2 -ax map-hifi -t $TH -Y --secondary=no -R "@RG\tID:${1}_$2\tSM:${1}_$2" $K/ref/hg38_parmask.mmi $HIFI/$1/$2.fastq.gz 2> $3.minimap2.log | samtools sort -@ 4 -m 2G -T $3.st -o $3.tmp - && mv $3.tmp $3 && samtools index -@ 8 $3; }
    export -f aligned; export TH REF HIFI K
    for p in ${*:-$LABELED}; do for s in n t; do
        measured $R/align/$p/$s "{\"experiment\": \"align\", \"tool\": \"minimap2+samtools\", \"pair\": \"$p\", \"sample\": \"$s\", \"bam\": \"$HIFI/$p/$s.bam\"}" 60G bash -c "aligned $p $s $HIFI/$p/$s.bam"
    done; done ;;
phase)
    IMG=hkubal/clair3:v2.0.3; FA=$REF/hg38_parmask.fa
    for p in ${*:-$CASTLE}; do
        W=$HIFI/$p; O=$R/phase/$p; T=$O/times.tsv; mkdir -p $O
        [ -s $T ] || echo -e "step\twall_s\tcpu_s\tpeak_gib\texit" > $T
        step() {
            local name=$1 sl=sb_phase_${p}_$1.slice; shift
            grep -q "^$name	" $T && return 0
            CG=$sl /usr/bin/time -f "%e %U %S %M %x" -o $O/$name.time "$@" > $O/$name.log 2>&1
            local c=/sys/fs/cgroup/$sl cpu=0 peak=0
            [ -d $c ] && { cpu=$(awk '/^usage_usec/{print $2/1e6}' $c/cpu.stat); peak=$(awk '/^anon /{print $2}' $c/memory.stat); }
            read wall us sy rss ex < $O/$name.time
            awk -v n=$name -v w=$wall -v c=$cpu -v u=$us -v y=$sy -v m=$peak -v r=$rss -v x=$ex 'BEGIN{r *= 1024; if (r > m) m = r; printf "%s\t%s\t%.1f\t%.2f\t%s\n", n, w, c + u + y, m / 1073741824, x}' >> $T
        }
        dkc() { docker run --rm --log-driver none --cpuset-cpus $CPUS --cgroup-parent=$CG -u "$(id -u):$(id -g)" -v $S:$S $IMG "$@"; }
        export -f dkc; export S IMG CPUS
        step clair3 bash -c "dkc /opt/bin/run_clair3.sh --bam_fn=$W/n.bam --ref_fn=$FA --threads=$TH --platform=hifi --model_path=/opt/models/hifi_revio --output=$W/clair3 --enable_phasing --longphase_for_phasing"
        for s in n t; do
            [ -s $W/$s.tagged ] && continue
            step haplotag_$s bash -c "dkc /opt/bin/longphase haplotag -s $W/clair3/phased_merge_output.vcf.gz -b $W/$s.bam -r $FA -t $TH --tagSupplementary -o $W/$s.hp"
            samtools view -@ 8 -b -o $W/$s.unmapped.bam $W/$s.bam '*'
            samtools cat -o $W/$s.hpu.bam $W/$s.hp.bam $W/$s.unmapped.bam && mv $W/$s.hpu.bam $W/$s.hp.bam
            step index_$s taskset -c $CPUS samtools index -@ $TH $W/$s.hp.bam
            a=$(samtools idxstats $W/$s.bam | awk '{m+=$3; u+=$4} END{print m, u}'); b=$(samtools idxstats $W/$s.hp.bam | awk '{m+=$3; u+=$4} END{print m, u}')
            [ "$a" = "$b" ] || { log "$p $s: read counts differ ($a vs $b); original kept"; continue; }
            mv $W/$s.hp.bam $W/$s.bam; mv $W/$s.hp.bam.bai $W/$s.bam.bai; discard $W/$s.unmapped.bam; echo "$a" > $W/$s.tagged
        done
    done ;;
call)
    for t in ${1:-$TOOLS_ALL}; do for p in ${2:-$CASTLE}; do
        mem=60G; th=$TH; [ $t = savana ] && { mem=150G; th=8; }; [ $t = severus_hp ] && [ ! -s $HIFI/$p/n.tagged ] && continue
        n=$t; [ $t = severus_hp ] && n=severus_phased
        measured $R/competitors/$n/$p "{\"experiment\": \"competitor\", \"tool\": \"$n\", \"pair\": \"$p\", \"threads\": $th}" $mem \
            bash -c "source $H/env.sh; caller $t $HIFI/$p/t.bam $HIFI/$p/n.bam $R/competitors/$n/$p/out $th"
    done; done ;;
*) echo "unknown step $X" >&2; exit 1 ;;
esac
