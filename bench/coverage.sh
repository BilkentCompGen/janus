#!/bin/bash
# Fig. 2b and the depth and purity tables: janus at 5-30x on every pair and on tumor/normal mixtures at 30x, then Severus,
# SAVANA and nanomonsv on the same reads, every call set scored against the consensus of the pair.
#   bench/coverage.sh depth [pair ...]     janus at 5, 10, 15, 20, 30x of both samples (runs/coverage/<pair>_<n>x)
#   bench/coverage.sh purity [pair ...]    janus on tumors of 20-100% purity at 30x (runs/purity/<pair>_<q>p)
#   bench/coverage.sh callers [pair ...]   Severus, SAVANA and nanomonsv on the same read subsets (runs/coverage_callers)
#   bench/coverage.sh score [pair ...]     every call set at every depth and purity against the consensus (runs/coverage_eval/scores.tsv)
# The subsets are taken from the packed reads shuffled once (janus pack --shuffle) with the same order for every depth.
source "$(dirname "$0")/env.sh"
X=${1:?depth, purity, callers or score}; shift
mixes() {   # the shuffled packed reads of a pair under $1
    mkdir -p $1
    for s in n t; do [ -s $1/$s.svb.idx ] || taskset -c $CPUS $J pack --reads $(src $2 $s) --out $1/$s.svb --threads $TH --shuffle > $1/pack_$s.log 2>&1; done
}
case $X in
depth)
    E=$R/coverage
    for p in ${*:-$LABELED}; do
        M=$E/mix/$p; mixes $M $p
        for n in 5 10 15 20 30; do
            [ $n -le $(depth $p) ] || continue
            D=$E/${p}_${n}x; [ -s $D/run.json ] && continue; mkdir -p $D/work
            for s in n t; do $TOOLS mix $D/work/$s.svb $M/$s.svb:0:$((n * G)) >> $Q 2>&1 || continue 2; done
            sheet $D/samples.tsv $D/work/n.svb $D/work/t.svb
            janus_run coverage $D $p $D/samples.tsv $n
        done
        discard $M
    done ;;
purity)
    E=$R/purity
    for p in ${*:-$LABELED}; do
        M=$E/mix/$p; mixes $M $p
        nb=$(python3 -c "import struct; print(struct.unpack_from('<Q', open('$M/n.svb', 'rb').read(32), 16)[0])")
        nm=$((nb - 24 * G - G / 10)); [ $nm -gt $((30 * G)) ] && nm=$((30 * G))
        [ -s $M/nm.svb.idx ] || $TOOLS mix $M/nm.svb $M/n.svb:0:$nm >> $Q 2>&1 || continue
        for q in 20 40 60 80 100; do
            [ -s $M/t$q.svb.idx ] || $TOOLS mix $M/t$q.svb $M/t.svb:0:$((q * 30 * G / 100)) $([ $q -lt 100 ] && echo $M/n.svb:$nm:$(((100 - q) * 30 * G / 100))) >> $Q 2>&1 || continue
            D=$E/${p}_${q}p; mkdir -p $D/work
            for f in svb svb.idx; do ln -sf $M/nm.$f $D/work/n.$f; ln -sf $M/t$q.$f $D/work/t.$f; done
            sheet $D/samples.tsv $D/work/n.svb $D/work/t.svb
            janus_run purity $D $p $D/samples.tsv ""
        done
        discard $M
    done ;;
callers)
    C=$R/coverage_callers; mkdir -p $C/names $C/lists
    sub() { samtools view -@ 12 -b -N $1 -o $3.tmp $2 && mv $3.tmp $3 && samtools index -@ 12 $3; }
    calls() {
        local set=$1 T=$2 N=$3 t n th
        for t in severus_hp savana nanomonsv; do
            n=$t; [ $t = severus_hp ] && n=severus; th=$TH; [ $t = savana ] && th=8
            local D=$C/$n/$set; [ -s $D/calls.vcf ] && continue; mkdir -p $D
            /usr/bin/time -f "%e %M %x" -o $D/time.txt taskset -c $CPUS bash -c "source $H/env.sh; caller $t $T $N $D/out $th" > $D/stdout.txt 2>&1
            [ -s $D/out/calls.vcf ] && { cp $D/out/calls.vcf $D/out/log.txt $D/; discard $D/out; }
        done
    }
    have() { local t; for t in severus savana nanomonsv; do [ -s $C/$t/$1/calls.vcf ] || return 1; done; }
    for p in ${*:-$CASTLE}; do
        W=$C/work/$p; B=$HIFI/$p; mkdir -p $W/30x; ln -sfn $B/clair3 $W/30x/clair3
        for s in n t; do [ -s $C/names/$p.$s.tsv ] || samtools view -@ 8 -F 0x900 $B/$s.bam | awk -F'\t' '{print $1 "\t" length($10)}' > $C/names/$p.$s.tsv; done
        $TOOLS select $p $C "$(echo $CASTLE | tr ' ' ,)"
        for x in 30 20 15 10 5; do
            [ -s $C/lists/$p.n.${x}x.txt ] || continue
            have ${p}_${x}x && continue
            [ -s $W/30x/t.bam.bai ] || sub $C/lists/$p.t.30x.txt $B/t.bam $W/30x/t.bam
            [ -s $W/30x/n.bam.bai ] || sub $C/lists/$p.n.30x.txt $B/n.bam $W/30x/n.bam
            if [ $x != 30 ]; then
                mkdir -p $W/${x}x; ln -sfn $B/clair3 $W/${x}x/clair3
                for s in n t; do [ -s $W/${x}x/$s.bam.bai ] || sub $C/lists/$p.$s.${x}x.txt $W/30x/$s.bam $W/${x}x/$s.bam; done
            fi
            calls ${p}_${x}x $W/${x}x/t.bam $W/${x}x/n.bam
            [ $x != 30 ] && discard $W/${x}x
        done
        [ -s $C/lists/$p.n.purity.txt ] && mkdir -p $W/pn && ln -sfn $B/clair3 $W/pn/clair3
        [ -s $C/lists/$p.n.purity.txt ] && { [ -s $W/pn/n.bam.bai ] || sub $C/lists/$p.n.purity.txt $B/n.bam $W/pn/n.bam; }
        for q in 100 80 60 40 20; do
            [ -s $C/lists/$p.t.${q}p.tumor_part.txt ] || continue
            have ${p}_${q}p && continue
            [ $q = 100 ] && { [ -s $W/30x/t.bam.bai ] || sub $C/lists/$p.t.30x.txt $B/t.bam $W/30x/t.bam; }
            D=$W/p$q; mkdir -p $D
            if [ $q = 100 ]; then ln -sfn $W/30x/t.bam $D/t.bam; ln -sfn $W/30x/t.bam.bai $D/t.bam.bai
            else samtools merge -@ 8 -l 1 -f -o $D/t.bam <(samtools view -@ 4 -u -N $C/lists/$p.t.${q}p.tumor_part.txt $B/t.bam) <(samtools view -@ 4 -u -N $C/lists/$p.t.${q}p.normal_part.txt $B/n.bam) && samtools index -@ 8 $D/t.bam
            fi
            ln -sfn $W/pn/n.bam $D/n.bam; ln -sfn $W/pn/n.bam.bai $D/n.bam.bai; ln -sfn $B/clair3 $D/clair3
            calls ${p}_${q}p $D/t.bam $D/n.bam
            discard $D
        done
        discard $W
    done ;;
score)
    C=$R/coverage_callers; O=${O:-$R/coverage_eval}; VD=${VD:-$H/verdict}; E=${E:-$H/eval_sets}; mkdir -p $O
    base() { case $1 in hg008) echo $HIFI/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz ;; colo829) echo $HIFI/truth/colo829/truth.vcf ;; *) echo $E/$1/all_with/None_minda_ensemble.vcf ;; esac; }
    one() {
        local d=$O/$2/$3; [ -s $d/minda/results/overall.tsv ] && return; mkdir -p $d
        case $3 in janus) grep -E '^#|	PASS	' $4 > $d/in.vcf ;; *) cp $4 $d/in.vcf ;; esac
        tset $(base $1) $d/in.vcf $d/minda
    }
    for p in ${*:-$LABELED}; do
        for s in ${p}_30x ${p}_20x ${p}_15x ${p}_10x ${p}_5x ${p}_100p ${p}_80p ${p}_60p ${p}_40p ${p}_20p; do
            v=$VD/coverage_matched/$s.vcf; [ ${s: -1} = p ] && v=$VD/purity_matched/$s.vcf
            [ -s $v ] && one $p $s janus $v &
            for t in severus savana nanomonsv; do [ -s $C/$t/$s/calls.vcf ] && one $p $s $t $C/$t/$s/calls.vcf & done
        done
        wait
    done
    for f in $O/*/*/minda/results/overall.tsv; do
        s=$(basename $(dirname $(dirname $(dirname $(dirname $f))))); t=$(basename $(dirname $(dirname $(dirname $f))))
        awk -F'\t' -v s=$s -v t=$t 'NR==2{printf "%s\t%s\t%s\t%.4f\t%.4f\t%.4f\n", s, t, $8, $5, $6, $7}' $f
    done | sort > $O/scores.tsv ;;
*) echo "unknown step $X" >&2; exit 1 ;;
esac
