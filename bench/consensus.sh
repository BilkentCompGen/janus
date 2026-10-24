#!/bin/bash
# Fig. 2a, 2c, 2d and the supplementary accuracy tables: the Minda consensus call sets and every call set scored against them.
#   bench/consensus.sh labels     the published CASTLE consensus of the 11 call sets regenerated with Minda (datasets/hifi/truth/castle11)
#   bench/consensus.sh without    the consensus of the 17 call sets other than janus, per pair (eval_sets/<pair>/all_without): the training labels
#   bench/consensus.sh with       the consensus of all 18 call sets and every call set scored against both (eval_sets/<pair>/all_with, ...)
#   bench/consensus.sh loo        every HiFi call set scored against the consensus of the other 17 (eval_loo)
#   bench/consensus.sh external   janus and every HiFi caller on HG008 and COLO829 against their truth sets (runs/external_eval)
#   bench/consensus.sh callsets   the 18 call sets of the five pairs as given to Minda, collected in bench/callsets/<pair>
# `with`, `loo` and `external` need the classified janus calls of classifier.sh; E, LOO, VD and O override the directories (sweep.sh).
source "$(dirname "$0")/env.sh"
X=${1:?see the header}; shift
E=${E:-$H/eval_sets}; VD=${VD:-$H/verdict}
case $X in
labels)
    T=$HIFI/truth/castle11; mkdir -p $T
    for p in ${*:-$CASTLE}; do
        castle $p
        $MINDA ensemble --vcfs $PH $ONT $SR --conditions "[[[caller_names[:4], caller_names[4:8], caller_names[8:]], '>=', 2],'&', [caller_names,'>=', 4]]" --min_size 50 --tolerance 500 --out_dir $T/$p > $T/$p.log 2>&1
    done ;;
without|with)
    one() {
        local p=$1 W=$E/$1 hifi="" t s v; castle $p; mkdir -p $W/in; : > $W/in/SOURCES
        for t in $HIFI_SETS; do
            v=$R/competitors/$t/$p/calls.vcf; [ -s $v ] || v=$R/competitors/${t%_phased}/$p/calls.vcf; echo "$t $v" >> $W/in/SOURCES
            norm $v $t $W/in/$t.vcf; hifi="$hifi $W/in/$t.vcf"
        done
        [ -s $W/all_without/None_minda_ensemble.vcf ] || ens $W/all_without 10 $hifi $ONT $SR
        [ $X = without ] && return
        grep -E '^#|PASS' $VD/$p.vcf > $W/in/janus.vcf
        ens $W/all_with 11 $hifi $W/in/janus.vcf $ONT $SR
        tset $W/all_without/None_minda_ensemble.vcf $W/in/janus.vcf $W/all_without_janus
        for t in sniffles2_mosaic delly; do
            [ -s $R/competitors/$t/$p/calls.vcf ] || continue
            norm $R/competitors/$t/$p/calls.vcf $t $W/in/$t.vcf
            for s in all_with all_without; do tset $W/$s/None_minda_ensemble.vcf $W/in/$t.vcf $W/${s}_$t; done
        done
    }
    for p in ${*:-$CASTLE}; do one $p & done; wait ;;
loo)
    LO=${LOO:-$H/eval_loo}
    one() {
        local p=$1 x t W hifi; castle $p
        for x in $HIFI_SETS janus; do
            W=$LO/$p/$x; hifi=""; mkdir -p $W; [ -s $W/score/results/overall.tsv ] && continue
            for t in $HIFI_SETS janus; do [ $t = $x ] || hifi="$hifi $E/$p/in/$t.vcf"; done
            [ -s $W/cons/None_minda_ensemble.vcf ] || ens $W/cons 10 $hifi $ONT $SR
            tset $W/cons/None_minda_ensemble.vcf $E/$p/in/$x.vcf $W/score
        done
    }
    for p in ${*:-$CASTLE}; do one $p & done; wait ;;
external)
    O=${O:-$R/external_eval}; mkdir -p $O
    declare -A TB=([hg008]=$HIFI/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz [colo829]=$HIFI/truth/colo829/truth.vcf)
    for p in hg008 colo829; do
        for t in janus $HIFI_SETS sniffles2_mosaic; do
            d=$O/$p/$t; mkdir -p $d
            case $t in
                janus) grep -E '^#|	PASS	' $VD/$p.vcf > $d/in.vcf ;;
                *) v=$R/competitors/$t/$p/calls.vcf; [ -s $v ] || v=$R/competitors/${t%_phased}/$p/calls.vcf; [ -s $v ] || continue; norm $v $t $d/in.vcf ;;
            esac
            [ -s $d/minda/results/overall.tsv ] || tset ${TB[$p]} $d/in.vcf $d/minda
            echo -e "$p\t$t\t$(awk -F'\t' 'NR==2{printf "%s\t%s\t%s\t%.3f\t%.3f\t%.3f", $8, $2, $4, $5, $6, $7}' $d/minda/results/overall.tsv)\t$($TOOLS score ${TB[$p]} $d/in.vcf)"
        done
    done > $O/scores.new; mv $O/scores.new $O/scores.tsv; column -t $O/scores.tsv ;;
callsets)
    for p in ${*:-$CASTLE}; do
        castle $p; d=$H/callsets/$p; mkdir -p $d; cp $E/$p/in/janus.vcf $d/janus.vcf
        for t in $HIFI_SETS; do cp $E/$p/in/$t.vcf $d/${t%_phased}.vcf; done
        set -- $ONT; cp $1 $d/severus_ont.vcf; cp $2 $d/savana_ont.vcf; cp $3 $d/nanomonsv_ont.vcf; cp $4 $d/sniffles2_ont.vcf
        set -- $SR; cp $1 $d/svaba.vcf; cp $2 $d/gripss.vcf; cp $3 $d/manta.vcf
        gzip -9 -f $d/*.vcf; echo "$p: $(ls $d | wc -l) call sets"
    done ;;
*) echo "unknown step $X" >&2; exit 1 ;;
esac
