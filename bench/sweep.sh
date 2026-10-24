#!/bin/bash
# Supplementary figure of the count-stage thresholds: tau_normal in (5, 10, 15, 20, 25) x tau_tumor in (20, 30, 40, 50, 60),
# plus (70, 40), (10, 5) and (70, 5) so that every combination of minimum read supports (m0, mt) occurs, on the CASTLE full,
# depth and purity runs. janus is rerun only for a pair of minimum read supports not seen before; every distinct
# configuration is classified (each pair held out) and scored like the paper's call sets.
#   bench/sweep.sh runs [pair ...]   the reruns (runs/sweep/runs/<experiment>/<run>/t<tau_normal>_<tau_tumor>)
#   bench/sweep.sh score             the configurations: classifier, consensus with and without janus, leave-one-out, every depth and purity (runs/sweep/conf/<config>, rescore.txt)
source "$(dirname "$0")/env.sh"
X=${1:?runs or score}; shift
W=$R/sweep; SW="$TOOLS sweep"
case $X in
runs)
    mkdir -p $W/runs
    one() {
        local k=$1 n=$2 a=$3 b=$4 p=${2%%_*} D=$W/runs/$1/$2/t$3_$4 MX=$W/mix/${2%%_*} x
        mkdir -p $D/work; log "[$k/$n tau $a/$b] $([ -e $D/work/done.extract_t ] && echo resume || echo start)"
        case $k in
            full) sheet $D/samples.tsv $(src $p n) $(src $p t) ;;
            coverage) x=${n##*_}; x=${x%x}; for s in n t; do [ -s $D/work/$s.svb.idx ] || $TOOLS mix $D/work/$s.svb $MX/$s.svb:0:$((x * G)) >> $Q 2>&1; done; sheet $D/samples.tsv $D/work/n.svb $D/work/t.svb ;;
            purity) x=${n##*_}; x=${x%p}; for f in svb svb.idx; do ln -sf $MX/nm.$f $D/work/n.$f; ln -sf $MX/t$x.$f $D/work/t.$f; done; sheet $D/samples.tsv $D/work/n.svb $D/work/t.svb ;;
        esac
        $J run --samples $D/samples.tsv --ref $REF/hg38_parmask --out $D/work --threads $TH --shallow $a --clonal $b > $D/pipeline.log 2>&1
        if [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ]; then
            mkdir -p $D/logs; mv $D/work/calls.vcf $D/work/n.cnt $D/work/t.cnt $D/work/t.cnt.hist $D/; mv $D/work/log.* $D/work/steps.tsv $D/logs/
            discard $D/work; log "$($SW check $k $n $D)"
        else log "[$k/$n tau $a/$b] FAILED: $(tail -1 $D/pipeline.log)"; [ -e $D/attempt ] && touch $D/failed || touch $D/attempt; fi
    }
    room() { df --output=avail -BG /home | tail -1 | tr -dc 0-9; }
    for p in ${*:-$CASTLE}; do
        MX=$W/mix/$p
        while :; do
            L=$($SW plan $p | while read k n a b m0 mt; do [ -e $W/runs/$k/$n/t${a}_$b/calls.vcf ] || [ -e $W/runs/$k/$n/t${a}_$b/failed ] || echo "$k $n $a $b"; done)
            [ -z "$L" ] && break
            NEW=$(echo "$L" | while read k n a b; do [ -e $W/runs/$k/$n/t${a}_$b/work/done.extract_t ] || echo "$k $n $a $b"; done)
            if echo "$NEW" | grep -q '^coverage\|^purity'; then
                mkdir -p $MX
                for s in n t; do [ -s $MX/$s.svb.idx ] || $J pack --reads $(src $p $s) --out $MX/$s.svb --threads $TH --shuffle > $MX/pack_$s.log 2>&1; done
            fi
            if echo "$NEW" | grep -q '^purity'; then
                nb=$(python3 -c "import struct; print(struct.unpack_from('<Q', open('$MX/n.svb', 'rb').read(32), 16)[0])")
                nm=$((nb - 24 * G - G / 10)); [ $nm -gt $((30 * G)) ] && nm=$((30 * G))
                [ -s $MX/nm.svb.idx ] || $TOOLS mix $MX/nm.svb $MX/n.svb:0:$nm >> $Q 2>&1
                for q in $(echo "$NEW" | awk '$1 == "purity" {sub(/.*_/, "", $2); sub(/p$/, "", $2); print $2}' | sort -u); do
                    [ -s $MX/t$q.svb.idx ] || $TOOLS mix $MX/t$q.svb $MX/t.svb:0:$((q * 30 * G / 100)) $([ $q -lt 100 ] && echo $MX/n.svb:$nm:$(((100 - q) * 30 * G / 100))) >> $Q 2>&1
                done
            fi
            echo "$L" | { while read k n a b; do while [ $(jobs -rp | wc -l) -ge 2 ] || { [ $(jobs -rp | wc -l) -gt 0 ] && [ $(room) -lt 150 ]; }; do wait -n; done; one $k $n $a $b & done; wait; }
        done
        discard $MX
    done ;;
score)
    $SW configs > $W/configs.txt 2> $W/configs.err
    for k in $(cat $W/configs.txt); do
        C=$W/conf/$k; mkdir -p $C/verdict
        [ -s $C/verdict/folds.json ] || ( cd $H && OMP_NUM_THREADS=$TH $PY $H/classifier.py verdict --runs $C/full --coverage $C/coverage --purity $C/purity --out $C/verdict > $C/verdict/stdout.txt 2>&1 )
        log "[config $k] $(tail -1 $C/verdict/log.txt)"
        E=$C/eval_sets VD=$C/verdict bash $H/consensus.sh with
        E=$C/eval_sets LOO=$C/eval_loo bash $H/consensus.sh loo
        E=$C/eval_sets VD=$C/verdict O=$C/coverage_eval bash $H/coverage.sh score $CASTLE
        $SW rescore $C > $C/rescore.tsv; log "[config $k] done"
    done
    $SW rescores > $W/rescore.txt ;;
*) echo "unknown step $X" >&2; exit 1 ;;
esac
