#!/bin/bash
# Fig. 3 and the simulation tables: the 20-tumor simulation (datasets/simulation) at five normal/tumor depths, as pairwise
# runs and as joint runs, classified and scored against the simulated SVs.
#   bench/simulation.sh truth     the simulated SVs of every tumor as reference junctions, and those on a lost chromosome copy (runs/sim/truth)
#   bench/simulation.sh reads     the packed read sets at 20x, 10x and 5x (runs/sim/reads/x<depth>)
#   bench/simulation.sh pairs     every tumor with the normal at the five depth settings, three at a time (runs/sim/settings/<n>_<t>/<tumor>)
#   bench/simulation.sh joint     the joint run of the 20 tumors at every setting (runs/sim/settings/<n>_<t>/joint)
#   bench/simulation.sh pooled    the 20 tumors at 5x pooled into one tumor sample with the normal at 20x (runs/sim/pooled/20_5)
#   bench/simulation.sh score     classification of the pairwise runs, the joint runs with every acceptance rule of the supplement, the scores (runs/sim/scores.json)
source "$(dirname "$0")/env.sh"
X=${1:?see the header}; shift
PH=$DS/simulation; E=$SIM/settings; M=$S/models/janus.model; RF=$REF/chm13v2
SETTINGS="5 5|10 10|20 10|20 5|20 20"
case $X in
truth)
    $TOOLS simtruth $PH/genomes $SIM/truth && $TOOLS lost $PH/genomes $SIM/truth $RF.fa.fai ;;
reads)
    # 20x: the tumors as packed by the simulation; the normal's reads trimmed of the adapters Badread adds (svb2fq) and packed again.
    # 10x and 5x: every second and every fourth read.
    [ -x $H/svb2fq ] || gcc -O2 -o $H/svb2fq $H/svb2fq.c
    mkdir -p $SIM/reads/x20
    [ -s $SIM/reads/x20/germline.svb.idx ] || { $H/svb2fq $PH/reads/germline/germline.svb germline 2> $SIM/reads/x20/germline.log | pigz -p $TH -3 > $SIM/reads/x20/germline.fastq.gz && $J pack --reads $SIM/reads/x20/germline.fastq.gz --out $SIM/reads/x20/germline.svb --threads $TH && discard $SIM/reads/x20/germline.fastq.gz; }
    for t in $TUMORS; do for e in svb svb.idx; do ln -sf $PH/reads/$t/$t.$e $SIM/reads/x20/$t.$e; done; done
    for x in 5 10; do for s in germline $TUMORS; do [ -s $SIM/reads/x$x/$s.svb.idx ] || { mkdir -p $SIM/reads/x$x; $TOOLS thin $PH/reads/$s/$s.svb $SIM/reads/x$x/$s.svb $((20 / x)) >> $SIM/reads/x$x.log; }; done; done ;;
pairs)
    one() {
        local nx=$1 tx=$2 t=$3 N=$E/${1}_${2} p=$E/${1}_${2}/$3
        [ -s $p/t.cnt.hist ] && return; mkdir -p $p/work
        sheet $p/samples.tsv $(simreads $nx germline) $(simreads $tx $t)
        if [ -s $N/n.svb.exc ] && [ -s $N/n.cnt ]; then   # the normal's cores and counts are the same for every tumor of a setting
            ln -sf $(simreads $nx germline) $p/work/n.svb; ln -sf $(simreads $nx germline).idx $p/work/n.svb.idx
            ln -sf $N/n.svb.exc $p/work/n.svb.exc; ln -sf $N/n.cnt $p/work/n.cnt; touch $p/work/done.extract_n $p/work/done.count_n
        fi
        /usr/bin/time -f "%e %M %U %S %x" -o $p/time.txt $J run --samples $p/samples.tsv --ref $RF --out $p/work --threads $TH > $p/pipeline.log 2>&1
        cp $p/work/steps.tsv $p/work/calls.vcf $p/work/t.cnt.hist $p/work/t.cnt $p/work/t.svb.exc $p/ 2> /dev/null
        for f in n.cnt n.svb.exc; do [ -s $N/$f ] || cp $p/work/$f $N/; done
        [ -s $p/calls.vcf ] && [ -s $p/t.svb.exc ] && discard $p/work
    }
    export -f one sheet simreads discard; export E J RF TH SIM TRASH
    echo "$SETTINGS" | tr '|' '\n' | while read nx tx; do
        N=$E/${nx}_${tx}; mkdir -p $N
        [ -s $N/n.cnt ] || for t in $TUMORS; do [ -s $N/$t/t.cnt.hist ] && continue; one $nx $tx $t; break; done
        echo $TUMORS | tr ' ' '\n' | xargs -P 3 -I{} bash -c "one $nx $tx {}"
        log "pairs $nx/$tx: $(ls $N/*/t.cnt.hist 2> /dev/null | wc -l) of 20"
    done ;;
joint)
    echo "$SETTINGS" | tr '|' '\n' | while read nx tx; do
        d=$E/${nx}_${tx}/joint; [ -s $d/calls.vcf ] && continue; mkdir -p $d/work
        jointsheet $d/samples.tsv $nx $tx
        $J run --samples $d/samples.tsv --ref $RF --out $d/work --threads $TH > $d/pipeline.log 2>&1
        cp $d/work/steps.tsv $d/work/calls.vcf $d/work/*.cnt $d/work/*.cnt.hist $d/work/*.svb.exc $d/
        for t in $TUMORS; do python3 -c "print(open('$d/$t.cnt', 'rb').read(27)[26])" > $d/$t.minimum; done
        [ -s $d/calls.vcf ] && discard $d/work; log "joint $nx/$tx: $(grep -vc '^#' $d/calls.vcf 2> /dev/null) records"
    done ;;
pooled)
    set -- ${1:-20} ${2:-5}; P=$SIM/pooled/${1}_${2}; mkdir -p $P/work
    [ -s $P/work/t.svb.idx ] || $TOOLS svbcat $P/work/t.svb $(for t in $TUMORS; do simreads $2 $t; done) > $P/svbcat.log
    sheet $P/samples.tsv $(simreads $1 germline) $P/work/t.svb
    /usr/bin/time -f "%e %M %U %S %x" -o $P/time.txt $J run --samples $P/samples.tsv --ref $RF --model $M --out $P/work --threads $TH > $P/pipeline.log 2>&1
    cp $P/work/steps.tsv $P/work/calls.vcf $P/work/janus.vcf $P/work/t.cnt.hist $P/; discard $P/work ;;
score)
    echo "$SETTINGS" | tr '|' '\n' | while read nx tx; do
        d=$E/${nx}_${tx}
        for t in $TUMORS; do [ -s $d/$t/janus.vcf ] || $J classify --calls $d/$t/calls.vcf --hist $d/$t/t.cnt.hist --model $M --out $d/$t/janus.vcf; done
        HS=$(for t in $TUMORS; do printf "%s/joint/%s.cnt.hist," $d $t; done | sed 's/,$//'); MS=$(for t in $TUMORS; do cat $d/joint/$t.minimum; done | paste -sd,)
        jt() { [ -s $1 ] || $J joint --calls $d/joint/calls.vcf --hist $HS --minimum $MS --model $M --out $1 "${@:2}"; }
        jt $d/joint/janus.vcf; jt $d/joint/first.vcf --no-clade
        for mu in 0.3 0.5 1.0 1.5 2.0; do jt $d/joint/mu_$mu.vcf --mu $mu; done
        if [ ${nx}_${tx} = 20_20 ] || [ ${nx}_${tx} = 20_5 ]; then for r in 0.0 0.1 0.2 0.3; do for f in 0.05 0.1 0.2 0.3; do jt $d/joint/rule_${r}_$f.vcf --no-clade --raise $r --floor $f; done; done; fi
    done
    $TOOLS simjoint | tee $SIM/scores.txt
    [ -s $SIM/pooled/20_5/janus.vcf ] && $TOOLS pooledscore 20_5 ;;
*) echo "unknown step $X" >&2; exit 1 ;;
esac
