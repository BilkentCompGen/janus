#!/bin/bash
# The classifier: the leave-one-pair-out classification of every janus run (the janus call sets of the paper, bench/verdict),
# the shipped models, and the supplementary experiments on the classifier.
#   bench/classifier.sh verdict    train on the 17-call-set consensus labels, classify the full, depth and purity runs with each CASTLE pair held out; models into models/
#   bench/classifier.sh labels     the same with the published 11-call-set consensus as labels (supp/labels11), scored
#   bench/classifier.sh supp       held-out AUC, all 3-train/2-test splits, the parameter ablation, the first-pass duplicates and the feature groups (supp)
#   bench/classifier.sh grid       the 3^5 parameter grid, four shards (supp/grid.tsv)
source "$(dirname "$0")/env.sh"
X=${1:?verdict, labels, supp or grid}; shift
C="$PY $H/classifier.py"
case $X in
verdict)
    for p in $CASTLE; do [ -s $H/eval_sets/$p/all_without/None_minda_ensemble.vcf ] || bash $H/consensus.sh without $p; done
    mkdir -p $H/verdict $S/models; ( cd $H && OMP_NUM_THREADS=$TH taskset -c $CPUS $C verdict --models $S/models > $H/verdict/stdout.txt 2>&1 ); cat $H/verdict/log.txt ;;
labels)
    V11=$H/supp/labels11; O=$V11/minda; E=$H/eval_sets; mkdir -p $O
    ( cd $H && OMP_NUM_THREADS=$TH taskset -c $CPUS $C verdict --out $V11 --labels $HIFI/truth/castle11 > $V11/stdout.txt 2>&1 )
    row() { awk -F'\t' -v k="$2" 'NR>1 && (k == "" || $1 == k) {printf "%.4f\t%.4f\t%.4f\n", $5, $6, $7; exit}' $1; }
    pair() {
        local p=$1 W=$O/$1 hifi="" t x; castle $p; mkdir -p $W
        grep -E '^#|	PASS	' $V11/$p.vcf > $W/janus.vcf
        for t in $HIFI_SETS; do hifi="$hifi $E/$p/in/$t.vcf"; done
        ens $W/with 11 $hifi $W/janus.vcf $ONT $SR
        tset $E/$p/all_without/None_minda_ensemble.vcf $W/janus.vcf $W/without
        for x in 30 20 15 10 5; do grep -E '^#|	PASS	' $V11/coverage_matched/${p}_${x}x.vcf > $W/d$x.vcf; tset $E/$p/all_with/None_minda_ensemble.vcf $W/d$x.vcf $W/d$x; done
        echo -e "$p\twith\t$(awk -F'\t' '$1 ~ /^K_/ {printf "%.4f\t%.4f\t%.4f", $5, $6, $7}' $W/with/results/overall.tsv)"
        echo -e "$p\twithout\t$(row $W/without/results/overall.tsv)"
        for x in 30 20 15 10 5; do echo -e "$p\t${x}x\t$(row $W/d$x/results/overall.tsv)"; done
    }
    ext() {
        local W=$O/$1; mkdir -p $W; grep -E '^#|	PASS	' $V11/$1.vcf > $W/janus.vcf; tset $2 $W/janus.vcf $W/truth
        echo -e "$1\ttruth\t$(row $W/truth/results/overall.tsv)\t$(grep -vc '^#' $W/janus.vcf)"
    }
    for p in $CASTLE; do pair $p > $O/$p.tsv & done
    ext hg008 $HIFI/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz > $O/hg008.tsv &
    ext colo829 $HIFI/truth/colo829/truth.vcf > $O/colo829.tsv &
    wait; cat $O/*.tsv > $V11/scores.tsv; cat $V11/scores.tsv ;;
supp)
    for j in auc splits ablation firstpass features; do ( cd $H && $C supp $j ) & done; wait ;;
grid)
    for i in 0 1 2 3; do ( cd $H && OMP_NUM_THREADS=4 $C supp grid --shard $i/4 > $H/supp/grid.$i.log 2>&1 ) & done; wait
    awk 'FNR == 1 && NR != 1 {next} !s[$0]++' $H/supp/grid.[0-9].tsv > $H/supp/grid.tsv ;;
*) echo "unknown step $X" >&2; exit 1 ;;
esac
