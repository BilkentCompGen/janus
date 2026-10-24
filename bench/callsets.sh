#!/bin/bash
# The call sets of every experiment of the paper, gzipped, in bench/callsets/<experiment>/...; janus's files are its complete
# output (every record with FILTER PASS, LOW or DUPL and the probability VP), the other callers' files their PASS calls.
#   full/<pair>/        the seven pairs: on the five CASTLE pairs the 18 call sets of the consensus (janus and ten HiFi callers,
#                       the published CASTLE ONT and Illumina sets; Minda was given janus's PASS records) plus Sniffles2 mosaic
#                       and Delly (HCC1954 only); on HG008 and COLO829
#                       janus and the HiFi callers as scored against their truth sets
#   depth/<pair>/<n>x/  janus, Severus, SAVANA and nanomonsv at 30, 20, 15, 10 and 5x (the other callers on the CASTLE pairs only)
#   purity/<pair>/<q>p/ the same on tumors of 100 to 20% purity at 30x
#   simulation/<n>_<t>/ the 20 pairwise runs (pairwise/T01..T20) and the joint run of the simulated cohort at each depth setting,
#                       and pooled_20_5: the 20 tumors at 5x pooled into one sample
#   bench/callsets.sh   (after consensus.sh, coverage.sh and simulation.sh)
source "$(dirname "$0")/env.sh"
O=$H/callsets; E=$H/eval_sets; VD=$H/verdict
put() { [ -s $1 ] && gzip -9 -c $1 > $2.vcf.gz; }
name() { case $1 in severus*|savana*|nanomonsv|sniffles2) echo ${1%_phased}_hifi ;; *) echo $1 ;; esac; }
for p in $CASTLE; do
    castle $p; d=$O/full/$p; mkdir -p $d
    put $VD/$p.vcf $d/janus
    for t in $HIFI_SETS sniffles2_mosaic delly; do put $E/$p/in/$t.vcf $d/$(name $t); done
    set -- $ONT; put $1 $d/severus_ont; put $2 $d/savana_ont; put $3 $d/nanomonsv_ont; put $4 $d/sniffles2_ont
    set -- $SR; put $1 $d/svaba; put $2 $d/gripss; put $3 $d/manta
done
for p in hg008 colo829; do
    d=$O/full/$p; mkdir -p $d
    put $VD/$p.vcf $d/janus
    for t in $HIFI_SETS sniffles2_mosaic; do put $R/external_eval/$p/$t/in.vcf $d/$(name $t); done
done
for p in $LABELED; do
    for s in 30x 20x 15x 10x 5x 100p 80p 60p 40p 20p; do
        k=depth; v=$VD/coverage_matched/${p}_$s.vcf; [ ${s: -1} = p ] && { k=purity; v=$VD/purity_matched/${p}_$s.vcf; }
        [ -s $v ] || continue; d=$O/$k/$p/$s; mkdir -p $d
        put $v $d/janus
        for t in severus savana nanomonsv; do put $R/coverage_callers/$t/${p}_$s/calls.vcf $d/${t}_hifi; done
    done
done
for s in 20_20 20_10 20_5 10_10 5_5; do
    d=$O/simulation/$s; mkdir -p $d/pairwise
    for t in $TUMORS; do put $SIM/settings/$s/$t/janus.vcf $d/pairwise/$t; done
    put $SIM/settings/$s/joint/janus.vcf $d/joint
done
put $SIM/pooled/20_5/janus.vcf $O/simulation/pooled_20_5
echo "$(find $O -name '*.vcf.gz' | wc -l) call sets, $(du -sh $O | cut -f1)"
