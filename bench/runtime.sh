#!/bin/bash
# Fig. 4 and the run time tables: the CASTLE pairs and the simulated runs once more, one at a time on an otherwise idle
# machine, measured (runs_clean/full/<pair>, runs_clean/sim/timing/<depths>/{pair/<tumor>,joint}).
#   bench/runtime.sh
source "$(dirname "$0")/env.sh"
Q=$RC/queue.log; mkdir -p $RC
for p in $CASTLE; do mkdir -p $RC/full; sheet $RC/full/$p.tsv $(src $p n) $(src $p t); janus_run full $RC/full/$p $p $RC/full/$p.tsv ""; done
source $H/simulation.sh reads   # the packed read sets of the simulation at every depth
tm() { sync; sleep 10; /usr/bin/time -f "%e %M %U %S %x" -o $1 taskset -c $CPUS "${@:2}"; }
for set in 20_20 20_10 20_5 10_10 5_5; do
    nx=${set%_*}; tx=${set#*_}
    for t in $(seq -f "T%02g" 1 20); do
        p=$RC/sim/timing/$set/pair/$t; [ -s $p/time.txt ] && continue; mkdir -p $p
        sheet $p/samples.tsv $(simreads $nx germline) $(simreads $tx $t)
        REFNAME=chm13v2 tm $p/time.part $J run --samples $p/samples.tsv --ref $REF/chm13v2 --out $p/work --threads $TH > $p/pipeline.log 2>&1
        cp $p/work/steps.tsv $p/work/calls.vcf $p/; [ -s $p/calls.vcf ] && mv $p/time.part $p/time.txt && discard $p/work
    done
    d=$RC/sim/timing/$set/joint; [ -s $d/time.txt ] && continue; mkdir -p $d
    jointsheet $d/samples.tsv $nx $tx
    tm $d/time.part $J run --samples $d/samples.tsv --ref $REF/chm13v2 --out $d/work --threads $TH > $d/pipeline.log 2>&1
    cp $d/work/steps.tsv $d/work/calls.vcf $d/; [ -s $d/calls.vcf ] && mv $d/time.part $d/time.txt && discard $d/work
done
# minimap2 on the 21 simulated read sets at 5/5 and 20/20, one at a time (runs/sim/mm2/<depths>/<sample>.time, Fig. 3b);
# figures.py uses the stored 5/5 measurement until these files exist
MMI=$K/ref/chm13v2.map-hifi.mmi; [ -s $MMI ] || minimap2 -x map-hifi -d $MMI $REF/chm13v2.fa
[ -x $H/svb2fq ] || gcc -O2 -o $H/svb2fq $H/svb2fq.c
mapped() { $H/svb2fq $1 $2 2> /dev/null | minimap2 -ax map-hifi -t $TH -Y --secondary=no $MMI - 2> $3.minimap2.log | samtools sort -@ 4 -m 2G -T $3.st -o $3.bam - && samtools index -@ 8 $3.bam; }
export -f mapped; export H TH MMI
for x in 5 20; do
    d=$SIM/mm2/${x}_${x}; mkdir -p $d
    for s in germline $TUMORS; do
        [ -s $d/$s.time ] && continue; sync; sleep 5
        /usr/bin/time -f "%e %M %U %S %x" -o $d/$s.time.part taskset -c $CPUS bash -lc "mapped $(readlink -f $(simreads $x $s)) $s $d/$s"
        discard $d/$s.bam; discard $d/$s.bam.bai; mv $d/$s.time.part $d/$s.time
    done
done
