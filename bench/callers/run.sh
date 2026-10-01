#!/bin/bash
# run.sh: the other callers from the BAMs (PAR-masked hg38), timed: Severus on all five pairs, Sniffles2 (mosaic mode, tumor only)
# on all five, SAVANA and nanomonsv on hg008; hg008 results scored against NIST. Outputs under bench/callers/.
S=/home/yuio/snape; R=$S/ref; H=$S/hifi; L=$S/lib/callers; K=$S/bench/callers; T=16
tm() { d=$1; shift; mkdir -p $d; /usr/bin/time -f "$d %e s wall %U s user %M KB exit %x" "$@" > $d/log 2> $d/time; tail -1 $d/time >> $K/times.txt; }
for p in hg008 hcc1954 hcc1937 h1437 h2009; do
  [ -e $K/severus/$p/done ] && continue
  tm $K/severus/$p $L/severus/bin/severus --target-bam $H/$p/t.bam --control-bam $H/$p/n.bam --out-dir $K/severus/$p/out --vntr-bed $R/vntr.trf.bed -t $T && touch $K/severus/$p/done
  [ $p = hg008 ] && [ -s $K/severus/hg008/out/somatic_SVs/severus_somatic.vcf ] && $S/bench/eval.sh hg008 $K/severus/hg008/out/somatic_SVs/severus_somatic.vcf severus > $K/severus/hg008/eval.txt 2>&1
done
for p in hg008 hcc1954 hcc1937 h1437 h2009; do
  [ -e $K/sniffles/$p/done ] && continue
  tm $K/sniffles/$p bash $K/sniffles_somatic.sh $p && touch $K/sniffles/$p/done
done
[ -s $K/sniffles/hg008/sniffles.vcf ] && $S/bench/eval.sh hg008 $K/sniffles/hg008/sniffles.vcf sniffles > $K/sniffles/hg008/eval.txt 2>&1
tm $K/savana/hg008 $L/savana/bin/savana -t $H/hg008/t.bam -n $H/hg008/n.bam --ref $R/hg38_parmask.fa --outdir $K/savana/hg008/out --threads $T --pb --sample hg008
v=$(ls $K/savana/hg008/out/*somatic*.vcf 2>/dev/null | head -1); [ -n "$v" ] && $S/bench/eval.sh hg008 $v savana > $K/savana/hg008/eval.txt 2>&1
d=$K/nanomonsv/hg008; mkdir -p $d
tm $d/parse_t $L/nanomonsv/bin/nanomonsv parse $H/hg008/t.bam $d/t --reference_fasta $R/hg38_parmask.fa
tm $d/parse_n $L/nanomonsv/bin/nanomonsv parse $H/hg008/n.bam $d/n --reference_fasta $R/hg38_parmask.fa
tm $d/get $L/nanomonsv/bin/nanomonsv get $d/t $H/hg008/t.bam $R/hg38_parmask.fa --control_prefix $d/n --control_bam $H/hg008/n.bam --processes $T
[ -s $d/t.nanomonsv.result.vcf ] && $S/bench/eval.sh hg008 $d/t.nanomonsv.result.vcf nanomonsv > $d/eval.txt 2>&1
cat $K/*/hg008/eval.txt 2>/dev/null | grep A_nist > $K/scores_hg008.txt
touch $K/done
