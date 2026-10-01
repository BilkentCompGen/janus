#!/bin/bash
# repro.sh: do our own runs of Severus and Sniffles2 (from our BAMs) reproduce the published CASTLE call sets? For every
# CASTLE pair whose run is finished: our calls and the published calls scored against the same consensus, and our calls
# scored against the published ones (recall = share of the published calls we also have, precision = share of ours that
# are published). Writes repro.txt; safe to re-run while the queue is still going.
S=/home/yuio/snape; K=$S/bench/callers; B=$S/hifi/castle-calls/variant_calls_and_benchmarks; T=$S/hifi/truth
sc() { python3 $S/bench/score.py $1 $2; }
{
for p in hcc1954 h2009 hcc1937 h1437; do
  P=${p^^}; D=$B/$P/$P
  v=$K/severus/$p/out/somatic_SVs/severus_somatic.vcf
  if [ -e $K/severus/$p/done ] && [ -s $v ]; then
    echo "$p severus own_vs_consensus $(sc $T/castle_without_severus/$p/None_minda_ensemble.vcf $v)"
    echo "$p severus published_vs_consensus $(sc $T/castle_without_severus/$p/None_minda_ensemble.vcf $D/severus_HIFI_vaf.vcf)"
    echo "$p severus own_vs_published $(sc $D/severus_HIFI_vaf.vcf $v)"
  fi
  v=$K/sniffles/$p/sniffles.vcf
  if [ -e $K/sniffles/$p/done ] && [ -s $v ]; then
    echo "$p sniffles2 own_vs_consensus $(sc $T/castle_without_snape/$p/None_minda_ensemble.vcf $v)"
    echo "$p sniffles2 published_vs_consensus $(sc $T/castle_without_snape/$p/None_minda_ensemble.vcf $D/sniffles2_HIFI_vaf.vcf)"
    echo "$p sniffles2 own_vs_published $(sc $D/sniffles2_HIFI_vaf.vcf $v)"
  fi
done
} | tee $K/repro.txt
