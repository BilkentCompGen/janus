#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/cue2:latest
GPU=1
SHM=16g
for s in t n; do
  b=$T; [ $s = n ] && b=$N
  mkdir -p $O/$s.work
  printf 'bam: "%s"\nfa: "%s"\nmodel_path: "/opt/models/cue2.v2.pt"\nchr_names: [%s]\nn_cpus: %s\ngpu_ids: [0]\nn_jobs_per_gpu: 3\n' $b $REF "$(printf '"chr%s",' {1..22} X Y | sed 's/,$//')" $P > $O/$s.work/call.yaml
  dk $IMG cue call --config $O/$s.work/call.yaml
  cp $O/$s.work/results/cue_svs.vcf $O/$s.vcf
done
python3 $K/subtract.py $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf
done_msg
