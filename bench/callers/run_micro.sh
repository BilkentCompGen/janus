#!/usr/bin/env bash
K=/home/yuio/snape/bench/callers; M=$K/micro; [ "${PAIR:-}" = demo ] && M=$K/micro/demo
for t in "$@"; do
  o=$M/out/$t
  if [ -e $o ]; then mkdir -p /home/yuio/trash/20261002-competitors/micro-out; mv $o /home/yuio/trash/20261002-competitors/micro-out/${PAIR:-micro}.$t.$(date +%H%M%S); fi
  mkdir -p $o
  ( while :; do docker stats --no-stream --format '{{.Name}} {{.MemUsage}}' 2>/dev/null; sleep 2; done ) > $o/mem.txt & sp=$!
  taskset -c 16-31 /usr/bin/time -f "%e %M %x" -o $o/time.txt $K/$t.sh $M/t.bam $M/n.bam $o 16
  kill $sp 2>/dev/null
  read wall kb ex < <(tail -1 $o/time.txt)
  dmem=$(awk '{v=$2; u=v; gsub(/[0-9.]/,"",u); gsub(/[^0-9.]/,"",v); m=(u=="GiB")?v*1048576:(u=="MiB")?v*1024:(u=="KiB")?v:0; if(m>x)x=m}END{printf "%d", x}' $o/mem.txt)
  [ "${dmem:-0}" -gt "$kb" ] && kb=$dmem
  n=$( [ -s $o/calls.vcf ] && (grep -vc "^#" $o/calls.vcf; true) || echo NA)
  printf "%s\t%s\t%s\t%s\t%s\n" $t $ex $wall $kb $n | tee -a $M/results.tsv
done
