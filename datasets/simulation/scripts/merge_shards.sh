#!/bin/bash
# Per sample: concatenate the gz shards into one <sample>.fastq.gz, verify the result is a
# byte-exact concatenation, then move the shards to trash. gzip members concatenate, so this
# is lossless and needs no recompression. Per-shard badread .log files stay where they are.
set -u
R=$(cd "$(dirname "$0")/.." && pwd)/reads
TRASH=$HOME/trash/$(date +%Y%m%d)-simulation-shards
mkdir -p "$TRASH"

s=$1
out=$R/$s/$s.fastq.gz
if [ -s "$out" ]; then echo "[$s] zaten var, atlandi"; exit 0; fi

mapfile -t sh < <(ls -1 "$R/$s/shards"/*.fastq.gz | sort)
n=${#sh[@]}
want=0
for f in "${sh[@]}"; do want=$((want + $(stat -c%s "$f"))); done

echo "[$s] $n shard, $want bayt -> $out"
t0=$(date +%s)
cat "${sh[@]}" > "$out.tmp" || { echo "[$s] CAT HATASI"; exit 1; }
got=$(stat -c%s "$out.tmp")
t1=$(date +%s)

if [ "$got" != "$want" ]; then
    echo "[$s] BOYUT UYUSMADI got=$got want=$want -- shard lar SILINMEDI"
    exit 1
fi
mv "$out.tmp" "$out"
mkdir -p "$TRASH/$s"
mv "${sh[@]}" "$TRASH/$s/"
echo "[$s] TAMAM $got bayt, $((t1-t0)) sn, $(( want / 1048576 / (t1-t0+1) )) MB/s, shard lar trash e alindi"
