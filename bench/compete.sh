#!/usr/bin/env bash
S=/home/yuio/snape; H=$S/bench; K=$H/callers; R=$H/runs; G=3100000000; IDX=$S/ref/hg38_parmask.mmi
FULL="severus severus_phased savana nanomonsv sniffles2 sniffles2_mosaic svision-pro delly minisv somasv longcalld svscope"
CORE="severus savana nanomonsv sniffles2 svision-pro"
LEVELS=""
ALL="hcc1954 h2009 hcc1937 h1437 hcc1395"
PAIRS="${*:-$ALL}"
export CPUS=0-15 DMEM=120g
Q=$R/competitors/queue.log; mkdir -p $R/competitors $R/align
log() { echo "$(date '+%F %T') $*" >> $Q; }

one() {
    local D=$1 P=$2; shift 2
    [ -s $D/run.json ] && return 0
    mkdir -p $D/out; echo "$P" > $D/params.json
    local sl="sb_$(basename $(dirname $D) | tr -c 'a-z0-9\n' '_')_$(basename $D | tr -c 'a-z0-9\n' '_').slice"
    log "[$(basename $(dirname $D))/$(basename $D)] start"
    CGPARENT=$sl taskset -c 0-15 /usr/bin/time -f "%e %M %U %S %x" -o $D/time.txt "$@" > $D/out/stdout.txt 2>&1
    local c=/sys/fs/cgroup/$sl
    [ -d $c ] && { awk '{print}' $c/cpu.stat | grep -E "^(usage|user|system)_usec"; echo "memory_peak $(cat $c/memory.peak)"; } > $D/cgroup.txt
    [ -d $c ] && rmdir $c 2>/dev/null
    python3 $H/measure_tool.py $D >> $Q 2>&1
    for f in calls.vcf log.txt stdout.txt; do [ -e $D/out/$f ] && mv $D/out/$f $D/; done
    rm -rf $D/out
}

align() {
    local p=$1 s=$2 B=$3
    bash -c "minimap2 -ax map-hifi -t 16 -Y --secondary=no -R '@RG\tID:${p}_$s\tSM:${p}_$s' $IDX $S/hifi/$p/$s.fastq.gz 2> $B.minimap2.log \
        | samtools sort -@ 4 -m 2G -T $B.st -o $B.tmp - && mv $B.tmp $B && samtools index -@ 8 $B"
}

names() { pigz -dc "$1" | awk -v lim=$2 'NR%4==1{n=substr($1,2)} NR%4==2{print n; s+=length($0); if(s>=lim) exit}' > "$3"; }

depth() { python3 -c "import json; r=json.load(open('$R/full/$1/run.json'))['reads']; print(int(min(r['n']['bases'], r['t']['bases']) // $G))"; }

for p in $PAIRS; do
    W=$R/competitors/work/$p; mkdir -p $W
    for s in n t; do
        A=$R/align/$p/$s
        [ -s $A/run.json ] || one $A "{\"experiment\": \"align\", \"tool\": \"minimap2+samtools\", \"pair\": \"$p\", \"sample\": \"$s\", \"bam\": \"$W/$s.bam\"}" \
            bash -c "$(declare -f align); S=$S IDX=$IDX align $p $s $W/$s.bam"
        if [ ! -s $W/$s.bam.bai ]; then log "[$p] alignment of $s missing or failed; pair skipped"; continue 2; fi
        want=$(python3 -c "import json; print(json.load(open('$R/full/$p/run.json'))['reads']['$s']['reads'])")
        got=$(taskset -c 0-15 samtools view -c -@ 8 -F 0x900 $W/$s.bam)
        [ "$want" = "$got" ] && log "[$p] $s: $got reads aligned" || { log "[$p] $s: $got reads in the BAM, $want in snape's pack: pair skipped"; continue 2; }
    done
    for t in $FULL; do
        one $R/competitors/$t/$p "{\"experiment\": \"competitor\", \"tool\": \"$t\", \"pair\": \"$p\", \"coverage_x\": null}" \
            $K/$t.sh $W/t.bam $W/n.bam $R/competitors/$t/$p/out 16
    done
    d=$(depth $p)
    for n in $LEVELS; do
        [ $n -le $d ] || continue
        for s in n t; do
            [ -s $W/${s}_${n}x.bam.bai ] && continue
            names $S/hifi/$p/$s.fastq.gz $((n * G)) $W/${s}_${n}x.names
            taskset -c 0-15 samtools view -@ 8 -b -N $W/${s}_${n}x.names -o $W/${s}_${n}x.bam $W/$s.bam && taskset -c 0-15 samtools index -@ 8 $W/${s}_${n}x.bam
        done
        for t in $CORE; do
            one $R/competitors/$t/${p}_${n}x "{\"experiment\": \"competitor\", \"tool\": \"$t\", \"pair\": \"$p\", \"coverage_x\": $n, \"bases_limit\": $((n * G))}" \
                $K/$t.sh $W/t_${n}x.bam $W/n_${n}x.bam $R/competitors/$t/${p}_${n}x/out 16
        done
        rm -f $W/n_${n}x.* $W/t_${n}x.*
    done
    for s in n t; do python3 - "$R/align/$p/$s/run.json" "$W/$s.bam" <<'EOF'
import json, os, sys
j = json.load(open(sys.argv[1])); j["bam_bytes"] = os.path.getsize(sys.argv[2]); json.dump(j, open(sys.argv[1], "w"), indent=1)
EOF
    done
    log "[$p] done, BAMs kept in $W"
done
log "queue done"
