#!/usr/bin/env bash
# compete.sh [pair ...]: the competitors on the same reads as snape, one pair after the other, on cores 0-15 with 16
# threads (as snape ran), measured like snape's runs. Per pair:
#   1. align normal and tumor (minimap2 -ax map-hifi -Y --secondary=no | samtools sort, one sample after the other),
#      check the primary reads against snape's pack counts      -> runs/align/<pair>/{n,t}/run.json
#   2. every tool on the full BAMs                               -> runs/competitors/<tool>/<pair>/
#   3. the core tools on the depth subsets 5x 10x 20x 40x (as deep as both samples go): the same reads snape's coverage
#      runs used (the first N x 3.1 Gb of each fastq.gz, by read name)   -> runs/competitors/<tool>/<pair>_<N>x/
#   4. the BAMs, the subsets and the tools' intermediates are deleted (approved by the user, 2026-10-02); each run keeps
#      calls.vcf, log.txt, run.json, files.tsv (sizes of everything before deletion).
# Docker tools run in a cgroup slice per run, so their CPU time and peak memory are measured (measure_tool.py).
# Resumes: a run with run.json is skipped. Progress: runs/competitors/queue.log.
S=/home/yuio/snape; H=$S/bench; K=$H/callers; R=$H/runs; G=3100000000; IDX=$S/ref/hg38_parmask.mmi
FULL="severus severus_phased savana nanomonsv sniffles2 sniffles2_mosaic svision-pro delly minisv somasv longcalld svscope"
CORE="severus savana nanomonsv sniffles2 svision-pro"
LEVELS="5 10 20 40"
ALL="hcc1954 h2009 hcc1937 h1437 hg008 colo829 hcc1395 hg009 h209 h2126"
PAIRS="${*:-$ALL}"
export CPUS=0-15 DMEM=120g
Q=$R/competitors/queue.log; mkdir -p $R/competitors $R/align
log() { echo "$(date '+%F %T') $*" >> $Q; }

# measured run: one <dir> <params json> <command ...>; the command writes into <dir>/out
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
    # keep the calls and the log, delete the tool's intermediates
    for f in calls.vcf log.txt stdout.txt; do [ -e $D/out/$f ] && mv $D/out/$f $D/; done
    rm -rf $D/out
}

# align <pair> <n|t> <bam>: minimap2 | samtools sort, read group per sample
align() {
    local p=$1 s=$2 B=$3
    bash -c "minimap2 -ax map-hifi -t 16 -Y --secondary=no -R '@RG\tID:${p}_$s\tSM:${p}_$s' $IDX $S/hifi/$p/$s.fastq.gz 2> $B.minimap2.log \
        | samtools sort -@ 4 -m 2G -T $B.st -o $B.tmp - && mv $B.tmp $B && samtools index -@ 8 $B"
}

# names <fastq.gz> <bases> <out>: names of the first reads reaching <bases> bases (pack -b keeps the same reads)
names() { pigz -dc "$1" | awk -v lim=$2 'NR%4==1{n=substr($1,2)} NR%4==2{print n; s+=length($0); if(s>=lim) exit}' > "$3"; }

depth() { python3 -c "import json; r=json.load(open('$R/full/$1/run.json'))['reads']; print(int(min(r['n']['bases'], r['t']['bases']) // $G))"; }

for p in $PAIRS; do
    W=$R/competitors/work/$p; mkdir -p $W
    # 1. alignment, measured per sample; the BAM is checked against snape's read counts
    for s in n t; do
        A=$R/align/$p/$s
        [ -s $A/run.json ] || one $A "{\"experiment\": \"align\", \"tool\": \"minimap2+samtools\", \"pair\": \"$p\", \"sample\": \"$s\", \"bam\": \"$W/$s.bam\"}" \
            bash -c "$(declare -f align); S=$S IDX=$IDX align $p $s $W/$s.bam"
        if [ ! -s $W/$s.bam.bai ]; then log "[$p] alignment of $s missing or failed; pair skipped"; continue 2; fi
        want=$(python3 -c "import json; print(json.load(open('$R/full/$p/run.json'))['reads']['$s']['reads'])")
        got=$(taskset -c 0-15 samtools view -c -@ 8 -F 0x900 $W/$s.bam)
        [ "$want" = "$got" ] && log "[$p] $s: $got reads aligned" || { log "[$p] $s: $got reads in the BAM, $want in snape's pack: pair skipped"; continue 2; }
    done
    # 2. every tool on all reads
    for t in $FULL; do
        one $R/competitors/$t/$p "{\"experiment\": \"competitor\", \"tool\": \"$t\", \"pair\": \"$p\", \"coverage_x\": null}" \
            $K/$t.sh $W/t.bam $W/n.bam $R/competitors/$t/$p/out 16
    done
    # 3. the core tools at the depths snape's coverage runs used, on the same reads
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
    # 4. the pair's BAMs go
    for s in n t; do python3 - "$R/align/$p/$s/run.json" "$W/$s.bam" <<'EOF'
import json, os, sys
j = json.load(open(sys.argv[1])); j["bam_bytes"] = os.path.getsize(sys.argv[2]); json.dump(j, open(sys.argv[1], "w"), indent=1)
EOF
    done
    rm -rf $W
    log "[$p] done, BAMs and intermediates deleted"
done
log "queue done"
