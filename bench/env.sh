# Shared settings and helpers of the benchmark scripts; every script sources this file.
S=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd); H=$S/bench; R=$H/runs; RC=$H/runs_clean; J=$S/build/janus; DS=$S/datasets; HIFI=$DS/hifi; REF=$S/ref
K=$H/callers; L=$K/env; PY=$S/.venv-bench/bin/python; TOOLS="$PY $H/tools.py"; MINDA=$HOME/.local/bin/minda
CASTLE="hcc1954 h2009 hcc1937 h1437 hcc1395"; LABELED="$CASTLE hg008 colo829"
TOOLS_ALL="severus_hp savana nanomonsv sniffles2 sniffles2_mosaic svision-pro cutesv debreak dysgu kled svdss delly"
HIFI_SETS="severus_phased savana_phased nanomonsv sniffles2 dysgu svision-pro cutesv debreak kled svdss"
G=3100000000; TH=16; CPUS=0-15
Q=${Q:-$R/queue.log}; mkdir -p $R
TRASH=$HOME/trash/$(date +%Y%m%d)-janus-bench   # intermediates of finished runs are moved here, never deleted
log() { echo "$(date '+%F %T') [$(basename $0 .sh)] $*" | tee -a $Q; }
discard() { mkdir -p $TRASH && mv "$1" $TRASH/$(basename $(dirname $1)).$(basename $1).$$; }

# the reads of a pair's sample: the BAM aligned by competitors.sh when it exists, else the fetched FASTQ
src() { [ -s $HIFI/$1/$2.bam ] && echo $HIFI/$1/$2.bam || echo $HIFI/$1/$2.fastq.gz; }

# a two-sample sheet for `janus run` (samples named n and t, as the result files are)
sheet() { printf "name\trole\tpath\nn\tnormal\t%s\nt\ttumor\t%s\n" $2 $3 > $1; }

# the published CASTLE call sets of a pair (HiFi, ONT and short-read) as given to Minda
castle() {
    local P=${1^^} SA=${1^^} GR=${1^^}_subsample.bam.gripss.filtered_vaf.vcf; CD=$HIFI/castle-calls/variant_calls_and_benchmarks/$P/$P
    [ $1 = hcc1395 ] && { SA=H1395; GR=H1395_realigned.bam.gripss.filtered_vaf.vcf; }
    PH="$CD/severus_HIFI_vaf.vcf $CD/savana_$SA.haplotagged.classified.somatic_HIFI_vaf.vcf $CD/tumor.nanomonsv.result_HIFI_vaf.vcf $CD/sniffles2_HIFI_vaf.vcf"
    ONT="$CD/severus_ont_vaf.vcf $CD/savana_$SA.haplotagged.classified.somatic_ONT_vaf.vcf $CD/tumor.nanomonsv.result_ONT_vaf.vcf $CD/sniffles2_ONT_vaf.vcf"
    SR="$CD/svaba_merged.vcf $CD/$GR $CD/manta_somaticSV_vaf.vcf"
}
norm() { case $2 in severus*|savana*|nanomonsv|sniffles2) cp $1 $3 ;; *) $TOOLS norm $1 $2 > $3 ;; esac; }
# Minda consensus of call sets: at least two of the three technologies and at least six call sets agree
ens() { local o=$1 k=$2; shift 2; $MINDA ensemble --vcfs "$@" --conditions "[[[caller_names[:$k], caller_names[$k:$((k + 4))], caller_names[$((k + 4)):]], '>=', 2],'&', [caller_names,'>=', 6]]" --min_size 50 --tolerance 500 --out_dir $o > $o.log 2>&1; }
tset() { $MINDA truthset --base $1 --vcfs $2 --min_size 50 --tolerance 500 --out_dir $3 > $3.log 2>&1; }
depth() { python3 -c "import json; r=json.load(open('$R/full/$1/run.json'))['reads']; print(int(min(r['n']['bases'], r['t']['bases']) // $G))"; }

# the simulation: the packed reads of a sample at a depth (simulation.sh reads) and the sheet of a joint run (normal at $2x, the 20 tumors at $3x)
SIM=$R/sim; TUMORS=$(seq -f "T%02g" 1 20 | tr '\n' ' ')
simreads() { echo $SIM/reads/x$1/$2.svb; }
jointsheet() { { printf "name\trole\tpath\ngermline\tnormal\t%s\n" $(simreads $2 germline); for t in $TUMORS; do printf "%s\ttumor\t%s\n" $t $(simreads $3 $t); done; } > $1; }
DEV=$(basename $(findmnt -no SOURCE /home) | sed 's/p[0-9]*$//'); HW=$(ls -d /sys/block/$DEV/device/hwmon* 2>/dev/null | head -1)

# janus on one sheet, measured: the run in $D/work with pidstat, iostat, vmstat and temperature monitors, then tools.py measure
# writes run.json, steps.tsv and files.tsv; calls.vcf, the count files and the logs are kept, the rest of work is discarded.
# janus_run EXPERIMENT DIR PAIR SHEET [coverage_x] [extra options of janus run]
janus_run() {
    local E=$1 D=$2 p=$3 sheet=$4 n=$5; shift 5
    [ -s $D/run.json ] && return
    mkdir -p $D/work $D/monitor $D/logs; log "[$(basename $D)] start"
    for c in $(awk 'NR > 1 {print $3}' $sheet) $D/work/n.svb $D/work/t.svb; do [ -e $c ] && dd if=$c iflag=nocache count=0 status=none; done
    sync; sleep 30
    python3 -c "import json; print(json.dumps(dict(experiment='$E', pair='$p', coverage_x=${n:-None}, bases_per_x=$G, threads=$TH, intermediates='deleted', samples={l.split()[0]: l.split()[2] for l in open('$sheet').read().split(chr(10))[1:] if l})))" > $D/monitor/params.json
    pidstat -h -r -u -d -C 'janus|minimap2' 2 >> $D/monitor/pidstat.txt 2>&1 & local m1=$!
    iostat -dxm -t $DEV 5 >> $D/monitor/iostat.txt 2>&1 & local m2=$!
    vmstat -t 5 >> $D/monitor/vmstat.txt 2>&1 & local m3=$!
    ( while :; do echo "$(date +%s) $(cat $HW/temp1_input) $(cat $HW/temp3_input) $(sensors 2>/dev/null | awk '/Tctl/{print $2}')"; sleep 5; done ) >> $D/monitor/temp.txt & local m4=$!
    date +%s.%N > $D/start.txt
    /usr/bin/time -f "total %e s wall %M KB peak %U s user %S s sys %I fs-in-blocks %O fs-out-blocks %x exit" \
        taskset -c $CPUS $J run --samples $sheet --ref $REF/${REFNAME:-hg38_parmask} --out $D/work --threads $TH "$@" > $D/pipeline.log 2> $D/total.txt
    kill $m1 $m2 $m3 $m4 2>/dev/null; wait $m1 $m2 $m3 $m4 2>/dev/null
    [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ] || { log "[$(basename $D)] FAILED: $(tail -1 $D/total.txt)"; return; }
    $TOOLS measure $D >> $Q 2>&1 || { log "[$(basename $D)] measure failed, intermediates left in place"; return; }
    mv $D/work/calls.vcf $D/work/*.cnt $D/work/*.cnt.hist $D/work/*.svb.exc $D/work/steps.tsv $D/ 2>/dev/null; mv $D/work/log.* $D/logs/
    discard $D/work; log "[$(basename $D)] done"
}

# another tool, measured in a cgroup slice (memory limit, peak memory of the container) with /usr/bin/time: measured DIR PARAMS_JSON MEM command...
measured() {
    local D=$1 P=$2 mem=$3; shift 3
    [ -s $D/run.json ] && return 0
    mkdir -p $D/out; echo "$P" > $D/params.json
    local sl="sb_$(basename $(dirname $D) | tr -c 'a-z0-9\n' '_')_$(basename $D | tr -c 'a-z0-9\n' '_').slice" c pid
    log "[$(basename $(dirname $D))/$(basename $D)] start"
    CGPARENT=$sl DMEM=$mem systemd-run --user --scope --quiet --slice=$sl -p MemoryMax=$mem -p MemorySwapMax=0 \
        /usr/bin/time -f "%e %M %U %S %x" -o $D/time.txt taskset -c $CPUS "$@" > $D/out/stdout.txt 2>&1 &
    pid=$!; c=/sys/fs/cgroup/$sl
    ( m=0; while kill -0 $pid 2> /dev/null; do x=$(awk '/^anon /{a=$2} /^shmem /{s=$2} END{print a+s}' $c/memory.stat 2> /dev/null); [ "${x:-0}" -gt $m ] && { m=$x; echo $m > $D/anon_peak.txt; }; sleep 2; done ) &
    wait $pid
    [ -d $c ] && { echo "memory_peak $(cat $c/memory.peak)"; [ -s $D/anon_peak.txt ] && echo "anon_peak $(cat $D/anon_peak.txt)"; } > $D/cgroup.txt
    $TOOLS tool $D >> $Q 2>&1 || log "[$(basename $D)] no measurement"
    for f in calls.vcf log.txt stdout.txt; do [ -e $D/out/$f ] && mv $D/out/$f $D/; done
    discard $D/out
}

# a container of a caller: the repository, the reads and the output directory mounted, the CPUs and memory limited
dk() {
    local img=$1 m=() d; shift; mkdir -p $O/.tmp
    for d in $S $(printf "%s\n" "$(dirname $T)" "$(dirname $N)" "$O" | grep -v "^$S\(/\|$\)" | sort -u); do m+=(-v "$d:$d"); done
    docker run --rm --log-driver none --cpuset-cpus $CPUS --memory ${DMEM:-40g} --memory-swap ${DMEM:-40g} --oom-score-adj 1000 ${CGPARENT:+--cgroup-parent=$CGPARENT} \
        -e OMP_NUM_THREADS=$P -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$O/.tmp:/tmp" -w "$O" "${m[@]}" "$img" "$@"
}

# one caller on one pair: caller TOOL TUMOR.bam NORMAL.bam OUTDIR THREADS; the PASS calls end in OUTDIR/calls.vcf
caller() {
    local tool=$1 T=$(readlink -f $2) N=$(readlink -f $3) P=$5; mkdir -p $4; local O=$(readlink -f $4)
    [ $P -gt 16 ] && P=16
    local FA=$REF/hg38_parmask.fa VNTR=$REF/hg38_parmask.trf.bed
    exec > $O/log.txt 2>&1
    echo "[$(date '+%F %T')] $tool tumor=$T normal=$N out=$O threads=$P"
    pass() { awk -F'\t' '/^#/ || $7=="PASS" || $7=="."' "$1" > "$2"; }
    sub() { $TOOLS subtract "$@"; }
    case $tool in
    severus_hp) $L/severus/bin/severus --target-bam $T --control-bam $N --phasing-vcf $(dirname $N)/clair3/phased_merge_output.vcf.gz --out-dir $O/severus --vntr-bed $VNTR -t $P
        pass $O/severus/somatic_SVs/severus_somatic.vcf $O/calls.vcf ;;
    savana) $L/savana/bin/python $H/tools.py savana --tumour $T --normal $N --ref $FA --outdir $O/savana --threads $P --pb --sample s --overwrite
        pass $O/savana/s.classified.somatic.vcf $O/calls.vcf ;;
    nanomonsv) local M=$L/nanomonsv/bin/nanomonsv
        sort -k1,1 -k2,2n $VNTR | bgzip > $O/repeat.bed.gz && tabix -p bed $O/repeat.bed.gz; mkdir -p $O/t $O/n
        $M parse $T $O/t/t --reference_fasta $FA; $M parse $N $O/n/n --reference_fasta $FA
        $M get $O/t/t $T $FA --control_prefix $O/n/n --control_bam $N --processes $P --qv25 --simple_repeat_bed $O/repeat.bed.gz
        pass $O/t/t.nanomonsv.result.vcf $O/calls.vcf ;;
    sniffles2) local F=$L/sniffles/bin/sniffles
        $F --input $N --snf $O/n.snf --threads $P --reference $FA --tandem-repeats $VNTR --allow-overwrite
        $F --input $T --snf $O/t.snf --threads $P --reference $FA --tandem-repeats $VNTR --allow-overwrite
        $F --input $O/n.snf $O/t.snf --vcf $O/joint.vcf --threads $P --reference $FA --allow-overwrite
        awk -F'\t' '/^#/ || ($7=="PASS" && $8 ~ /SUPP_VEC=01(;|$)/)' $O/joint.vcf > $O/calls.vcf ;;
    sniffles2_mosaic) local F=$L/sniffles/bin/sniffles
        $F --input $T --vcf $O/t.mosaic.vcf --mosaic --threads $P --reference $FA --tandem-repeats $VNTR --allow-overwrite
        $F --input $N --vcf $O/n.vcf --threads $P --reference $FA --tandem-repeats $VNTR --allow-overwrite --minsupport 1
        sub $O/t.mosaic.vcf $O/n.vcf --dist 1000 > $O/calls.vcf ;;
    svision-pro) local I=janus-bench/svision-pro:2.5-pysam0.21 RR=$K/ref/svision-pro NP=$(( P / 2 )); [ $NP -lt 1 ] && NP=1; mkdir -p $O/svp
        dk $I SVision-pro --target_path $T --base_path $N --genome_path $FA --model_path $RR/src/pre_process/model_liteunet_256_8_16_32_32_32.pth \
           --access_path $RR/src/pre_process/hg38.access.10M.bed --out_path $O/svp --sample_name tumor --preset hifi --detect_mode somatic --process_num $NP --unet_cpu_num 2 --device cpu
        local v=$(ls $O/svp/*.vcf | grep -v somatic_s | head -1)
        dk $I python $RR/extract_op.py --input_vcf $v --extract somatic; pass ${v%.vcf}.somatic_s1.vcf $O/calls.vcf ;;
    cutesv) local I=quay.io/biocontainers/cutesv:2.1.4--pyhdfd78af_1 HI="--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"
        mkdir -p $O/t.work $O/n.work
        dk $I cuteSV $T $FA $O/t.vcf $O/t.work $HI -l 50 -S tumor -t $P; dk $I cuteSV $N $FA $O/n.vcf $O/n.work $HI -l 50 -s 2 -S normal -t $P
        sub $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf ;;
    debreak) local I=quay.io/biocontainers/debreak:1.3--h9ee0642_0
        dk $I debreak --bam $T -o $O/t -t $P --min_size 50; dk $I debreak --bam $N -o $O/n -t $P --min_size 50 --min_support 2
        sub $O/t/debreak.vcf $O/n/debreak.vcf --dist 500 > $O/calls.vcf ;;
    dysgu) local I=quay.io/biocontainers/dysgu:1.8.8--py313h5d164f8_0
        dk $I dysgu call --mode pacbio-revio --diploid False -p $P -x -o $O/t.vcf $FA $O/tmp_t $T
        dk $I dysgu call --mode pacbio-revio -p $P -x -o $O/n.vcf $FA $O/tmp_n $N
        dk $I dysgu filter -p $P --normal-vcf $O/n.vcf -o $O/somatic.vcf $O/t.vcf $N
        bcftools view -f PASS -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="TRA" || abs(INFO/SVLEN)>=50' $O/somatic.vcf -Ov -o $O/calls.vcf ;;
    kled) local I=quay.io/biocontainers/kled:1.2.11--he5eddf3_0
        dk $I kled -R $FA --CCS -t $P $T > $O/t.vcf; dk $I kled -R $FA --CCS -t $P --F2 $N > $O/n.vcf
        sub $O/t.vcf $O/n.vcf --dist 1000 > $O/calls.vcf ;;
    svdss) local I=quay.io/biocontainers/svdss:2.1.1--he17396a_0 FMD=$K/ref/svdss/hg38_parmask.fmd
        one() { dk $I SVDSS smooth --reference $FA --bam $1 --threads $P > $O/$2.smoothed.bam; samtools index -@ $P $O/$2.smoothed.bam
            dk $I SVDSS search --index $FMD --bam $O/$2.smoothed.bam --threads $P > $O/$2.sfs.txt
            dk $I SVDSS call --reference $FA --bam $O/$2.smoothed.bam --sfs $O/$2.sfs.txt --threads $P --min-sv-length 50 --min-cluster-weight $3 > $O/$2.vcf; }
        one $T t 4; one $N n 2; sub $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf ;;
    delly) local I=quay.io/biocontainers/delly:2.6.0--h3752d28_0
        sm() { samtools view -H $1 | awk -F'\t' '/^@RG/{for(i=2;i<=NF;i++) if($i~/^SM:/){print substr($i,4); exit}}'; }
        cp $K/ref/human.hg38.excl.tsv $O/excl.tsv; printf "%s\ttumor\n%s\tcontrol\n" "$(sm $T)" "$(sm $N)" > $O/samples.tsv
        dk $I delly lr -y pb -g $FA -x $O/excl.tsv -o $O/joint.bcf $T $N
        dk $I delly filter -f somatic -o $O/somatic.bcf -s $O/samples.tsv $O/joint.bcf
        bcftools view -f PASS $O/somatic.bcf -Ov -o $O/calls.vcf ;;
    *) echo "unknown caller $tool"; exit 2 ;;
    esac
    echo "[$(date '+%F %T')] done: $(grep -vc '^#' $O/calls.vcf) calls in $O/calls.vcf"
}
