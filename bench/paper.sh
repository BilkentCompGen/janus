#!/usr/bin/env bash
# usage: bench/paper.sh <step> [pairs or tools]
#   ref          LCP core and minimap2 indexes of the references (hg38 with chrY PAR masked, CHM13 v2), SVDSS index, Delly exclusion list, SVision-pro model
#   fetch        download the read sets of the paper's pairs into hifi/<pair>
#   align        minimap2 of both read sets of each CASTLE pair, measured (runs/align)
#   full         janus on each pair, measured (runs/full)
#   truth        the published CASTLE consensus regenerated with Minda, the training labels of the labels step (hifi/truth/castle_without_snape)
#   coverage     janus at 5-30x on every pair, then verdict
#   purity       janus on tumor/normal mixtures at 30x, then verdict
#   verdict      one classifier over the full, coverage and purity runs, each CASTLE pair held out, trained on the consensus of the 17 call sets other than janus (verdict/; builds that consensus first if missing)
#   phase        Clair3 + LongPhase on each normal, haplotag both BAMs, measured (runs/phase)
#   competitors  every other caller on the CASTLE pairs, measured (runs/competitors/<tool>/<pair>); competitors <tool> hg008 for the other pairs
#   covcallers   Severus, SAVANA and nanomonsv on the same reads as the coverage and purity runs of janus (runs/coverage_callers)
#   consensus    Minda consensus of the CASTLE call sets with and without janus, every call set scored against both (eval_sets)
#   loo          every HiFi call set scored against a consensus built without it (eval_loo)
#   external     janus and every HiFi caller on HG008 and COLO829 against their truth sets (runs/external_eval)
#   covscore     janus, Severus, SAVANA and nanomonsv at every depth and purity (runs/coverage_eval)
#   labels       the classifier trained on the published CASTLE consensus of 11 call sets, scored (supp/labels11)
#   supp         held-out AUC, all 3-train/2-test splits, the classifier ablation and the feature groups (supp)
#   grid         full grids of the gradient-boosted trees (HistGradientBoosting) and of XGBoost, 3^5 parameter sets each, every set run like the paper's classifier with each CASTLE pair held out (supp/grid.tsv, supp/xgb.tsv)
#   sim          the 20-tumor simulation at five normal/tumor depths: pairwise and joint runs, classified and scored (runs/sim)
#   simmap       minimap2 on the 21 simulated read sets at 5/5 and 20/20, one at a time, measured (runs/sim/mm2)
#   sweep        the count-stage thresholds theta0 (5, 10, 15, 20, 25) x theta_t (20, 30, 40, 50, 60), plus (70, 40), (10, 5) and (70, 5) for m0 = 2 or mt = 3 on every run, on the CASTLE full, coverage and purity runs: janus reruns only a minimum read support pair (m0, mt) not seen before, one classifier and Minda score per distinct configuration (runs/sweep)
#   sweepscore   each sweep configuration scored like the paper: consensus with and without janus, leave-one-out, every depth and purity (runs/sweep/conf/<config>/rescore.tsv, runs/sweep/rescore.txt)
#   timing       the CASTLE pairs and the simulated runs once more, one at a time on an otherwise idle machine, measured (runs_clean)
#   figures      Figure 1, the other figures, the tables and the quoted numbers (bench/figures.py, bench/out), and the paper
#   guard        memory guard, disk guard and load log (run as services)
#   caller       one caller on one pair (used by competitors): caller <tool> <tumor.bam> <normal.bam> <out dir> <threads>
S=/home/yuio/janus; H=$S/bench; R=$H/runs; K=$H/callers; L=$S/lib/callers; G=3100000000; V=$S/.venv-ml/bin/python
CASTLE="hcc1954 h2009 hcc1937 h1437 hcc1395"; LABELED="$CASTLE hg008 colo829"
TOOLS="severus_hp savana nanomonsv sniffles2 sniffles2_mosaic svision-pro cutesv debreak dysgu kled svdss delly"
HIFI="severus_phased savana_phased nanomonsv sniffles2 dysgu svision-pro cutesv debreak kled svdss"; M=$HOME/.local/bin/minda
X=${1:?experiment, see the header of this script}; shift
Q=$R/queue.log; mkdir -p $R
log() { echo "$(date '+%F %T') [$X] $*" >> $Q; }
py() { python3 $H/paper.py "$@"; }
read -r -d '' SVB2FQ <<'CEOF'
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    FILE* f = fopen(argv[1], "rb");
    const char* name = argc > 2 ? argv[2] : "r";
    char hdr[32];
    if (!f || fread(hdr, 1, 32, f) != 32 || memcmp(hdr, "SVB1", 4)) { fprintf(stderr, "not an SVB1 file\n"); return 1; }
    static char out[1 << 20];
    setvbuf(stdout, out, _IOFBF, sizeof out);
    const char B[4] = {'A', 'C', 'G', 'T'};
    size_t cap = 1 << 20;
    uint8_t* p = malloc(cap);
    char* s = malloc(cap * 4 + 1);
    char* q = malloc(cap * 4 + 1);
    uint64_t n = 0;
    uint32_t tag;
    while (fread(&tag, 4, 1, f) == 1) {
        uint32_t L = tag & 0x7fffffffu, raw = tag >> 31, nb = raw ? L : (L + 3) >> 2;
        if (nb > cap || L > cap * 4) {
            cap = (nb > L ? nb : L) * 2;
            p = realloc(p, cap); s = realloc(s, cap * 4 + 1); q = realloc(q, cap * 4 + 1);
        }
        if (fread(p, 1, nb, f) != nb) { fprintf(stderr, "truncated at read %lu\n", (unsigned long)n); return 1; }
        if (raw) memcpy(s, p, L);
        else for (uint32_t i = 0; i < L; i++) s[i] = B[(p[i >> 2] >> (6 - 2 * (i & 3))) & 3];
        char* r = s; uint32_t n = L;
        const char* A = "AATGTACTTCGTTCAGTTACGTATTGCT";
        const char* E = "GCAATACGTAACTGAACGAAGTACATT";
        uint32_t cut = 0;
        for (uint32_t j = 0; j + 21 <= n && j < 100; j++) if (!memcmp(r + j, A + 7, 21)) { cut = j + 21; break; }
        if (!cut) for (uint32_t k = 20; k >= 8; k--) if (k <= n && !memcmp(r, A + 28 - k, k)) { cut = k; break; }
        r += cut; n -= cut;
        uint32_t end = n;
        for (uint32_t j = n > 100 ? n - 100 : 0; j + 21 <= n; j++) if (!memcmp(r + j, E, 21)) { end = j; break; }
        if (end == n) for (uint32_t k = 20; k >= 8; k--) if (k <= n && !memcmp(r + n - k, E, k)) { end = n - k; break; }
        n = end;
        if (n < 100) continue;
        memmove(s, r, n); L = n;
        memset(q, 'I', L);
        printf("@%s_%lu\n", name, (unsigned long)n++);
        fwrite(s, 1, L, stdout); fputs("\n+\n", stdout);
        fwrite(q, 1, L, stdout); fputc('\n', stdout);
    }
    fprintf(stderr, "%lu reads\n", (unsigned long)n);
    return 0;
}
CEOF
src() { [ -s $S/hifi/$1/$2.bam ] && echo $S/hifi/$1/$2.bam || echo $S/hifi/$1/$2.fastq.gz; }
castle() {
    local P=${1^^} SA=${1^^} GR=${1^^}_subsample.bam.gripss.filtered_vaf.vcf; CD=$S/hifi/castle-calls/variant_calls_and_benchmarks/$P/$P
    [ $1 = hcc1395 ] && { SA=H1395; GR=H1395_realigned.bam.gripss.filtered_vaf.vcf; }
    PH="$CD/severus_HIFI_vaf.vcf $CD/savana_$SA.haplotagged.classified.somatic_HIFI_vaf.vcf $CD/tumor.nanomonsv.result_HIFI_vaf.vcf $CD/sniffles2_HIFI_vaf.vcf"
    ONT="$CD/severus_ont_vaf.vcf $CD/savana_$SA.haplotagged.classified.somatic_ONT_vaf.vcf $CD/tumor.nanomonsv.result_ONT_vaf.vcf $CD/sniffles2_ONT_vaf.vcf"
    SR="$CD/svaba_merged.vcf $CD/$GR $CD/manta_somaticSV_vaf.vcf"
}
norm() { case $2 in severus*|savana*|nanomonsv|sniffles2) cp $1 $3 ;; *) py norm $1 $2 > $3 ;; esac; }
ens() { local o=$1 k=$2; shift 2; $M ensemble --vcfs "$@" --conditions "[[[caller_names[:$k], caller_names[$k:$((k + 4))], caller_names[$((k + 4)):]], '>=', 2],'&', [caller_names,'>=', 6]]" --min_size 50 --tolerance 500 --out_dir $o > $o.log 2>&1; }
tset() { $M truthset --base $1 --vcfs $2 --min_size 50 --tolerance 500 --out_dir $3 > $3.log 2>&1; }
depth() { python3 -c "import json; r=json.load(open('$R/full/$1/run.json'))['reads']; print(int(min(r['n']['bases'], r['t']['bases']) // $G))"; }
DEV=$(basename $(findmnt -no SOURCE /home) | sed 's/p[0-9]*$//'); HW=$(ls -d /sys/block/$DEV/device/hwmon* 2>/dev/null | head -1)

janus_run() {
    local E=$1 D=$2 p=$3 n=$4 th=16
    [ -s $D/run.json ] && return
    mkdir -p $D/work $D/monitor $D/logs; log "[$(basename $D)] start"
    for c in $(src $p n) $(src $p t) $D/work/n.svb $D/work/t.svb; do [ -e $c ] && dd if=$c iflag=nocache count=0 status=none; done
    sync; sleep 30
    printf '{"experiment": "%s", "pair": "%s", "coverage_x": %s, "bases_per_x": %s, "threads": %s, "intermediates": "deleted"}\n' \
        $E $p "${n:-null}" $G $th > $D/monitor/params.json
    pidstat -h -r -u -d -C 'janus|minimap2' 2 >> $D/monitor/pidstat.txt 2>&1 & m1=$!
    iostat -dxm -t $DEV 5 >> $D/monitor/iostat.txt 2>&1 & m2=$!
    vmstat -t 5 >> $D/monitor/vmstat.txt 2>&1 & m3=$!
    ( while :; do echo "$(date +%s) $(cat $HW/temp1_input) $(cat $HW/temp3_input) $(sensors 2>/dev/null | awk '/Tctl/{print $2}')"; sleep 5; done ) >> $D/monitor/temp.txt & m4=$!
    date +%s.%N > $D/start.txt
    /usr/bin/time -f "total %e s wall %M KB peak %U s user %S s sys %I fs-in-blocks %O fs-out-blocks %x exit" \
        taskset -c 0-$((th - 1)) bash $S/pipeline.sh $D/work $(src $p n) $(src $p t) $th 2> $D/total.txt
    kill $m1 $m2 $m3 $m4 2>/dev/null; wait $m1 $m2 $m3 $m4 2>/dev/null
    [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ] || { log "[$(basename $D)] FAILED: $(tail -1 $D/total.txt)"; return; }
    py measure $D >> $Q 2>&1 || { log "[$(basename $D)] measure failed, intermediates left in place"; return; }
    mv $D/work/calls.vcf $D/work/n.cnt $D/work/t.cnt $D/work/t.cnt.hist $D/work/n.svb.exc $D/work/t.svb.exc $D/; mv $D/work/log.* $D/logs/
    rm -rf $D/work; log "[$(basename $D)] done"
}

measured() {
    local D=$1 P=$2 mem=$3; shift 3
    [ -s $D/run.json ] && return 0
    mkdir -p $D/out; echo "$P" > $D/params.json
    local sl="sb_$(basename $(dirname $D) | tr -c 'a-z0-9\n' '_')_$(basename $D | tr -c 'a-z0-9\n' '_').slice" c pid
    log "[$(basename $(dirname $D))/$(basename $D)] start"
    CGPARENT=$sl DMEM=$mem systemd-run --user --scope --quiet --slice=$sl -p MemoryMax=$mem -p MemorySwapMax=0 \
        /usr/bin/time -f "%e %M %U %S %x" -o $D/time.txt taskset -c 0-15 "$@" > $D/out/stdout.txt 2>&1 &
    pid=$!; c=/sys/fs/cgroup/$sl
    ( m=0; while kill -0 $pid 2> /dev/null; do x=$(awk '/^anon /{a=$2} /^shmem /{s=$2} END{print a+s}' $c/memory.stat 2> /dev/null); [ "${x:-0}" -gt $m ] && { m=$x; echo $m > $D/anon_peak.txt; }; sleep 2; done ) &
    wait $pid
    [ -d $c ] && { echo "memory_peak $(cat $c/memory.peak)"; [ -s $D/anon_peak.txt ] && echo "anon_peak $(cat $D/anon_peak.txt)"; } > $D/cgroup.txt
    py tool $D >> $Q 2>&1 || log "[$(basename $D)] no measurement"
    for f in calls.vcf log.txt stdout.txt; do [ -e $D/out/$f ] && mv $D/out/$f $D/; done
    rm -rf $D/out
}

aligned() {
    local p=$1 s=$2 B=$3 ref=$4
    minimap2 -ax map-hifi -t 16 -Y --secondary=no -R "@RG\tID:${p}_$s\tSM:${p}_$s" $ref $(src $p $s) 2> $B.minimap2.log \
        | samtools sort -@ 4 -m 2G -T $B.st -o $B.tmp - && mv $B.tmp $B && samtools index -@ 8 $B
}

caller() {
    local tool=$1 T=$(readlink -f $2) N=$(readlink -f $3) P=$5; mkdir -p $4; local O=$(readlink -f $4)
    [ $P -gt 16 ] && P=16
    local REF=${REF:-$S/ref/hg38_parmask.fa} VNTR=${VNTR:-$S/ref/vntr.trf.bed} CPUS=${CPUS:-0-15}
    exec > $O/log.txt 2>&1
    echo "[$(date '+%F %T')] $tool tumor=$T normal=$N out=$O threads=$P"
    dk() {
        local img=$1 m=() d; shift; mkdir -p $O/.tmp
        for d in $S $(printf "%s\n" "$(dirname $T)" "$(dirname $N)" "$O" | grep -v "^$S\(/\|$\)" | sort -u); do m+=(-v "$d:$d"); done
        docker run --rm --log-driver none --cpuset-cpus $CPUS --memory ${DMEM:-40g} --memory-swap ${DMEM:-40g} --oom-score-adj 1000 ${CGPARENT:+--cgroup-parent=$CGPARENT} \
            -e OMP_NUM_THREADS=$P -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$O/.tmp:/tmp" -w "$O" "${m[@]}" "$img" "$@"
    }
    pass() { awk -F'\t' '/^#/ || $7=="PASS" || $7=="."' "$1" > "$2"; }
    sub() { py subtract "$@"; }
    case $tool in
    severus_hp) $L/severus/bin/severus --target-bam $T --control-bam $N --phasing-vcf $(dirname $N)/clair3/phased_merge_output.vcf.gz --out-dir $O/severus --vntr-bed $VNTR -t $P
        pass $O/severus/somatic_SVs/severus_somatic.vcf $O/calls.vcf ;;
    savana) $L/savana/bin/python $H/paper.py savana --tumour $T --normal $N --ref $REF --outdir $O/savana --threads $P --pb --sample s --overwrite
        pass $O/savana/s.classified.somatic.vcf $O/calls.vcf ;;
    nanomonsv) local M=$L/nanomonsv/bin/nanomonsv
        sort -k1,1 -k2,2n $VNTR | bgzip > $O/repeat.bed.gz && tabix -p bed $O/repeat.bed.gz; mkdir -p $O/t $O/n
        $M parse $T $O/t/t --reference_fasta $REF; $M parse $N $O/n/n --reference_fasta $REF
        $M get $O/t/t $T $REF --control_prefix $O/n/n --control_bam $N --processes $P --qv25 --simple_repeat_bed $O/repeat.bed.gz
        pass $O/t/t.nanomonsv.result.vcf $O/calls.vcf ;;
    sniffles2) local F=$L/sniffles/bin/sniffles
        $F --input $N --snf $O/n.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
        $F --input $T --snf $O/t.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
        $F --input $O/n.snf $O/t.snf --vcf $O/joint.vcf --threads $P --reference $REF --allow-overwrite
        awk -F'\t' '/^#/ || ($7=="PASS" && $8 ~ /SUPP_VEC=01(;|$)/)' $O/joint.vcf > $O/calls.vcf ;;
    sniffles2_mosaic) local F=$L/sniffles/bin/sniffles
        $F --input $T --vcf $O/t.mosaic.vcf --mosaic --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
        $F --input $N --vcf $O/n.vcf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite --minsupport 1
        sub $O/t.mosaic.vcf $O/n.vcf --dist 1000 > $O/calls.vcf ;;
    svision-pro) local I=janus-bench/svision-pro:2.5-pysam0.21 RR=$K/ref/svision-pro NP=$(( P / 2 )); [ $NP -lt 1 ] && NP=1; mkdir -p $O/svp
        dk $I SVision-pro --target_path $T --base_path $N --genome_path $REF --model_path $RR/src/pre_process/model_liteunet_256_8_16_32_32_32.pth \
           --access_path $RR/src/pre_process/hg38.access.10M.bed --out_path $O/svp --sample_name tumor --preset hifi --detect_mode somatic \
           --process_num $NP --unet_cpu_num 2 --device cpu
        local v=$(ls $O/svp/*.vcf | grep -v somatic_s | head -1)
        dk $I python $RR/extract_op.py --input_vcf $v --extract somatic; pass ${v%.vcf}.somatic_s1.vcf $O/calls.vcf ;;
    cutesv) local I=quay.io/biocontainers/cutesv:2.1.4--pyhdfd78af_1 HI="--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"
        mkdir -p $O/t.work $O/n.work
        dk $I cuteSV $T $REF $O/t.vcf $O/t.work $HI -l 50 -S tumor -t $P; dk $I cuteSV $N $REF $O/n.vcf $O/n.work $HI -l 50 -s 2 -S normal -t $P
        sub $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf ;;
    debreak) local I=quay.io/biocontainers/debreak:1.3--h9ee0642_0
        dk $I debreak --bam $T -o $O/t -t $P --min_size 50; dk $I debreak --bam $N -o $O/n -t $P --min_size 50 --min_support 2
        sub $O/t/debreak.vcf $O/n/debreak.vcf --dist 500 > $O/calls.vcf ;;
    dysgu) local I=quay.io/biocontainers/dysgu:1.8.8--py313h5d164f8_0
        dk $I dysgu call --mode pacbio-revio --diploid False -p $P -x -o $O/t.vcf $REF $O/tmp_t $T
        dk $I dysgu call --mode pacbio-revio -p $P -x -o $O/n.vcf $REF $O/tmp_n $N
        dk $I dysgu filter -p $P --normal-vcf $O/n.vcf -o $O/somatic.vcf $O/t.vcf $N
        bcftools view -f PASS -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="TRA" || abs(INFO/SVLEN)>=50' $O/somatic.vcf -Ov -o $O/calls.vcf ;;
    kled) local I=quay.io/biocontainers/kled:1.2.11--he5eddf3_0
        dk $I kled -R $REF --CCS -t $P $T > $O/t.vcf; dk $I kled -R $REF --CCS -t $P --F2 $N > $O/n.vcf
        sub $O/t.vcf $O/n.vcf --dist 1000 > $O/calls.vcf ;;
    svdss) local I=quay.io/biocontainers/svdss:2.1.1--he17396a_0 FMD=$K/ref/svdss/hg38_parmask.fmd
        one() { dk $I SVDSS smooth --reference $REF --bam $1 --threads $P > $O/$2.smoothed.bam; samtools index -@ $P $O/$2.smoothed.bam
            dk $I SVDSS search --index $FMD --bam $O/$2.smoothed.bam --threads $P > $O/$2.sfs.txt
            dk $I SVDSS call --reference $REF --bam $O/$2.smoothed.bam --sfs $O/$2.sfs.txt --threads $P --min-sv-length 50 --min-cluster-weight $3 > $O/$2.vcf; }
        one $T t 4; one $N n 2; sub $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf ;;
    delly) local I=quay.io/biocontainers/delly:2.6.0--h3752d28_0
        sm() { samtools view -H $1 | awk -F'\t' '/^@RG/{for(i=2;i<=NF;i++) if($i~/^SM:/){print substr($i,4); exit}}'; }
        cp ${EXCL:-$K/ref/human.hg38.excl.tsv} $O/excl.tsv; printf "%s\ttumor\n%s\tcontrol\n" "$(sm $T)" "$(sm $N)" > $O/samples.tsv
        dk $I delly lr -y pb -g $REF -x $O/excl.tsv -o $O/joint.bcf $T $N
        dk $I delly filter -f somatic -o $O/somatic.bcf -s $O/samples.tsv $O/joint.bcf
        bcftools view -f PASS $O/somatic.bcf -Ov -o $O/calls.vcf ;;
    *) echo "unknown caller $tool"; exit 2 ;;
    esac
    echo "[$(date '+%F %T')] done: $(grep -vc '^#' $O/calls.vcf) calls in $O/calls.vcf"
}

simreads() {
    local SP=$S/simulation/phylo/reads O=$R/sim B=$S/build/janus TUM=$(seq -f "T%02g" 1 20 | tr '\n' ' ') x s t
    for x in 5 10; do for s in germline $TUM; do [ -s $O/x$x/reads/$s.svb.idx ] || { mkdir -p $O/x$x/reads; py thin $SP/$s/$s.svb $O/x$x/reads/$s.svb $((20 / x)) >> $O/x$x/sub.log; }; done; done
    [ -x $O/svb2fq ] || { echo "$SVB2FQ" > $O/svb2fq.c; gcc -O2 -o $O/svb2fq $O/svb2fq.c; }
    mkdir -p $O/reads20; [ -s $O/reads20/germline.svb.idx ] || { $O/svb2fq $SP/germline/germline.svb germline 2> $O/reads20/germline.log | pigz -p 16 -3 > $O/reads20/germline.fastq.gz && $B pack -i $O/reads20/germline.fastq.gz -o $O/reads20/germline.svb -t 16; rm -f $O/reads20/germline.fastq.gz; }
    for t in $TUM; do for e in svb svb.idx; do ln -sf $SP/$t/$t.$e $O/reads20/$t.$e; done; done
}
reads() { [ $1 = 20 ] && echo $R/sim/reads20/$2 || echo $R/sim/x$1/reads/$2; }
joint() {
    local nx=$1 tx=$2 B=$S/build/janus RF=$S/ref/chm13v2 TUM=$(seq -f "T%02g" 1 20 | tr '\n' ' '); local d=${3:-$R/sim/settings/${1}_${2}/joint}; [ -s $d/calls.vcf ] && return; mkdir -p $d/work; cd $d/work
    st() { n=$1; shift; [ -e done.$n ] && return; /usr/bin/time -f "$n %e s %M KB %U user %S sys" "$@" > log.$n 2>&1 || exit 1; tail -1 log.$n >> steps.txt; touch done.$n; }
    local SM="germline $TUM" s; local svb=$(for s in $SM; do printf "%s.svb," $s; done | sed 's/,$//') cnt=$(for s in $SM; do printf "%s.cnt," $s; done | sed 's/,$//')
    for e in svb svb.idx; do ln -sf $(readlink -f $(reads $nx germline.$e)) germline.$e; for s in $TUM; do ln -sf $(readlink -f $(reads $tx $s.$e)) $s.$e; done; done
    for s in $SM; do st extract_$s $B extract -i $s.svb -o $s.svb.exc -t 16; done
    st count_germline $B count -i germline.svb -o germline.cnt -t 16 -m 0
    for s in $TUM; do st count_$s $B count -i $s.svb -o $s.cnt -t 16 -d germline.cnt; done
    local M=$(python3 -c "import sys; print(min(open(f'{t}.cnt', 'rb').read(27)[26] for t in sys.argv[1:]))" $TUM)
    st graph $B graph -i $cnt -o graph.gph -m $M
    st prune $B prune -i graph.gph -r $svb -o prune -t 16 -m $M
    st select $B select -i graph.gph -d prune.sites -r $svb -o ev -t 16
    st align bash -c "awk '/^>/ { split(substr(\$0, 3), a, \"_\"); f = \"p\" a[1] \".fa\" } { print > f }' ev.fa && seq 1 20 | xargs -P 4 -I{} bash -c '$B place -i p{}.fa -r $RF.cidx -d $RF.fa -o p{} -t 4 > log.place{} 2>&1 && minimap2 -cx map-hifi --secondary=no -t 4 -o local{}.paf p{}.win.fa p{}.on.fa 2> log.local{} && $B lift -i local{}.paf -d $RF.fa -o aln{}.paf > log.lift{} 2>&1' && cat \$(seq -f 'aln%g.paf' 1 20) > aln.paf"
    st clip $B clip -i ev -r aln.paf -o ev.gaps.fa
    st gaps minimap2 -cx sr --secondary=no -t 16 -o ev.gaps.raw.paf $RF.sr.mmi ev.gaps.fa
    st patch $B patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf
    st call $B call -i ev -r $svb -d $RF.fa -g graph.gph -o calls.vcf -v $RF.trf.bed -t 16 -m $M
    cp steps.txt calls.vcf $d/
    for s in $TUM; do python3 -c "print(open('$s.cnt', 'rb').read(27)[26])" > $d/$s.minimum; done
    [ -n "$3" ] || { cp germline.cnt germline.svb.exc $d/; for s in $TUM; do cp $s.cnt $s.cnt.hist $s.svb.exc $d/; done; }
    cd $d; [ -s calls.vcf ] && rm -rf work
}

case $X in
ref)
    cd $S/ref
    $S/build/janus refidx -i hg38_parmask.fa -o hg38_parmask.cidx -t 16
    minimap2 -x map-hifi -d hg38_parmask.mmi hg38_parmask.fa; minimap2 -x sr -d hg38_parmask.sr.mmi hg38_parmask.fa
    $S/build/janus refidx -i chm13v2.fa -o chm13v2.cidx -t 16; minimap2 -x sr -d chm13v2.sr.mmi chm13v2.fa
    [ -s $K/ref/human.hg38.excl.tsv ] || curl -sL -o $K/ref/human.hg38.excl.tsv https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv
    [ -d $K/ref/svision-pro ] || { git clone -q https://github.com/songbowang125/SVision-pro.git $K/ref/svision-pro && git -C $K/ref/svision-pro checkout -q v2.5; }
    mkdir -p $K/ref/svdss; docker run --rm --log-driver none -u "$(id -u):$(id -g)" -v $S:$S quay.io/biocontainers/svdss:2.1.1--he17396a_0 \
        SVDSS index --reference $S/ref/hg38_parmask.fa --index $K/ref/svdss/hg38_parmask.fmd --threads 16 && touch $K/ref/svdss/index.done ;;
fetch)
    bash -s <<'FETCHEOF'
#!/usr/bin/env bash
set -uo pipefail
H=/home/yuio/janus/hifi; RAW=$H/raw; T=$H/truth; TK=$H/tools/sratoolkit/bin; cd "$H"; mkdir -p "$RAW/tmp" "$T"
SDL="https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?location=s3.us-east-1&type=src&acc="
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic
PB=https://downloads.pacbcloud.com/public/revio
DL=4; CV=3; RESERVE=$((150 * 1000 * 1000 * 1000))

JOBS="
sra hcc1954 t SRR28305163
sra hcc1954 n SRR28305160
sra h2009 t SRR28305172
sra h2009 n SRR28305169
sra hcc1937 t SRR28305185
sra hcc1937 n SRR28305182
sra h1437 t SRR28305179
sra h1437 n SRR28305175
giab hg008 t $GIAB/HG008/Liss_lab/BCM_Revio_20240313/HG008-T_ubams $GIAB/HG008/Liss_lab/BCM_Revio_20240313/checksum.md5 m84101_240226_212100_s2.hifi_reads.bam m84101_240226_232018_s3.hifi_reads.bam
giab hg008 n $GIAB/HG008/Liss_lab/BCM_Revio_20240313/HG008-N-D_ubams $GIAB/HG008/Liss_lab/BCM_Revio_20240313/checksum.md5 m84059_240304_183205_s3.hifi_reads.bam m84059_240304_203144_s4.hifi_reads.bam
pacb hcc1395 t $PB/2023Q2/HCC1395/HCC1395/m84039_230414_235240_s2.hifi_reads.default.bam $PB/2023Q2/HCC1395/HCC1395/m84039_230415_002321_s3.hifi_reads.default.bam
pacb hcc1395 n $PB/2023Q2/HCC1395/HCC1395-BL/m84039_230412_183900_s2.hifi_reads.bc2002.bam $PB/2023Q2/HCC1395/HCC1395-BL/m84039_230415_005427_s4.hifi_reads.bc2002.bam
pacb colo829 t $PB/2023Q2/COLO829/COLO829/m84039_230312_025934_s1.hifi_reads.bc2026.bam $PB/2023Q2/COLO829/COLO829/m84039_230328_000836_s3.hifi_reads.bc2025.bam
pacb colo829 n $PB/2023Q2/COLO829/COLO829-BL/m84039_230327_230708_s1.hifi_reads.bc2007.bam $PB/2023Q2/COLO829/COLO829-BL/m84039_230327_233730_s2.hifi_reads.bc2007.bam
"

free() { df --output=avail -B1 "$H" | tail -1; }
log() { echo "$(date '+%F %T') $*"; }
slot() {
    local k=$1 n=$2 i fd; [ $k = dl ] && fd=8 || fd=7
    while :; do for i in $(seq 0 $((n - 1))); do eval "exec $fd> $RAW/tmp/$k.$i.lock"; flock -n $fd && return 0; done; sleep 20; done
}
unslot() { [ $1 = dl ] && { flock -u 8; exec 8>&-; } || { flock -u 7; exec 7>&-; }; }

etag_ok() {
    python3 - "$1" "$2" <<'EOF'
import hashlib, math, os, sys
f, tag = sys.argv[1], sys.argv[2].strip('"'); want, n = tag.split('-'); n = int(n); size = os.path.getsize(f)
for mib in (8, 16, 32, 64, 128, 256, 512):
    ps = mib << 20
    if math.ceil(size / ps) != n: continue
    d = b''
    with open(f, 'rb') as h:
        while True:
            b = h.read(ps)
            if not b: break
            d += hashlib.md5(b).digest()
    sys.exit(0 if hashlib.md5(d).hexdigest() == want else 1)
sys.exit(2)
EOF
}

get() {
    local f=$1 url=$2 bytes=$3 md5=$4 tries=0 ck=()
    [ -s "$f" ] && [ ! -e "$f.aria2" ] && return 0
    [ -n "$md5" ] && ck=(--checksum=md5="$md5")
    [ -e "$f.aria2" ] || while [ "$(free)" -lt $(( bytes + RESERVE )) ]; do log "[$f] waiting for disk"; sleep 600; done
    until aria2c -x16 -s16 -k20M --file-allocation=falloc --console-log-level=warn --summary-interval=60 \
          --show-console-readout=false --max-tries=0 --retry-wait=20 --auto-file-renaming=false "${ck[@]}" \
          -d "$(dirname "$f")" -o "$(basename "$f")" "$url" >> "$RAW/tmp/$(basename "$f").dl.log" 2>&1; do
        tries=$((tries + 1)); log "[$f] aria2c failed ($tries), resuming from the control file"
        [ $tries -ge 30 ] && { log "[FATAL] $f: 30 attempts"; return 1; }; sleep 30
    done
}

bams_to_fastq() {
    local p=$1 s=$2 g=$H/$1/$2.fastq.gz n=0 got b; shift 2
    slot cv $CV
    while [ "$(free)" -lt $(( $(du -cB1 "$@" | tail -1 | cut -f1) * 12 / 10 + RESERVE / 3 )) ]; do log "[$p $s] waiting for disk to convert"; sleep 600; done
    log "[$p $s] converting"
    for b in "$@"; do samtools fastq -@ 4 "$b" 2>> "$RAW/$p/$s.convert.log"; done | pigz -p 8 > "$g.tmp"
    unslot cv
    for b in "$@"; do n=$((n + $(samtools view -c -@ 4 "$b"))); done
    got=$(( $(pigz -dc "$g.tmp" | wc -l) / 4 ))
    if [ "$n" -gt 0 ] && [ "$n" = "$got" ]; then
        mv "$g.tmp" "$g"; rm -f -- "$@" "${@/%/.ok}"; log "[$p $s] done: $got reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $p $s: $n reads in the BAMs, $got in the fastq.gz; BAMs kept"; return 1; fi
}

pacb() {
    local p=$1 s=$2 u f tag bytes files=(); shift 2
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"
    for u in "$@"; do
        f=$RAW/$p/$s.$(basename "$u"); files+=("$f"); [ -e "$f.ok" ] && continue
        tag=$(curl -sSI --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="etag:"{print $2}')
        bytes=$(curl -sSI --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
        slot dl $DL; get "$f" "$u" "$bytes" "" || { unslot dl; return 1; }; unslot dl
        etag_ok "$f" "$tag" && { touch "$f.ok"; log "[$f] downloaded, ETag ok"; } || { log "[FATAL] $f: ETag $tag does not match"; return 1; }
    done
    bams_to_fastq "$p" "$s" "${files[@]}"
}

giab() {
    local p=$1 s=$2 dir=$3 sums m f md5 bytes files=(); sums=$(curl -sS --retry 5 "$4"); shift 4
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"
    for m in "$@"; do
        f=$RAW/$p/$s.$m; files+=("$f")
        md5=$(printf '%s\n' "$sums" | grep -F "$m" | grep -v "\.pbi" | awk '{print $1}' | head -1)
        bytes=$(curl -sSI --retry 5 "$dir/$m" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
        [ -n "$md5" ] && [ -n "$bytes" ] || { log "[FATAL] $p $s $m: no md5/size from GIAB"; return 1; }
        slot dl $DL; get "$f" "$dir/$m" "$bytes" "$md5" || { unslot dl; return 1; }; unslot dl
        log "[$f] downloaded, md5 ok"
    done
    bams_to_fastq "$p" "$s" "${files[@]}"
}

sra() {
    local p=$1 s=$2 acc=$3 f g=$H/$1/$2.fastq.gz meta url md5 bytes r w
    [ -s "$g" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"; f=$RAW/$p/$s.$acc.sra
    meta=$(curl -sS --retry 5 "$SDL$acc")
    url=$(printf '%s' "$meta" | grep -o '"link": "[^"]*"' | head -1 | cut -d'"' -f4)
    md5=$(printf '%s' "$meta" | grep -o '"md5": "[^"]*"' | head -1 | cut -d'"' -f4)
    bytes=$(printf '%s' "$meta" | grep -o '"size": [0-9]*' | head -1 | grep -o '[0-9]*$')
    [ -n "$url" ] && [ -n "$md5" ] && [ -n "$bytes" ] || { log "[FATAL] $acc: no link/md5/size from NCBI"; return 1; }
    slot dl $DL; get "$f" "$url" "$bytes" "$md5" || { unslot dl; return 1; }; unslot dl
    log "[$f] downloaded, md5 ok"
    slot cv $CV
    while [ "$(free)" -lt $(( bytes * 13 / 10 + RESERVE / 3 )) ]; do log "[$p $s] waiting for disk to convert"; sleep 600; done
    log "[$p $s] converting"
    "$TK/fastq-dump" --stdout --split-spot "$f" 2> "$RAW/$p/$s.convert.log" | pigz -p 8 > "$g.tmp"
    unslot cv
    r=$(grep -o '^Read [0-9]* spots' "$RAW/$p/$s.convert.log" | grep -o '[0-9]*')
    w=$(grep -o '^Written [0-9]* spots' "$RAW/$p/$s.convert.log" | grep -o '[0-9]*')
    if [ -n "$r" ] && [ "$r" = "$w" ] && pigz -t "$g.tmp"; then
        mv "$g.tmp" "$g"; rm -f -- "$f"; log "[$p $s] done: $r reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $p $s: conversion check failed (read $r, written $w); $f kept"; return 1; fi
}

zen() {
    mkdir -p "$T/$2"
    curl -s "https://zenodo.org/api/records/$1" | python3 -c "
import json,sys
for f in json.load(sys.stdin)['files']: print(f['key'], f['checksum'].split(':')[1], f['links']['self'])" |
    while read -r k md5 url; do
        [ -s "$T/$2/$k" ] && [ "$(md5sum "$T/$2/$k" | cut -d' ' -f1)" = "$md5" ] && continue
        aria2c -x16 -s16 -k1M --file-allocation=none --console-log-level=warn --summary-interval=0 --max-tries=0 \
               --retry-wait=10 --auto-file-renaming=false --allow-overwrite=true -d "$T/$2" -o "$k" "$url" > /dev/null
        [ "$(md5sum "$T/$2/$k" | cut -d' ' -f1)" = "$md5" ] && log "[truth $2] $k md5 ok" || log "[FATAL] truth $2 $k md5 mismatch"
    done
}
web() {
    local d=$1 u; shift; mkdir -p "$T/$d"
    for u in "$@"; do [ -s "$T/$d/$(basename "$u")" ] || curl -sSfL --retry 5 -o "$T/$d/$(basename "$u")" "$u" || log "[FATAL] truth $d $u"; done
}
truth() {
    zen 5019338 colo829
    web hcc1395 $(for n in 3 4 5 15; do echo "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-022-02816-6/MediaObjects/13059_2022_2816_MOESM${n}_ESM.xlsx"; done)
    log "[truth] done"
}

{
log "########## start ${*:-all} ##########"
[ $# -eq 0 ] && truth &
while read -r fn p rest; do
    [ -z "$fn" ] && continue
    [ $# -gt 0 ] && [[ " $* " != *" $p "* ]] && continue
    $fn $p $rest &
    sleep 3
done <<< "$JOBS"
wait
log "########## end ##########"
} 2>&1 | tee -a "$H/fetch.log"
FETCHEOF
    ;;
align)
    for p in ${*:-$CASTLE}; do for s in n t; do
        measured $R/align/$p/$s "{\"experiment\": \"align\", \"tool\": \"minimap2+samtools\", \"pair\": \"$p\", \"sample\": \"$s\", \"bam\": \"$S/hifi/$p/$s.bam\"}" 60G \
            bash -c "$(declare -f src aligned); S=$S aligned $p $s $S/hifi/$p/$s.bam $S/ref/hg38_parmask.mmi"
    done; done ;;
full)
    for p in ${*:-$LABELED}; do janus_run full $R/full/$p $p ""; done ;;
verdict)
    for p in $CASTLE; do [ -s $H/eval_sets/$p/all_without/None_minda_ensemble.vcf ] || NOJANUS=1 bash $H/paper.sh consensus $p; done
    mkdir -p $H/verdict; ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 $V $H/paper.py verdict > $H/verdict/stdout.txt 2>&1 ); cat $H/verdict/log.txt ;;
truth)
    T=$S/hifi/truth/castle_without_snape; mkdir -p $T
    for p in ${*:-$CASTLE}; do
        castle $p
        $M ensemble --vcfs $PH $ONT $SR --conditions "[[[caller_names[:4], caller_names[4:8], caller_names[8:]], '>=', 2],'&', [caller_names,'>=', 4]]" \
            --min_size 50 --tolerance 500 --out_dir $T/$p > $T/$p.log 2>&1
    done ;;
external)
    O=${O:-$R/external_eval}; mkdir -p $O
    declare -A TB=([hg008]=$S/hifi/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz [colo829]=$S/hifi/truth/colo829/truth.vcf)
    for p in hg008 colo829; do
        for t in snape severus_phased savana_phased nanomonsv sniffles2 sniffles2_mosaic dysgu svision-pro cutesv debreak kled svdss; do
            d=$O/$p/$t; mkdir -p $d
            case $t in
                snape) grep -E '^#|	PASS	' $H/verdict/$p.vcf > $d/in.vcf ;;
                *) v=$R/competitors/$t/$p/calls.vcf; [ -s $v ] || v=$R/competitors/${t%_phased}/$p/calls.vcf; [ -s $v ] || continue
                   norm $v $t $d/in.vcf ;;
            esac
            [ -s $d/minda/results/overall.tsv ] || $M truthset --base ${TB[$p]} --vcfs $d/in.vcf --min_size 50 --tolerance 500 --out_dir $d/minda > $d/minda.log 2>&1
            echo -e "$p\t$t\t$(awk -F'\t' 'NR==2{printf "%s\t%s\t%s\t%.3f\t%.3f\t%.3f", $8, $2, $4, $5, $6, $7}' $d/minda/results/overall.tsv)\t$(py score ${TB[$p]} $d/in.vcf)"
        done
    done > $O/scores.new; mv $O/scores.new $O/scores.tsv; column -t $O/scores.tsv ;;
coverage)
    E=$R/coverage
    for p in ${*:-$CASTLE}; do
        M=$E/mix/$p; mkdir -p $M
        for s in n t; do [ -s $M/$s.svb.idx ] || taskset -c 0-15 $S/build/janus pack -i $(src $p $s) -o $M/$s.svb -t 16 -s 1 > $M/pack_$s.log 2>&1; done
        for n in 5 10 15 20 30; do
            [ $n -le $(depth $p) ] || continue
            D=$E/${p}_${n}x; [ -s $D/run.json ] && continue; mkdir -p $D/work
            for s in n t; do py mix $D/work/$s.svb $M/$s.svb:0:$((n * G)) >> $Q 2>&1 || continue 2; done
            touch $D/work/done.pack_n $D/work/done.pack_t
            janus_run coverage $D $p $n
        done
        rm -rf $M
    done
    bash $H/paper.sh verdict ;;
purity)
    E=$R/purity
    for p in ${*:-hcc1954 hcc1395}; do
        M=$E/mix/$p; mkdir -p $M
        for s in n t; do [ -s $M/$s.svb.idx ] || taskset -c 0-15 $S/build/janus pack -i $(src $p $s) -o $M/$s.svb -t 16 -s 1 > $M/pack_$s.log 2>&1; done
        nb=$(python3 -c "import struct; print(struct.unpack_from('<Q', open('$M/n.svb', 'rb').read(32), 16)[0])")
        nm=$((nb - 24 * G - G / 10)); [ $nm -gt $((30 * G)) ] && nm=$((30 * G))
        [ -s $M/nm.svb.idx ] || py mix $M/nm.svb $M/n.svb:0:$nm >> $Q 2>&1 || continue
        for q in 20 40 60 80 100; do
            [ -s $M/t$q.svb.idx ] || py mix $M/t$q.svb $M/t.svb:0:$((q * 30 * G / 100)) $([ $q -lt 100 ] && echo $M/n.svb:$nm:$(((100 - q) * 30 * G / 100))) >> $Q 2>&1 || continue
            D=$E/${p}_${q}p; mkdir -p $D/work
            for f in svb svb.idx; do ln -sf $M/nm.$f $D/work/n.$f; ln -sf $M/t$q.$f $D/work/t.$f; done
            touch $D/work/done.pack_n $D/work/done.pack_t
            janus_run purity $D $p ""
        done
        rm -rf $M
    done
    bash $H/paper.sh verdict ;;
phase)
    IMG=hkubal/clair3:v2.0.3; RF=$S/ref/hg38_parmask.fa
    for p in ${*:-$CASTLE}; do
        W=$S/hifi/$p; O=$R/phase/$p; T=$O/times.tsv; mkdir -p $O
        [ -s $T ] || echo -e "step\twall_s\tcpu_s\tpeak_gib\texit" > $T
        step() {
            local name=$1 sl=sb_phase_${p}_$1.slice; shift
            grep -q "^$name	" $T && return 0
            CG=$sl /usr/bin/time -f "%e %U %S %M %x" -o $O/$name.time "$@" > $O/$name.log 2>&1
            local c=/sys/fs/cgroup/$sl cpu=0 peak=0
            [ -d $c ] && { cpu=$(awk '/^usage_usec/{print $2/1e6}' $c/cpu.stat); peak=$(awk '/^anon /{print $2}' $c/memory.stat); }
            read wall us sy rss ex < $O/$name.time
            awk -v n=$name -v w=$wall -v c=$cpu -v u=$us -v y=$sy -v m=$peak -v r=$rss -v x=$ex 'BEGIN{r *= 1024; if (r > m) m = r; printf "%s\t%s\t%.1f\t%.2f\t%s\n", n, w, c + u + y, m / 1073741824, x}' >> $T
        }
        dkc() { docker run --rm --log-driver none --cpuset-cpus 0-15 --cgroup-parent=$CG -u "$(id -u):$(id -g)" -v $S:$S $IMG "$@"; }
        export -f dkc; export S IMG
        step clair3 bash -c "dkc /opt/bin/run_clair3.sh --bam_fn=$W/n.bam --ref_fn=$RF --threads=16 --platform=hifi --model_path=/opt/models/hifi_revio --output=$W/clair3 --enable_phasing --longphase_for_phasing"
        for s in n t; do
            [ -s $W/$s.tagged ] && continue
            step haplotag_$s bash -c "dkc /opt/bin/longphase haplotag -s $W/clair3/phased_merge_output.vcf.gz -b $W/$s.bam -r $RF -t 16 --tagSupplementary -o $W/$s.hp"
            samtools view -@ 8 -b -o $W/$s.unmapped.bam $W/$s.bam '*'
            samtools cat -o $W/$s.hpu.bam $W/$s.hp.bam $W/$s.unmapped.bam && mv $W/$s.hpu.bam $W/$s.hp.bam
            step index_$s taskset -c 0-15 samtools index -@ 16 $W/$s.hp.bam
            a=$(samtools idxstats $W/$s.bam | awk '{m+=$3; u+=$4} END{print m, u}'); b=$(samtools idxstats $W/$s.hp.bam | awk '{m+=$3; u+=$4} END{print m, u}')
            [ "$a" = "$b" ] || { log "$p $s: read counts differ ($a vs $b); original kept"; continue; }
            mv $W/$s.hp.bam $W/$s.bam; mv $W/$s.hp.bam.bai $W/$s.bam.bai; rm -f $W/$s.unmapped.bam; echo "$a" > $W/$s.tagged
        done
    done ;;
competitors)
    for t in ${1:-$TOOLS}; do for p in ${2:-$CASTLE}; do
        mem=60G; th=16; [ $t = savana ] && { mem=150G; th=8; }; [ $t = severus_hp ] && [ ! -s $S/hifi/$p/n.tagged ] && continue
        n=$t; [ $t = severus_hp ] && n=severus_phased
        measured $R/competitors/$n/$p "{\"experiment\": \"competitor\", \"tool\": \"$n\", \"pair\": \"$p\", \"threads\": $th}" $mem \
            bash $H/paper.sh caller $t $S/hifi/$p/t.bam $S/hifi/$p/n.bam $R/competitors/$n/$p/out $th
    done; done ;;
sim)
    SP=$S/simulation/phylo/reads; O=$R/sim; E=$O/settings; B=$S/build/janus; RF=$S/ref/chm13v2; TUM=$(seq -f "T%02g" 1 20 | tr '\n' ' ')
    [ -s $O/truth/tree.json ] || py simtruth $O/truth
    simreads
    one() {
        local nx=$1 tx=$2 t=$3; local N=$E/${1}_${2}; local p=$N/$t
        [ -s $p/t.cnt.hist ] && return; mkdir -p $p/work
        for e in svb svb.idx; do ln -sf $(readlink -f $(reads $nx germline.$e)) $p/work/n.$e; ln -sf $(readlink -f $(reads $tx $t.$e)) $p/work/t.$e; done
        touch $p/work/done.pack_n $p/work/done.pack_t
        if [ -s $N/n.svb.exc ] && [ -s $N/n.cnt ]; then
            ln -sf $N/n.svb.exc $p/work/n.svb.exc; ln -sf $N/n.cnt $p/work/n.cnt; touch $p/work/done.extract_n $p/work/done.count_n
        fi
        REF=chm13v2 /usr/bin/time -f "%e %M %U %S %x" -o $p/time.txt bash $S/pipeline.sh $p/work $p/work/n.svb $p/work/t.svb 16 > $p/pipeline.log 2>&1
        cp $p/work/steps.txt $p/work/calls.vcf $p/work/t.cnt.hist $p/work/t.cnt $p/work/t.svb.exc $p/ 2> /dev/null
        for f in n.cnt n.svb.exc; do [ -s $N/$f ] || cp $p/work/$f $N/; done
        [ -s $p/calls.vcf ] && [ -s $p/t.svb.exc ] && rm -rf $p/work
    }
    export -f one reads; export S O E TUM R
    for set in "5 5" "10 10" "20 10" "20 5" "20 20"; do
        set -- $set; N=$E/${1}_${2}
        [ -s $N/n.cnt ] || for t in $TUM; do [ -s $N/$t/t.cnt.hist ] && continue; one $1 $2 $t; break; done
        echo $TUM | tr ' ' '\n' | xargs -P 3 -I{} bash -c "one $1 $2 {}"
        log "pairs $1/$2: $(ls $N/*/t.cnt.hist 2> /dev/null | wc -l) of 20"
        ( joint $1 $2 ); log "joint $1/$2: $(grep -vc '^#' $N/joint/calls.vcf 2> /dev/null) records"
    done
    $V $H/paper.py simjoint | tee $O/scores.txt ;;
covcallers)
    C=$R/coverage_callers; mkdir -p $C/names $C/lists
    sub() { samtools view -@ 12 -b -N $1 -o $3.tmp $2 && mv $3.tmp $3 && samtools index -@ 12 $3; }
    calls() {
        local set=$1 T=$2 N=$3 t n th
        for t in severus_hp savana nanomonsv; do
            n=$t; [ $t = severus_hp ] && n=severus; th=16; [ $t = savana ] && th=8
            local D=$C/$n/$set; [ -s $D/calls.vcf ] && continue; mkdir -p $D
            /usr/bin/time -f "%e %M %x" -o $D/time.txt taskset -c ${CPUS:-0-15} bash $H/paper.sh caller $t $T $N $D/out $th > $D/stdout.txt 2>&1
            [ -s $D/out/calls.vcf ] && { cp $D/out/calls.vcf $D/out/log.txt $D/; rm -rf $D/out; }
        done
    }
    have() { local t; for t in severus savana nanomonsv; do [ -s $C/$t/$1/calls.vcf ] || return 1; done; }
    for p in ${*:-$CASTLE}; do
        W=$C/work/$p; B=$S/hifi/$p; mkdir -p $W/30x; ln -sfn $B/clair3 $W/30x/clair3
        for s in n t; do [ -s $C/names/$p.$s.tsv ] || samtools view -@ 8 -F 0x900 $B/$s.bam | awk -F'\t' '{print $1 "\t" length($10)}' > $C/names/$p.$s.tsv; done
        py select $p $C
        for x in 30 20 15 10 5; do
            [ -s $C/lists/$p.n.${x}x.txt ] || continue
            have ${p}_${x}x && continue
            [ -s $W/30x/t.bam.bai ] || sub $C/lists/$p.t.30x.txt $B/t.bam $W/30x/t.bam
            [ -s $W/30x/n.bam.bai ] || sub $C/lists/$p.n.30x.txt $B/n.bam $W/30x/n.bam
            if [ $x != 30 ]; then
                mkdir -p $W/${x}x; ln -sfn $B/clair3 $W/${x}x/clair3
                for s in n t; do [ -s $W/${x}x/$s.bam.bai ] || sub $C/lists/$p.$s.${x}x.txt $W/30x/$s.bam $W/${x}x/$s.bam; done
            fi
            calls ${p}_${x}x $W/${x}x/t.bam $W/${x}x/n.bam
            [ $x != 30 ] && rm -rf $W/${x}x
        done
        [ -s $C/lists/$p.n.purity.txt ] && mkdir -p $W/pn && ln -sfn $B/clair3 $W/pn/clair3
        [ -s $C/lists/$p.n.purity.txt ] && { [ -s $W/pn/n.bam.bai ] || sub $C/lists/$p.n.purity.txt $B/n.bam $W/pn/n.bam; }
        for q in 100 80 60 40 20; do
            [ -s $C/lists/$p.t.${q}p.tumor_part.txt ] || continue
            have ${p}_${q}p && continue
            [ $q = 100 ] && { [ -s $W/30x/t.bam.bai ] || sub $C/lists/$p.t.30x.txt $B/t.bam $W/30x/t.bam; }
            while [ $(jobs -rp | wc -l) -ge ${MIXJOBS:-1} ]; do wait -n; done
            (
            D=$W/p$q; mkdir -p $D
            if [ $q = 100 ]; then ln -sfn $W/30x/t.bam $D/t.bam; ln -sfn $W/30x/t.bam.bai $D/t.bam.bai
            else
                samtools merge -@ 8 -l 1 -f -o $D/t.bam <(samtools view -@ 4 -u -N $C/lists/$p.t.${q}p.tumor_part.txt $B/t.bam) <(samtools view -@ 4 -u -N $C/lists/$p.t.${q}p.normal_part.txt $B/n.bam) && samtools index -@ 8 $D/t.bam
            fi
            ln -sfn $W/pn/n.bam $D/n.bam; ln -sfn $W/pn/n.bam.bai $D/n.bam.bai; ln -sfn $B/clair3 $D/clair3
            calls ${p}_${q}p $D/t.bam $D/n.bam
            rm -rf $D
            ) &
        done
        wait
        rm -rf $W
    done ;;
consensus)
    E=${E:-$H/eval_sets}
    one() {
        local p=$1 W=$E/$1 hifi="" t s v; castle $p; mkdir -p $W/in; : > $W/in/SOURCES
        for t in $HIFI; do
            v=$R/competitors/$t/$p/calls.vcf; [ -s $v ] || v=$R/competitors/${t%_phased}/$p/calls.vcf; echo "$t $v" >> $W/in/SOURCES
            norm $v $t $W/in/$t.vcf; hifi="$hifi $W/in/$t.vcf"
        done
        [ -s $W/all_without/None_minda_ensemble.vcf ] || ens $W/all_without 10 $hifi $ONT $SR
        [ -n "$NOJANUS" ] && return
        grep -E '^#|PASS' ${VD:-$H/verdict}/$p.vcf > $W/in/snape.vcf
        ens $W/all_with 11 $hifi $W/in/snape.vcf $ONT $SR
        tset $W/all_without/None_minda_ensemble.vcf $W/in/snape.vcf $W/all_without_snape
        for t in sniffles2_mosaic delly; do
            [ -s $R/competitors/$t/$p/calls.vcf ] || continue
            norm $R/competitors/$t/$p/calls.vcf $t $W/in/$t.vcf
            for s in all_with all_without; do tset $W/$s/None_minda_ensemble.vcf $W/in/$t.vcf $W/${s}_$t; done
        done
    }
    for p in ${*:-$CASTLE}; do one $p & done; wait ;;
loo)
    E=${E:-$H/eval_sets}; LO=${LOO:-$H/eval_loo}
    one() {
        local p=$1 x t W hifi; castle $p
        for x in $HIFI snape; do
            W=$LO/$p/$x; hifi=""; mkdir -p $W; [ -s $W/score/results/overall.tsv ] && continue
            for t in $HIFI snape; do [ $t = $x ] || hifi="$hifi $E/$p/in/$t.vcf"; done
            [ -s $W/cons/None_minda_ensemble.vcf ] || ens $W/cons 10 $hifi $ONT $SR
            tset $W/cons/None_minda_ensemble.vcf $E/$p/in/$x.vcf $W/score
        done
    }
    for p in ${*:-$CASTLE}; do one $p & done; wait ;;
covscore)
    C=$R/coverage_callers; O=${O:-$R/coverage_eval}; mkdir -p $O
    base() { case $1 in hg008) echo $S/hifi/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz ;; colo829) echo $S/hifi/truth/colo829/truth.vcf ;; *) echo ${E:-$H/eval_sets}/$1/all_with/None_minda_ensemble.vcf ;; esac; }
    one() {
        local d=$O/$2/$3; [ -s $d/minda/results/overall.tsv ] && return; mkdir -p $d
        case $3 in snape) grep -E '^#|	PASS	' $4 > $d/in.vcf ;; *) cp $4 $d/in.vcf ;; esac
        tset $(base $1) $d/in.vcf $d/minda
    }
    for p in ${*:-$LABELED}; do
        for s in ${p}_30x ${p}_20x ${p}_15x ${p}_10x ${p}_5x ${p}_100p ${p}_80p ${p}_60p ${p}_40p ${p}_20p; do
            v=${VD:-$H/verdict}/coverage_matched/$s.vcf; [ ${s: -1} = p ] && v=${VD:-$H/verdict}/purity_matched/$s.vcf
            [ -s $v ] && one $p $s snape $v &
            for t in severus savana nanomonsv; do [ -s $C/$t/$s/calls.vcf ] && one $p $s $t $C/$t/$s/calls.vcf & done
        done
        wait
    done
    for f in $O/*/*/minda/results/overall.tsv; do
        s=$(basename $(dirname $(dirname $(dirname $(dirname $f))))); t=$(basename $(dirname $(dirname $(dirname $f))))
        awk -F'\t' -v s=$s -v t=$t 'NR==2{printf "%s\t%s\t%s\t%.4f\t%.4f\t%.4f\n", s, t, $8, $5, $6, $7}' $f
    done | sort > $O/scores.tsv ;;
labels)
    V11=${V11:-$H/supp/labels11}; O=$V11/minda; E=$H/eval_sets; mkdir -p $O
    ( cd $H && OUTV=$V11 OMP_NUM_THREADS=16 taskset -c 0-15 $V $H/paper.py labels > $V11/stdout.txt 2>&1 )
    row() { awk -F'\t' -v k="$2" 'NR>1 && (k == "" || $1 == k) {printf "%.4f\t%.4f\t%.4f\n", $5, $6, $7; exit}' $1; }
    pair() {
        local p=$1 W=$O/$1 hifi="" t x; castle $p; mkdir -p $W
        grep -E '^#|	PASS	' $V11/$p.vcf > $W/snape.vcf
        for t in $HIFI; do hifi="$hifi $E/$p/in/$t.vcf"; done
        ens $W/with 11 $hifi $W/snape.vcf $ONT $SR
        tset $E/$p/all_without/None_minda_ensemble.vcf $W/snape.vcf $W/without
        for x in 30 20 15 10 5; do grep -E '^#|	PASS	' $V11/coverage_matched/${p}_${x}x.vcf > $W/d$x.vcf; tset $E/$p/all_with/None_minda_ensemble.vcf $W/d$x.vcf $W/d$x; done
        echo -e "$p\twith\t$(awk -F'\t' '$1 ~ /^K_/ {printf "%.4f\t%.4f\t%.4f", $5, $6, $7}' $W/with/results/overall.tsv)"
        echo -e "$p\twithout\t$(row $W/without/results/overall.tsv)"
        for x in 30 20 15 10 5; do echo -e "$p\t${x}x\t$(row $W/d$x/results/overall.tsv)"; done
    }
    ext() {
        local W=$O/$1; mkdir -p $W; grep -E '^#|	PASS	' $V11/$1.vcf > $W/snape.vcf; tset $2 $W/snape.vcf $W/truth
        echo -e "$1\ttruth\t$(row $W/truth/results/overall.tsv)\t$(grep -vc '^#' $W/snape.vcf)"
    }
    for p in $CASTLE; do pair $p > $O/$p.tsv & done
    ext hg008 $S/hifi/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz > $O/hg008.tsv &
    ext colo829 $S/hifi/truth/colo829/truth.vcf > $O/colo829.tsv &
    wait; cat $O/*.tsv > $V11/scores.tsv; cat $V11/scores.tsv ;;
supp)
    for j in auc splits ablation features; do ( cd $H && $V $H/paper.py supp $j ) & done; wait ;;
grid)
    O=${SUPPOUT:-$H/supp}; mkdir -p $O
    for j in grid xgb; do for i in 0 1 2 3; do ( cd $H && SHARD=$i/4 OMP_NUM_THREADS=4 $V $H/paper.py supp $j > $O/$j.$i.log 2>&1 ) & done; done; wait
    for j in grid xgb; do awk 'FNR == 1 && NR != 1 {next} !s[$0]++' $O/$j.[0-9].tsv > $O/$j.tsv; done ;;
simmap)
    O=$R/sim; MMI=$S/simulation/phylo/ref/chm13.map-hifi.mmi; TUM=$(seq -f "T%02g" 1 20 | tr '\n' ' '); simreads
    [ -s $MMI ] || minimap2 -x map-hifi -d $MMI $S/ref/chm13v2.fa
    for x in 5 20; do
        d=$O/mm2/${x}_${x}; mkdir -p $d
        for s in germline $TUM; do
            [ -s $d/$s.time ] && continue; sync; sleep 5
            /usr/bin/time -f "%e %M %U %S %x" -o $d/$s.time.part taskset -c 0-15 bash -c "$O/svb2fq $(readlink -f $(reads $x $s.svb)) $s 2> /dev/null | minimap2 -ax map-hifi -t 16 -Y --secondary=no $MMI - 2> $d/$s.minimap2.log | samtools sort -@ 4 -m 2G -T $d/$s.st -o $d/$s.bam - && samtools index -@ 8 $d/$s.bam"
            rm -f $d/$s.bam $d/$s.bam.bai; mv $d/$s.time.part $d/$s.time
        done
    done ;;
sweep)
    W=$R/sweep; mkdir -p $W/runs
    one() {
        local k=$1 n=$2 a=$3 b=$4 p=${2%%_*} D=$W/runs/$1/$2/t$3_$4 MX=$W/mix/${2%%_*} x
        mkdir -p $D/work; log "[$k/$n theta $a/$b] $([ -e $D/work/done.pack_t ] && echo resume || echo start)"
        [ -e $D/work/done.pack_t ] || case $k in
            coverage) x=${n##*_}; x=${x%x}; for s in n t; do py mix $D/work/$s.svb $MX/$s.svb:0:$((x * G)) >> $Q 2>&1; done; touch $D/work/done.pack_n $D/work/done.pack_t ;;
            purity) x=${n##*_}; x=${x%p}; for f in svb svb.idx; do ln -sf $MX/nm.$f $D/work/n.$f; ln -sf $MX/t$x.$f $D/work/t.$f; done; touch $D/work/done.pack_n $D/work/done.pack_t ;;
        esac
        SHALLOW=$a CLONAL=$b bash $S/pipeline.sh $D/work $(src $p n) $(src $p t) 16 > $D/pipeline.log 2>&1
        if [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ]; then
            mkdir -p $D/logs; mv $D/work/calls.vcf $D/work/n.cnt $D/work/t.cnt $D/work/t.cnt.hist $D/; mv $D/work/log.* $D/work/steps.txt $D/logs/
            rm -rf $D/work; log "$(py sweep check $k $n $D)"
        else log "[$k/$n theta $a/$b] FAILED: $(tail -1 $D/pipeline.log) $(tail -qn 1 $D/work/log.* 2>/dev/null | grep -m 1 -i 'cannot\|error\|exit')"; [ -e $D/attempt ] && touch $D/failed || touch $D/attempt; fi
    }
    room() { df --output=avail -BG /home | tail -1 | tr -dc 0-9; }
    for p in ${*:-$CASTLE}; do
        MX=$W/mix/$p
        while :; do
            L=$(py sweep plan $p | while read k n a b m0 mt; do [ -e $W/runs/$k/$n/t${a}_$b/calls.vcf ] || [ -e $W/runs/$k/$n/t${a}_$b/failed ] || echo "$k $n $a $b"; done)
            [ -z "$L" ] && break
            NEW=$(echo "$L" | while read k n a b; do [ -e $W/runs/$k/$n/t${a}_$b/work/done.pack_t ] || echo "$k $n $a $b"; done)
            if echo "$NEW" | grep -q '^coverage\|^purity'; then
                mkdir -p $MX
                for s in n t; do [ -s $MX/$s.svb.idx ] || $S/build/janus pack -i $(src $p $s) -o $MX/$s.svb -t 16 -s 1 > $MX/pack_$s.log 2>&1; done
            fi
            if echo "$NEW" | grep -q '^purity'; then
                nb=$(python3 -c "import struct; print(struct.unpack_from('<Q', open('$MX/n.svb', 'rb').read(32), 16)[0])")
                nm=$((nb - 24 * G - G / 10)); [ $nm -gt $((30 * G)) ] && nm=$((30 * G))
                [ -s $MX/nm.svb.idx ] || py mix $MX/nm.svb $MX/n.svb:0:$nm >> $Q 2>&1
                for q in $(echo "$NEW" | awk '$1 == "purity" {sub(/.*_/, "", $2); sub(/p$/, "", $2); print $2}' | sort -u); do
                    [ -s $MX/t$q.svb.idx ] || py mix $MX/t$q.svb $MX/t.svb:0:$((q * 30 * G / 100)) $([ $q -lt 100 ] && echo $MX/n.svb:$nm:$(((100 - q) * 30 * G / 100))) >> $Q 2>&1
                done
            fi
            echo "$L" | { while read k n a b; do while [ $(jobs -rp | wc -l) -ge 2 ] || { [ $(jobs -rp | wc -l) -gt 0 ] && [ $(room) -lt 150 ]; }; do wait -n; done; one $k $n $a $b & done; wait; }
        done
        rm -rf $MX
    done
    ( cd $H && py sweep configs ) > $W/configs.txt 2> $W/configs.err
    for k in $(cat $W/configs.txt); do
        C=$W/conf/$k; mkdir -p $C/verdict
        [ -s $C/verdict/folds.json ] || ( cd $H && RUNS=$C/full COV=$C/coverage PUR=$C/purity OUTV=$C/verdict OMP_NUM_THREADS=16 $V $H/paper.py verdict > $C/verdict/stdout.txt 2>&1 )
        log "[config $k] $(tail -1 $C/verdict/log.txt)"
        sc() { local d=$C/score/$1/$2; [ -s $d/minda/results/overall.tsv ] && return; mkdir -p $d; grep -E '^#|	PASS	' $3 > $d/in.vcf; tset $H/eval_sets/${2%%_*}/all_$1/None_minda_ensemble.vcf $d/in.vcf $d/minda; }
        for p in $CASTLE; do
            sc without $p $C/verdict/$p.vcf & sc with $p $C/verdict/$p.vcf &
            for f in $C/verdict/coverage_matched/${p}_*x.vcf $C/verdict/purity_matched/${p}_*p.vcf; do [ -s $f ] && sc without $(basename $f .vcf) $f & done
            wait
        done
    done
    py sweep table > $W/summary.txt; log "sweep done" ;;
sweepscore)
    W=$R/sweep
    ev() {
        local k=$1 C=$W/conf/$1; log "[config $k] start"
        E=$C/eval_sets VD=$C/verdict bash $H/paper.sh consensus
        E=$C/eval_sets LOO=$C/eval_loo bash $H/paper.sh loo
        E=$C/eval_sets VD=$C/verdict O=$C/coverage_eval bash $H/paper.sh covscore $CASTLE
        py sweep rescore $C > $C/rescore.tsv; log "[config $k] done"
    }
    while :; do
        a=$(systemctl --user is-active janus-sweep)
        for k in $(printf '%s\n' $* $(cat $W/configs.txt) | awk '!s[$0]++'); do
            [ -s $W/conf/$k/rescore.tsv ] || [ ! -s $W/conf/$k/verdict/folds.json ] && continue
            while [ $(jobs -rp | wc -l) -ge 2 ]; do wait -n; done; ev $k &
        done; wait
        [ "$a" = active ] || break; sleep 600
    done
    py sweep rescores > $W/rescore.txt ;;
timing)
    RC=$H/runs_clean; Q=$RC/queue.log; mkdir -p $RC; B=$S/build/janus; TUM=$(seq -f "T%02g" 1 20 | tr '\n' ' ')
    for p in $CASTLE; do janus_run full $RC/full/$p $p ""; done
    O=$R/sim; RF=$S/ref/chm13v2; simreads
    tm() { sync; sleep 10; /usr/bin/time -f "%e %M %U %S %x" -o $1 taskset -c 0-15 "${@:2}"; }
    for set in 20_20 20_10 20_5 10_10 5_5; do
        nx=${set%_*}; tx=${set#*_}
        for t in $TUM; do
            p=$RC/sim/timing/$set/pair/$t; [ -s $p/time.txt ] && continue; mkdir -p $p/work
            for e in svb svb.idx; do ln -sf $(readlink -f $(reads $nx germline.$e)) $p/work/n.$e; ln -sf $(readlink -f $(reads $tx $t.$e)) $p/work/t.$e; done
            touch $p/work/done.pack_n $p/work/done.pack_t
            REF=chm13v2 tm $p/time.part bash $S/pipeline.sh $p/work $p/work/n.svb $p/work/t.svb 16 > $p/pipeline.log 2>&1
            cp $p/work/steps.txt $p/work/calls.vcf $p/; [ -s $p/calls.vcf ] && mv $p/time.part $p/time.txt && rm -rf $p/work
        done
        d=$RC/sim/timing/$set/joint; [ -s $d/time.txt ] || { tm $d.time bash $H/paper.sh joint $nx $tx $d; mv $d.time $d/time.txt; }
    done ;;
figures)
    for i in latex-tikz phangorn; do n=$i; [ $i = phangorn ] && n=janus-bench/phangorn; docker image inspect $n > /dev/null 2>&1 || docker build -q -t $n $K/images/$i > /dev/null; done
    tex() { docker run --rm --log-driver none -v $1:/w -w /w -u "$(id -u):$(id -g)" -e HOME=/tmp latex-tikz sh -c "$2"; }
    F=$S/paper/overleaf; mkdir -p $F/figures $H/out
    tex $S/paper/figures/system "pdflatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp 20_system.tex > /dev/null && cp /tmp/20_system.pdf /w/" && cp $S/paper/figures/system/20_system.pdf $F/figures/fig_system.pdf
    $V $H/figures.py > $H/out/figures.log
    tex $F "pdflatex -interaction=nonstopmode main.tex > /dev/null; bibtex main > /dev/null; pdflatex -interaction=nonstopmode main.tex > /dev/null; pdflatex -interaction=nonstopmode main.tex | grep -E '^!|undefined'" ;;
joint)
    joint "$@" ;;
guard)
    case $1 in
    memory) while :; do a=$(awk '/MemAvailable/{print int($2 / 1048576)}' /proc/meminfo)
            [ "$a" -lt 14 ] && { echo "$(date '+%F %T') memguard: $a GB available; killing containers" >> $R/night/INCIDENTS.md; c=$(docker ps -q); [ -n "$c" ] && docker kill $c > /dev/null; sleep 90; }
            sleep 4; done ;;
    disk) held=0; while :; do a=$(df --output=avail -BG /home | tail -1 | tr -dc 0-9)
            if [ $held = 0 ] && [ $a -lt 60 ]; then c=$(docker ps -q); [ -n "$c" ] && docker pause $c > /dev/null; echo "$(date '+%F %T') DISK GUARD: $a GB free; containers paused" >> $R/night/INCIDENTS.md; held=1
            elif [ $held = 1 ] && [ $a -ge 120 ]; then [ -n "$c" ] && docker unpause $c > /dev/null; echo "$(date '+%F %T') DISK GUARD: $a GB free, resumed" >> $R/night/INCIDENTS.md; held=0; fi
            sleep 60; done ;;
    load) O=$R/load.tsv; [ -s $O ] || echo -e "time\tload1\tall_cores" > $O
          cpu() { awk '/^cpu /{print $2 + $3 + $4 + $7 + $8}' /proc/stat; }; a=$(cpu); t=$(date +%s.%N); h=$(getconf CLK_TCK)
          while :; do sleep 30; b=$(cpu); u=$(date +%s.%N)
              awk -v a=$a -v b=$b -v t=$t -v u=$u -v h=$h -v l="$(cut -d' ' -f1 /proc/loadavg)" -v d="$(date '+%F %T')" 'BEGIN{printf "%s\t%s\t%.1f\n", d, l, (b - a) / h / (u - t)}' >> $O
              a=$b; t=$u; done ;;
    esac ;;
caller)
    caller "$@" ;;
*)
    echo "unknown experiment $X; see the header of $0" >&2; exit 1 ;;
esac
