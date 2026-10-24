#!/usr/bin/env bash
# The simulated cohort of the paper: a diploid normal genome carrying HG00733's germline variants (HGSVC3, T2T-CHM13),
# a tree of 20 tumor clones with somatic SNVs, indels, SVs and copy-number changes (spec/), and Badread HiFi reads at
# 20x of the normal and of every tumor. The steps run in order and resume where they stopped; every number of the
# simulation comes from spec/somatic.yaml and spec/parameters.yaml, the seed is 20260818 throughout.
#   datasets/simulate.sh env        the python environment: insilicoSV 0.1.9 with the patch in patches_insilicosv_v0.1.9.patch, Badread 0.4.2, and the COSMIC v3.4 signatures
#   datasets/simulate.sh ref        ref/chm13.primary.fa (chr1-22 and chrX of ../../ref/chm13v2.fa) and the CHM13 centromere annotation
#   datasets/simulate.sh germline   genomes/germline: HG00733's phased SNVs, indels and SVs (datasets/hifi/truth/hgsvc3) applied per chromosome with insilicoSV, verified base by base
#   datasets/simulate.sh tree       genomes/tumor: the 39-node tree of spec/tree.json realized level by level (scripts/realize_tree.py)
#   datasets/simulate.sh cnv        the copy-number profile of the tree and the leaf genomes with it applied (genomes/tumor/<leaf>/genome_cnv.fa)
#   datasets/simulate.sh reads      Badread reads of the normal and the 20 leaves at 20x, merged into reads/<sample>/<sample>.fastq.gz and packed for janus (<sample>.svb)
#   datasets/simulate.sh all
set -uo pipefail
PH=$(cd "$(dirname "$0")/simulation" && pwd); S=$(cd "$PH/../.." && pwd); V=$PH/.venv; PY=$V/bin/python; cd "$PH"
DEPTH=20; NP_GERMLINE=11; NP_TREE=3; NT_READS=30
say() { echo "[$(date '+%F %T')] $*" | tee -a logs/simulate.log; }
mkdir -p logs
step=${1:-all}

if [ $step = env ] || [ $step = all ]; then
    if [ ! -x $V/bin/insilicosv ]; then
        uv venv --python 3.12 $V && uv pip install --python $PY numpy pyyaml scipy git+https://github.com/rrwick/Badread.git@v0.4.2 git+https://github.com/PopicLab/insilicoSV.git@v0.1.9
        ( cd $($PY -c "import insilicosv, os; print(os.path.dirname(os.path.dirname(insilicosv.__file__)))") && patch -p1 < $PH/patches_insilicosv_v0.1.9.patch )
        bash scripts/check_patches.sh || exit 1
    fi
    mkdir -p spec/cosmic
    for f in COSMIC_v3.4_SBS_GRCh38.txt COSMIC_v3.4_ID_GRCh37.txt; do [ -s spec/cosmic/$f ] || curl -sSfL -o spec/cosmic/$f https://cog.sanger.ac.uk/cosmic-signatures-production/documents/$f; done
    say "environment ready"
fi
if [ $step = ref ] || [ $step = all ]; then
    mkdir -p ref
    [ -s ref/chm13.primary.fa.fai ] || { samtools faidx $S/ref/chm13v2.fa $(seq -f chr%g 1 22) chrX > ref/chm13.primary.fa.tmp && mv ref/chm13.primary.fa.tmp ref/chm13.primary.fa && samtools faidx ref/chm13.primary.fa; }
    [ -s ref/chm13_censat_v2.1.bed ] || curl -sSfL -o ref/chm13_censat_v2.1.bed https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed
    $PY scripts/refinfo.py | head -3; say "reference ready"
fi
if [ $step = germline ] || [ $step = all ]; then
    if [ ! -s genomes/germline/germline.fa.fai ]; then
        bash scripts/build_germline.sh HG00733 $S/datasets/hifi/truth/hgsvc3 genomes/germline $NP_GERMLINE >> logs/build_germline.log 2>&1
        n=$(grep -h 'TAMAM' genomes/germline/logs/*.log 2>/dev/null | wc -l); [ "$n" -eq 23 ] || { say "only $n of 23 chromosomes finished; see genomes/germline/logs"; exit 1; }
        bash scripts/assemble_germline.sh genomes/germline >> logs/assemble_germline.log 2>&1 || exit 1
    fi
    say "germline genome: $(awk '{s+=$2} END{printf "%d contigs, %.3f Gb", NR, s/1e9}' genomes/germline/germline.fa.fai)"
fi
if [ $step = tree ] || [ $step = all ]; then
    NP_POST=$NP_TREE bash scripts/run_tree_parallel.sh genomes/germline/germline.fa genomes/tumor $NP_TREE >> logs/tree_build.log 2>&1 || exit 1
    say "tree: $(ls -d genomes/tumor/*/ | wc -l) node directories"
fi
if [ $step = cnv ] || [ $step = all ]; then
    [ -s genomes/tumor/cn_profile.json ] || $PY scripts/build_cn_profile.py --fai genomes/tumor/MRCA/genome.fa.fai --tree genomes/tumor/plan.json --spec spec/somatic.yaml --out genomes/tumor/cn_profile.json >> logs/cn_profile.log 2>&1
    for t in $(ls genomes/tumor | grep -E '^T[0-9]+$'); do [ -s genomes/tumor/$t/cnv_report.json ] || echo $t; done | xargs -P 4 -I{} $PY scripts/apply_cnv.py --node {} --tumor-dir genomes/tumor --profile genomes/tumor/cn_profile.json >> logs/apply_cnv.log 2>&1
    say "copy number applied to $(ls genomes/tumor/T*/cnv_report.json | wc -l) leaves"
fi
if [ $step = reads ] || [ $step = all ]; then
    mkdir -p reads
    md5sum genomes/germline/germline.fa genomes/tumor/T*/genome_cnv.fa > reads/genome_checksums.md5
    for s in germline $(ls genomes/tumor | grep -E '^T[0-9]+$'); do
        g=genomes/germline/germline.fa; [ $s = germline ] || g=genomes/tumor/$s/genome_cnv.fa
        [ -s reads/$s/summary.txt ] || bash scripts/sim_reads.sh $g reads/$s $DEPTH $NT_READS >> logs/reads_$s.log 2>&1 || { say "reads of $s failed"; exit 1; }
        [ -s reads/$s/$s.fastq.gz ] || bash scripts/merge_shards.sh $s >> logs/reads_$s.log 2>&1
        [ -s reads/$s/$s.svb.idx ] || $S/build/janus pack --reads reads/$s/$s.fastq.gz --out reads/$s/$s.svb --threads 16 > reads/$s/pack.log 2>&1
        say "$s: $(grep 'toplam baz' reads/$s/summary.txt)"
    done
    bash scripts/bgzip_genomes.sh >> logs/bgzip_genomes.log 2>&1
    say "reads ready: $(ls reads/*/*.svb | wc -l) samples"
fi
