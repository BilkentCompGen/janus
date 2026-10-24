#!/bin/bash
# The environment of the benchmark: the python environment of the scripts, Minda, the other callers (pip installs in
# bench/callers/env, container images for the rest), their references and the docker images of the figures.
#   bench/setup.sh
# Needs: uv, docker, minimap2, samtools, bcftools, tabix/bgzip, R (ape, phangorn via docker), aria2c, pigz, sensors/pidstat/iostat.
source "$(dirname "$0")/env.sh"
set -e
cd $S
# python of the benchmark scripts (scikit-learn, matplotlib, xgboost for the figures and the classifier experiments)
[ -x $PY ] || { uv venv --python 3.14 .venv-bench && uv pip install --python $PY scikit-learn==1.9.1 numpy==2.5.3 scipy==1.18.1 matplotlib==3.11.2 xgboost==3.4.1; }
# Minda (consensus and scoring)
[ -x $MINDA ] || uv tool install --from git+https://github.com/KolmogorovLab/minda minda --with pandas --with numpy --with pysam --with pybedtools
# the callers that run outside containers, each in its own environment
mkdir -p $L
[ -x $L/severus/bin/severus ] || { uv venv --python 3.12 $L/severus && uv pip install --python $L/severus/bin/python severus==1.7; }
[ -x $L/savana/bin/savana ] || { uv venv --python 3.11 $L/savana && uv pip install --python $L/savana/bin/python savana==1.3.8; }
[ -x $L/sniffles/bin/sniffles ] || { uv venv --python 3.12 $L/sniffles && uv pip install --python $L/sniffles/bin/python sniffles==2.8.0; }
[ -x $L/nanomonsv/bin/nanomonsv ] || { uv venv --python 3.12 $L/nanomonsv && uv pip install --python $L/nanomonsv/bin/python nanomonsv==0.9.0; }
# container images: the callers come from biocontainers at the versions named in env.sh; these three are built here
for i in latex-tikz phangorn svision-pro; do n=janus-bench/$i; [ $i = latex-tikz ] && n=latex-tikz; [ $i = svision-pro ] && n=janus-bench/svision-pro:2.5-pysam0.21; docker image inspect $n > /dev/null 2>&1 || docker build -q -t $n $K/images/$i; done
# the references of the callers: the minimap2 index for the alignments, SVDSS's FMD index, Delly's exclusion list, SVision-pro's model
mkdir -p $K/ref
[ -s $K/ref/hg38_parmask.mmi ] || minimap2 -x map-hifi -d $K/ref/hg38_parmask.mmi $REF/hg38_parmask.fa
[ -s $K/ref/human.hg38.excl.tsv ] || curl -sL -o $K/ref/human.hg38.excl.tsv https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv
[ -d $K/ref/svision-pro ] || { git clone -q https://github.com/songbowang125/SVision-pro.git $K/ref/svision-pro && git -C $K/ref/svision-pro checkout -q v2.5; }
[ -s $K/ref/svdss/index.done ] || { mkdir -p $K/ref/svdss; docker run --rm --log-driver none -u "$(id -u):$(id -g)" -v $S:$S quay.io/biocontainers/svdss:2.1.1--he17396a_0 SVDSS index --reference $REF/hg38_parmask.fa --index $K/ref/svdss/hg38_parmask.fmd --threads $TH && touch $K/ref/svdss/index.done; }
echo "benchmark environment ready"
