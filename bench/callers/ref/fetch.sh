#!/usr/bin/env bash
# fetch.sh: extra references for the competitor wrappers: T2T-CHM13v2.0 (UCSC hs1) and the CHM13 minigraph pangenome
# (HPRC 464 haplotypes, zenodo 16728828) for minisv; delly's hg38 exclude list; the SVision-pro v2.5 repo (model + access bed)
cd /home/yuio/snape/bench/callers/ref
[ -s hs1.done ] || { curl -sL -o hs1.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz && echo ok > hs1.done; }
[ -s gfa.done ] || { curl -sL -o CHM13-464.gfa.gz https://zenodo.org/api/records/16728828/files/CHM13-464.gfa.gz/content && echo ok > gfa.done; }
[ -s human.hg38.excl.tsv ] || curl -sL -o human.hg38.excl.tsv https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg38.excl.tsv
[ -d svision-pro ] || { git clone -q https://github.com/songbowang125/SVision-pro.git svision-pro && git -C svision-pro checkout -q v2.5; }
