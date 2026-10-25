#!/usr/bin/env bash
# The references of janus, built from public sources into ref/:
#   hg38_parmask.fa    UCSC hg38 without the alternate haplotypes (194 sequences), upper case, the chrY pseudoautosomal
#                      regions (PAR1 10001-2781479, PAR2 56887903-57217415) hard-masked so that PAR reads map to chrX
#   hg38_parmask.trf.bed   the tandem repeats of GRCh38 (the VNTR annotation shipped with Severus)
#   chm13v2.fa         T2T-CHM13v2.0 (the reference of the simulation)
#   chm13v2.trf.bed    its simple repeats and satellites (RepeatMasker annotation of the T2T consortium)
#   <ref>.cidx         the LCP core index of janus (janus refidx) and <ref>.sr.mmi, the minimap2 short-read index, of each
#   ref/build.sh [hg38|chm13] (both by default); needs curl, samtools, minimap2 and build/janus
set -euo pipefail
R=$(cd "$(dirname "$0")" && pwd); J=$R/../build/janus; cd $R
T2T=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies
index() { [ -s $1.cidx ] || $J refidx --ref $1.fa --out $1.cidx --threads 16; [ -s $1.sr.mmi ] || minimap2 -x sr -d $1.sr.mmi $1.fa; }
for ref in ${*:-hg38 chm13}; do
case $ref in
hg38)
    if [ ! -s hg38_parmask.fa.fai ]; then
        [ -s hg38.fa.gz ] || curl -sSfL -o hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        python3 - <<'EOF'
import gzip, re
PAR = [(10001, 2781479), (56887903, 57217415)]
MAIN = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]]
seqs = {}; name = None; seq = []
def flush():
    if name is None or name.endswith("_alt"): return
    s = "".join(seq).upper()
    if name == "chrY":
        s = list(s)
        for a, b in PAR: s[a - 1:b] = "N" * (b - a + 1)
        s = "".join(s)
    seqs[name] = s
for l in gzip.open("hg38.fa.gz", "rt"):
    if l[0] == ">": flush(); name = l[1:].split()[0]; seq = []
    else: seq.append(l.strip())
flush()
# the primary chromosomes in order, then the unplaced and unlocalized contigs by their GenBank accession (chrUn_GL000195v1 -> GL000195.1)
acc = {n: "{}.{}".format(*re.search(r"((?:GL|KI)\d+)v(\d+)", n).groups()) for n in seqs if n not in MAIN}
out = open("hg38_parmask.fa.tmp", "w")
for n in MAIN + sorted(acc, key=acc.get):
    s = seqs[n]; head = f"{n} {'MT' if n == 'chrM' else n[3:]}" if n in MAIN else f"{acc[n]} {acc[n]}"
    out.write(">" + head + "\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
out.close()
EOF
        mv hg38_parmask.fa.tmp hg38_parmask.fa && samtools faidx hg38_parmask.fa
    fi
    [ -s hg38_parmask.trf.bed ] || curl -sSfL -o hg38_parmask.trf.bed https://raw.githubusercontent.com/KolmogorovLab/Severus/main/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
    index hg38_parmask ;;
chm13)
    [ -s chm13v2.fa.fai ] || { curl -sSfL -o chm13v2.fa.gz $T2T/analysis_set/chm13v2.0.fa.gz && gzip -d chm13v2.fa.gz && samtools faidx chm13v2.fa; }
    [ -s chm13v2.trf.bed ] || curl -sSfL $T2T/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed | awk -F'\t' '$7 == "Simple_repeat" || $7 == "Satellite" {print $1 "\t" $2 "\t" $3}' > chm13v2.trf.bed
    index chm13v2 ;;
*) echo "unknown reference $ref" >&2; exit 1 ;;
esac
done
echo "references ready: $(ls *.cidx | tr '\n' ' ')"
