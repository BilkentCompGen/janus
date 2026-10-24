#!/usr/bin/env bash
# insilicoSV yamalarinin YERINDE oldugunu dogrular. Faz C/D/E'nin her kosusundan ONCE calistir.
#
# NEDEN: yamalar calisma agacinda COMMIT EDILMEMIS durumda. Tek bir `pip install -e .`,
# `git checkout` veya `git stash` ikisini de geri alir ve:
#   - yama 1 kaybolursa: import edilen TUM insersiyonlar SESSIZCE 1 bp'ye kirpilir
#   - yama 2 kaybolursa: output_adjacencies her SV kosusunda cokertir
# Ilki UYARI VERMEZ; ancak en sondaki spektrum dogrulamasinda fark edilirdi.
set -euo pipefail
PHYLO=$(cd "$(dirname "$0")/.." && pwd); F=$($PHYLO/.venv/bin/python -c "import insilicosv, os; print(os.path.join(os.path.dirname(insilicosv.__file__), 'variant_set.py'))")
n=$(grep -c "PATCH (snape/phylo 2026-08-18)" "$F" || true)
if [ "$n" -lt 2 ]; then
  echo "HATA: insilicoSV yamalari EKSIK (bulunan: $n/2)." >&2
  echo "apply: cd \$(dirname \$(dirname $F)) && patch -p1 < $PHYLO/patches_insilicosv_v0.1.9.patch" >&2
  exit 1
fi
# ampirik teyit: 5 bp'lik bir INS gercekten 5 bp mi ekliyor
W=$(mktemp -d)
python3 - "$W" <<'PY'
import sys, random
w=sys.argv[1]; random.seed(1)
s="".join(random.choice("ACGT") for _ in range(2000))
open(f"{w}/ref.fa","w").write(">c\n"+"\n".join(s[i:i+60] for i in range(0,2000,60))+"\n")
open(f"{w}/in.vcf","w").write(
 "##fileformat=VCFv4.2\n##contig=<ID=c,length=2000>\n"
 '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
 '##INFO=<ID=END,Number=1,Type=Integer,Description="e">\n'
 '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="l">\n'
 '##INFO=<ID=INSSEQ,Number=1,Type=String,Description="s">\n'
 '##ALT=<ID=INS,Description="i">\n'
 '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
 "c\t1000\ti0\tN\t<INS>\t100\tPASS\tSVTYPE=INS;END=1000;SVLEN=5;INSSEQ=ACGTA\tGT\t1|1\n")
open(f"{w}/config.yaml","w").write(f'reference: "{w}/ref.fa"\nrandom_seed: 1\nvariant_sets:\n    - import: "{w}/in.vcf"\n')
PY
( cd "$W" && $PHYLO/.venv/bin/insilicosv -c config.yaml >/dev/null 2>&1 )
L=$(grep -v '^>' "$W/sim.hapA.fa" | tr -d '\n' | wc -c)
if [ "$L" -ne 2005 ]; then
  echo "HATA: 5 bp INS testi basarisiz (uzunluk $L, beklenen 2005). Yama 1 etkin degil." >&2
  exit 1
fi
echo "insilicoSV yamalari YERINDE ve etkin (5 bp INS -> +5 bp dogrulandi)"
