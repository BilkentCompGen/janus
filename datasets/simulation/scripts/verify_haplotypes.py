#!/usr/bin/env python3
"""Uretilen haplotipleri BAGIMSIZ olarak yeniden insa edip karsilastirir.

    verify_haplotypes.py --ref ref.fa --hapA sim.hapA.fa --hapB sim.hapB.fa \
                         --import-vcf import.vcf --chrom chr20

FIKIR
-----
insilicoSV'nin kendi raporuna (sim.vcf, sim.stats.txt) guvenmek dairesel bir dogrulamadir:
uretici "yaptim" der, biz de ona inaniriz. Bunun yerine burada **varyant listesinden beklenen
haplotip dizisini sifirdan kendimiz kuruyoruz** ve insilicoSV'nin urettigiyle bayt bayt
karsilastiriyoruz. Esitlik, su dortunun AYNI ANDA dogru oldugunun kanitidir:

    konumlar · koordinat sozlesmesi (off-by-one) · faz (hangi haplotip) · zigosite

Esit degilse ilk farkli konum ve cevresi basilir; ayrica INS/DEL icin alternatif koordinat
sozlesmeleri denenip hangisinin tuttugu RAPORLANIR (varsayim yerine olcum).
"""
import argparse
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def parse_import_vcf(path):
    """-> [(pos, end, svtype, insseq, hapA, hapB)] 1-tabanli, insilicoSV sozlesmesi"""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos = int(f[1])
            ref, alt = f[3], f[4]
            info = {}
            if f[7] != ".":
                for kv in f[7].split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        info[k] = v
            gt = f[9].split(":")[0].replace("/", "|").split("|")
            hA, hB = gt[0] != "0", gt[1] != "0"
            svt = info.get("SVTYPE")
            if svt is None:
                svt = "SNP"
                end = pos
            else:
                end = int(info.get("END", pos))
            out.append({"pos": pos, "end": end, "type": svt, "ref": ref, "alt": alt,
                        "insseq": info.get("INSSEQ"), "hapA": hA, "hapB": hB})
    return out


def build_hap(ref, variants, which, ins_offset=0):
    """Beklenen haplotipi kur. which: 'hapA' | 'hapB'.
    ins_offset: INS icin POS yorumu denemesi (0 => POS bazinin ONUNE ekle)."""
    vs = [v for v in variants if v[which]]
    vs.sort(key=lambda v: (v["pos"], v["end"]))
    out = []
    cur = 0                       # 0-tabanli, ref uzerinde islenmis konum
    for v in vs:
        p, e, t = v["pos"], v["end"], v["type"]
        if t == "SNP":
            s = p - 1
            if s < cur:
                return None, f"ortusme SNP@{p} (cur={cur})"
            out.append(ref[cur:s])
            out.append(v["alt"])
            cur = s + 1
        elif t == "DEL":
            s = p - 1
            if s < cur:
                return None, f"ortusme DEL@{p} (cur={cur})"
            out.append(ref[cur:s])
            cur = e                       # e = silinen son baz (1-tabanli) -> 0-tabanli e
        elif t == "INS":
            s = p - 1 + ins_offset
            if s < cur:
                return None, f"ortusme INS@{p} (cur={cur})"
            out.append(ref[cur:s])
            out.append(v["insseq"] or "")
            cur = s
        elif t == "DUP":
            s = p - 1
            if s < cur:
                return None, f"ortusme DUP@{p} (cur={cur})"
            out.append(ref[cur:e])        # ref'i e'ye kadar yaz
            out.append(ref[s:e])          # segmenti tekrar yaz (tandem)
            cur = e
        elif t == "INV":
            s = p - 1
            if s < cur:
                return None, f"ortusme INV@{p} (cur={cur})"
            out.append(ref[cur:s])
            out.append(ref[s:e][::-1].translate(COMP))
            cur = e
        else:
            return None, f"bilinmeyen tip {t}@{p}"
    out.append(ref[cur:])
    return "".join(out), None


def first_diff(a, b):
    n = min(len(a), len(b))
    lo, hi = 0, n
    if a[:n] == b[:n]:
        return n if len(a) != len(b) else -1
    while lo < hi:
        mid = (lo + hi) // 2
        if a[:mid] == b[:mid]:
            lo = mid + 1
        else:
            hi = mid
    return lo - 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--hapA", required=True)
    ap.add_argument("--hapB", required=True)
    ap.add_argument("--import-vcf", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--sample-n", type=int, default=0)
    a = ap.parse_args()

    ref = list(read_fasta(a.ref).values())[0].upper()
    A = list(read_fasta(a.hapA).values())[0].upper()
    B = list(read_fasta(a.hapB).values())[0].upper()
    variants = parse_import_vcf(a.import_vcf)
    print(f"ref={len(ref):,}  hapA={len(A):,}  hapB={len(B):,}  variants={len(variants):,}")
    nA = sum(1 for v in variants if v["hapA"])
    nB = sum(1 for v in variants if v["hapB"])
    nHom = sum(1 for v in variants if v["hapA"] and v["hapB"])
    print(f"beklenen: hapA'da {nA:,}  hapB'de {nB:,}  homozigot {nHom:,}  "
          f"het/hom={(nA+nB-2*nHom)/nHom if nHom else float('nan'):.3f}")

    ok = True
    for which, obs in (("hapA", A), ("hapB", B)):
        best = None
        # ONCE ins_offset=0 denenir ve TUTARSA digerleri atlanir. Sozlesme
        # scripts/test_ins_convention.sh ile ampirik olarak sabitlendi (offset 0); kalan iki
        # deneme yalnizca bir uyusmazlikta TANILAMA icindir. Kromozom olceginde her deneme
        # 250 MB'lik bir dize kurdugu icin bu kisayol ucte bir maliyet demek.
        for off in (0, 1, -1):
            exp, err = build_hap(ref, variants, which, ins_offset=off)
            if exp is None:
                print(f"  {which}: ins_offset={off:+d} -> kurulamadi ({err})")
                continue
            same = exp == obs
            d = -1 if same else first_diff(exp, obs)
            print(f"  {which}: ins_offset={off:+d}  beklenen_uzunluk={len(exp):,} "
                  f"gozlenen={len(obs):,}  {'ESIT' if same else f'ilk fark @{d:,}'}")
            if same and best is None:
                best = off
                break                      # sozlesme tuttu; kalan denemelere gerek yok
        if best is None:
            ok = False
            # tanilama: hicbir sozlesme tutmadi -> ilk farkin cevresini bas
            exp, _ = build_hap(ref, variants, which, ins_offset=0)
            if exp:
                d = first_diff(exp, obs)
                lo = max(0, d - 40)
                print(f"    !! {which} hicbir sozlesmede esitlenmedi. ilk fark @{d:,}")
                print(f"    beklenen : ...{exp[lo:d]}[{exp[d:d+20]}]...")
                print(f"    gozlenen : ...{obs[lo:d]}[{obs[d:d+20]}]...")
        else:
            print(f"    OK {which}: ins_offset={best:+d} sozlesmesiyle BAYT BAYT ESIT")

    print("\nSONUC:", "GECTI — uretilen haplotipler bagimsiz yeniden insa ile birebir ayni"
          if ok else "KALDI — uretilen haplotipler beklenenle uyusmuyor")
    sys.exit(0 if ok else 2)


if __name__ == "__main__":
    main()
