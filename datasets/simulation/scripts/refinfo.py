#!/usr/bin/env python3
"""Referansa bagli TUM sabitlerin TEK kaynagi. Hicbir yerde kromozom uzunlugu GOMULMEZ.

    import refinfo
    L    = refinfo.lengths()          # {"chr1": 248387328, ...}
    CEN  = refinfo.centromeres()      # {"chr1": 124100716, ...}  (bp, HOR dizisi ortasi)
    HAPB = refinfo.hap_ref_bp()       # chr1-22 + chrX toplami

NEDEN VAR
---------
Onceki surumde uzunluk tablosu (HG38 = {...}) DORT ayri dosyaya elle kopyalanmisti ve
sentromer konumlari benim yuvarladigim degerlerdi. Referans degistiginde dortunu birden
guncellemek gerekiyordu; birini unutmak sessiz bir hata demekti. Simdi hepsi .fai'den ve
gercek censat anotasyonundan okunuyor.

REFERANS: T2T-CHM13v2.0 (chr1-22 + chrX, "chr" onekli).
hg38 DEGIL: hg38'in %3.95'i N (assembly boslugu) ve simule edilen genom bu N bloklarini
fiziksel olarak tasidigi icin badread onlardan da okuma uretiyordu -- olculdu: hizalanmayan
okumalarin %86'si poli-N, chr21'de okumalarin %15'i. Gercek veride poli-N okuma OLMAZ.
CHM13v2.0'da N sayisi SIFIR (olculdu).
"""
import os

PHYLO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(PHYLO, "ref", "chm13.primary.fa")
CENSAT = os.path.join(PHYLO, "ref", "chm13_censat_v2.1.bed")
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMS = AUTOSOMES + ["chrX"]

_L = None
_C = None


def lengths(fai=None):
    """kromozom -> uzunluk (bp), referansin .fai'sinden"""
    global _L
    if fai is None and _L is not None:
        return _L
    path = (fai or REF + ".fai")
    d = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            d[f[0]] = int(f[1])
    if fai is None:
        _L = d
    return d


def centromeres(bed=None):
    """kromozom -> sentromer orta noktasi (bp).

    Aktif HOR (higher-order repeat) alfa-satellit dizilerinin kapsadigi araligin ortasi.
    'dhor' (divergent HOR) satirlari haric tutulur: onlar aktif kinetokor bolgesinin
    disindaki eski, bozunmus dizilerdir ve araligi yapay olarak genisletirler.
    """
    global _C
    if bed is None and _C is not None:
        return _C
    path = bed or CENSAT
    span = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 4 or "hor" not in f[3] or "dhor" in f[3]:
                continue
            c, s, e = f[0], int(f[1]), int(f[2])
            if c not in span:
                span[c] = [s, e]
            else:
                span[c][0] = min(span[c][0], s)
                span[c][1] = max(span[c][1], e)
    d = {c: (v[0] + v[1]) // 2 for c, v in span.items()}
    if bed is None:
        _C = d
    return d


def hap_ref_bp(fai=None):
    """haploid referans uzunlugu (chr1-22 + chrX) -- okuma derinliginin paydasi"""
    L = lengths(fai)
    return sum(L[c] for c in CHROMS if c in L)


def centromere_fraction(chrom, fai=None, bed=None):
    """sentromerin kromozom uzunluguna orani (surukleme-dayanikli kol siniri icin)"""
    L = lengths(fai)
    C = centromeres(bed)
    if chrom not in L or chrom not in C:
        return 0.5
    return C[chrom] / L[chrom]


if __name__ == "__main__":
    L, C = lengths(), centromeres()
    print(f"referans     : {REF}")
    print(f"kontig       : {len(L)}")
    print(f"haploid (1-22+X): {hap_ref_bp():,} bp")
    print(f"{'chrom':8s} {'uzunluk':>13s} {'sentromer':>13s} {'oran':>7s}")
    for c in CHROMS:
        if c in L:
            print(f"{c:8s} {L[c]:13,} {C.get(c,0):13,} {centromere_fraction(c):7.3f}")
