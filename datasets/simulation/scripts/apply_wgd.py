#!/usr/bin/env python3
"""Tum-genom duplikasyonu: her kontigi ikiye kopyalar (46 -> 92 kontig).

    apply_wgd.py --in <MRCA_early>/genome.fa --out <MRCA_wgd>/genome.fa

NEDEN AYRI VE NEDEN AGACIN ORTASINDA
------------------------------------
WGD govde molekuler zamaninin %80'inde olur (spec: Roerink 2018 P2 %84, CI 80.8-87.6;
Bielski 2018: timeable vakalarin %97.3'unde TP53 mutasyonundan sonra). Bunun tek onemli
sonucu ploidi degil, **mutasyon yasi**dir:

    WGD'den ONCE olusan klonal mutasyon  -> ikilenir  -> CN=4'te 2 kopya -> VAF ~0.50
    WGD'den SONRA olusan klonal mutasyon -> ikilenmez -> CN=4'te 1 kopya -> VAF ~0.25

Gerstung 2020'nin "erken klonal / gec klonal" ayrimi tam olarak bu sinyaldir. WGD'yi agacin
SONUNDA uygularsak butun klonal mutasyonlar WGD'den once olmus sayilir, hepsi ayni VAF'a
duser ve sinyal YOK OLUR. Bu yuzden WGD, govdeyi ikiye bolen ayri bir adimdir.

Kontig adlandirma: `chr1_hapA` -> `chr1_hapA_c1`, `chr1_hapA_c2`. Sonraki dugumlerdeki
mutasyonlar bu kopyalardan BIRINE duser; kopya secimi ornekleyici tarafindan dogal olarak
yapilir (kontig = fiziksel kromatid).

Akis halinde calisir (samtools faidx), boylece 12 Gb'lik cikti icin tepe bellek tek kontig
mertebesinde kalir.
"""
import argparse
import os
import subprocess


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    a = ap.parse_args()

    if not os.path.exists(a.inp + ".fai"):
        subprocess.run(["samtools", "faidx", a.inp], check=True)
    contigs = []
    with open(a.inp + ".fai") as fh:
        for line in fh:
            f = line.split("\t")
            contigs.append((f[0], int(f[1])))

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    total = 0
    with open(a.out, "w") as fo:
        for name, ln in contigs:
            r = subprocess.run(["samtools", "faidx", a.inp, name],
                               capture_output=True, text=True, check=True)
            body = r.stdout.split("\n", 1)[1]
            for cp in (1, 2):
                fo.write(f">{name}_c{cp}\n")
                fo.write(body if body.endswith("\n") else body + "\n")
                total += ln
    subprocess.run(["samtools", "faidx", a.out], check=True)
    print(f"WGD: {len(contigs)} -> {len(contigs)*2} kontig, {total/1e9:.3f} Gb")
    print(f"     her kontig X -> X_c1 + X_c2 (ozdes); bundan SONRAKI mutasyonlar")
    print(f"     kopyalardan yalnizca BIRINE duser -> CN=4'te VAF ~0.25")


if __name__ == "__main__":
    main()
