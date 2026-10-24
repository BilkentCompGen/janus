#!/usr/bin/env python3
"""Bir yaprak genomuna kopya-sayisi profilini uygular (kontig cerrahisi, akis halinde).

    apply_cnv.py --node T01 --tumor-dir genomes/tumor \
                 --profile genomes/tumor/cn_profile.json [--out-name genome_cnv.fa]

MODEL
-----
Kopya sayisi = o lokusu kapsayan KONTIG sayisi. Okuma simulatoru kontiglerden uzunlukla
orantili okuma urettigi icin CN kendiliginden dogru derinlige donusur; ayri bir "CN modeli"
gerekmez.

Girdi genom, yapragin WGD durumuna gore 2N ya da 4N'dir:

    WGD-siz yaprak   chrN_hapA              chrN_hapB           (haplotip basina 1 kaynak)
    WGD'li yaprak    chrN_hapA_c1/_c2       chrN_hapB_c1/_c2    (haplotip basina 2 kaynak)

Cikti her iki durumda da "iz" (track) kontigleridir: haplotip A icin cn_A tane, B icin cn_B
tane iz; her iz, o haplotipte kendisine yer olan segmentlerin birlesimi.

KAYNAK KOPYA SECIMI -- DUZELTILEN HATA
--------------------------------------
Onceki surumde iz t'nin kaynagi t'nin PARITESINDEN geliyordu (t0->c1, t1->c2, t2->c1...).
Sonuc: CN=1'e dusen her bolgede hayatta kalan HEP c1 oluyordu. Olculdu: ciktinin %71'i c1,
%29'u c2 -- yani WGD'den SONRA c2 uzerinde olusan mutasyonlar sistematik olarak siliniyor,
c1'dekiler hep kurtuluyordu. Gercekte hangi kopyanin kaybedildigi segment segment rastgeledir.
Simdi kaynaklar HER SEGMENT ICIN ayri karistiriliyor: cn=1 ise c1/c2'den biri rastgele,
cn=2 ise ikisi de, cn=3 ise ikisi + rastgele biri tekrar.

LOH = bir haplotipin o segmentteki TUM izlerinin yoklugu (cn_A==0 veya cn_B==0). O bolgedeki
germline heterozigot varyantlar okuma duzeyinde homozigot gorunur -- gercek LOH mekanizmasi
tam olarak budur.

KOORDINAT SAGLAMLIGI: profil sinirlari KROMOZOM KESRI olarak tanimli; her kaynak kontigin KENDI
uzunluguna gore cozulur. Boylece dugumde birikmis indel/SV suruklemesi (~yuz kb) Mb olcegindeki
segment sinirlarinda anlamsiz kalir.

BELLEK: dizi hicbir zaman tumuyle bellege alinmaz; samtools faidx ciktisi parca parca okunup
sabit genislikte yeniden sarilarak yaziliyor (aksi halde kontig basina ~250 MB).
"""
import argparse
import collections
import json
import os
import random
import re
import subprocess
import sys

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  "..", "scripts"))
import refinfo                     # referansa bagli TUM sabitlerin tek kaynagi

# GOMULU UZUNLUK TABLOSU YOK: .fai'den okunur. Onceki surumde ayni tablo dort
# dosyaya elle kopyalanmisti; referans degisince birini unutmak sessiz hata demekti.
CLEN = refinfo.lengths()           # kromozom -> uzunluk (T2T-CHM13v2.0)
HAP_REF = refinfo.hap_ref_bp()


MAXCN = 8
WIDTH = 60


class FastaWriter:
    """sabit genislikte satirlarla yazan akis sarici (samtools faidx indeksleyebilsin diye)"""

    def __init__(self, fh):
        self.fh = fh
        self.carry = ""

    def open_contig(self, name):
        self.flush()
        self.fh.write(f">{name}\n")

    def feed(self, chunk):
        s = self.carry + chunk
        k = (len(s) // WIDTH) * WIDTH
        if k:
            self.fh.write("\n".join(s[i:i + WIDTH] for i in range(0, k, WIDTH)))
            self.fh.write("\n")
        self.carry = s[k:]

    def flush(self):
        if self.carry:
            self.fh.write(self.carry + "\n")
            self.carry = ""


def stream_region(src, region, writer):
    """samtools faidx ciktisini basligi atarak sariciya akitir; yazilan baz sayisini doner"""
    p = subprocess.Popen(["samtools", "faidx", src, region],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                         bufsize=1 << 20)
    wrote, first = 0, True
    while True:
        buf = p.stdout.read(1 << 20)
        if not buf:
            break
        if first:
            buf = buf[buf.find("\n") + 1:]
            first = False
        buf = buf.replace("\n", "")
        writer.feed(buf)
        wrote += len(buf)
    err = p.stderr.read()
    if p.wait() != 0:
        raise RuntimeError(f"samtools faidx {region} basarisiz: {err}")
    return wrote


def path_to_root(n, parent):
    o = []
    while n is not None:
        o.append(n)
        n = parent.get(n)
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--node", required=True)
    ap.add_argument("--tumor-dir", required=True)
    ap.add_argument("--profile", required=True)
    ap.add_argument("--out-name", default="genome_cnv.fa")
    ap.add_argument("--seed", type=int, default=20260818)
    a = ap.parse_args()

    tumor = os.path.abspath(a.tumor_dir)
    plan = json.load(open(os.path.join(tumor, "plan.json")))
    prof = json.load(open(os.path.abspath(a.profile)))
    parent = plan["parent"]
    wgd_nodes = set(plan["wgd_nodes"])
    node_dir = os.path.join(tumor, a.node)
    src = os.path.join(node_dir, "genome.fa")
    if not os.path.exists(src + ".fai"):
        subprocess.run(["samtools", "faidx", src], check=True)
    L = {}
    for line in open(src + ".fai"):
        f = line.split("\t")
        L[f[0]] = int(f[1])

    # --- kaynak kopyalar: (kromozom, haplotip) -> [kontig, ...] ---
    # WGD-siz yaprakta 1, WGD'li yaprakta 2 kaynak olur.
    srcs = collections.defaultdict(list)
    for c in L:
        m = re.match(r"(chr[0-9XY]+)_hap([AB])(?:_c(\d+))?$", c)
        if m:
            srcs[(m.group(1), m.group(2))].append(c)
    for k in srcs:
        srcs[k].sort()

    # --- bu yapragin CN durumu: kok yolu boyunca WGD ikilemesi + degisimler ---
    # build_cn_profile.compose() ile BIREBIR AYNI mantik; ikisi ayrilirsa yazilan genom
    # profile uymaz ve bunu validate/cnv_gate.py diskteki dosyalardan yakalar.
    segs = prof["segments"]
    changes = prof["changes"]
    state = [(s["cn_A"], s["cn_B"]) for s in segs]
    n_wgd = 0
    for n in reversed(path_to_root(a.node, parent)):
        if n in wgd_nodes:
            state = [(min(MAXCN, x * 2), min(MAXCN, y * 2)) for x, y in state]
            n_wgd += 1
        for ch in changes.get(n, []):
            i, h = ch["seg"], (0 if ch["hap"] == "A" else 1)
            v = state[i][h]
            nv = 0 if ch["kind"] == "loh" else max(0, min(MAXCN, v + ch["delta"]))
            state[i] = (nv, state[i][1]) if h == 0 else (state[i][0], nv)

    by_chrom = {}
    for i, s in enumerate(segs):
        by_chrom.setdefault(s["chrom"], []).append(i)
    for c in by_chrom:
        by_chrom[c].sort(key=lambda i: segs[i]["frac_start"])

    # --- segment basina kaynak kopya atamasi (parite DEGIL, rastgele) ---
    rng = random.Random(a.seed ^ (abs(hash(a.node)) % 999983))
    pick = {}
    for i, s in enumerate(segs):
        for hi, hap in ((0, "A"), (1, "B")):
            avail = srcs.get((s["chrom"], hap), [])
            need = state[i][hi]
            if not avail or need == 0:
                pick[(i, hap)] = []
                continue
            perm = list(avail)
            rng.shuffle(perm)
            pick[(i, hap)] = [perm[t % len(perm)] for t in range(need)]

    out_fa = os.path.join(node_dir, a.out_name)
    bed = os.path.join(node_dir, "cn_segments.bed")
    total_bp = 0
    written_cn = [0] * len(segs)
    src_bp = collections.Counter()

    with open(out_fa, "w") as fo, open(bed, "w") as fb:
        w = FastaWriter(fo)
        fb.write("#out_contig\tout_start\tout_end\tchrom\tfrac_start\tfrac_end\t"
                 "src_contig\tsrc_start\tsrc_end\tcn_A\tcn_B\n")
        for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
            idxs = by_chrom[chrom]
            for hi, hap in ((0, "A"), (1, "B")):
                mx = max(state[i][hi] for i in idxs)
                for t in range(mx):
                    name = f"{chrom}_hap{hap}_t{t}"
                    members = [i for i in idxs if state[i][hi] > t]
                    if not members:
                        continue
                    w.open_contig(name)
                    pos = 0
                    for i in members:
                        srcc = pick[(i, hap)][t]
                        Ls = L[srcc]
                        s0 = int(segs[i]["frac_start"] * Ls)
                        s1 = int(segs[i]["frac_end"] * Ls)
                        if s1 <= s0:
                            continue
                        n = stream_region(src, f"{srcc}:{s0+1}-{s1}", w)
                        fb.write(f"{name}\t{pos}\t{pos+n}\t{chrom}\t"
                                 f"{segs[i]['frac_start']}\t{segs[i]['frac_end']}\t"
                                 f"{srcc}\t{s0}\t{s1}\t{state[i][0]}\t{state[i][1]}\n")
                        pos += n
                        total_bp += n
                        written_cn[i] += 1
                        src_bp["c" + srcc.split("_c")[-1] if "_c" in srcc else "tek"] += n
        w.flush()

    # ---------------- BAGIMSIZ DOGRULAMA ----------------
    # Planlanan CN ile CIKTIYA GERCEKTEN YAZILAN kopya sayisini karsilastir. Bu, profilden
    # degil yazma dongusunun sayacindan gelir; atlanmis/bos kalan bir iz burada yakalanir.
    mism = [(segs[i]["chrom"], segs[i]["frac_start"], state[i][0] + state[i][1], written_cn[i])
            for i in range(len(segs)) if written_cn[i] != state[i][0] + state[i][1]]

    auto_bp = loh_bp = hd_bp = 0.0
    for i, s in enumerate(segs):
        if s["chrom"] == "chrX":
            continue
        base = (s["frac_end"] - s["frac_start"]) * CLEN[s["chrom"]]
        auto_bp += base
        a_, b_ = state[i]
        if a_ + b_ == 0:
            hd_bp += base
        elif min(a_, b_) == 0:
            loh_bp += base

    ls = prof.get("leaf_stats", {}).get(a.node, {})
    rep = {"node": a.node, "family": ls.get("family"),
           "wgd": n_wgd > 0, "n_wgd_on_path": n_wgd,
           "n_segments": sum(1 for x, y in state if x + y > 0),
           "ploidy_written": round(total_bp / HAP_REF, 4),
           "ploidy_planned": ls.get("ploidy"),
           "loh_fraction": round(loh_bp / auto_bp, 4),
           "homdel_fraction": round(hd_bp / auto_bp, 5),
           "total_bp": total_bp,
           "src_copy_balance": {k: round(v / max(1, total_bp), 4) for k, v in src_bp.items()},
           "cn_mismatches": len(mism)}
    if mism:
        rep["cn_mismatch_examples"] = mism[:5]
    subprocess.run(["samtools", "faidx", out_fa], check=True)
    rep["n_out_contigs"] = sum(1 for _ in open(out_fa + ".fai"))
    json.dump(rep, open(os.path.join(node_dir, "cnv_report.json"), "w"), indent=1)

    bal = " ".join(f"{k}={v:.0%}" for k, v in sorted(rep["src_copy_balance"].items()))
    print(f"{a.node}: aile={str(rep['family']):8s} WGD={rep['wgd']} "
          f"kontig={rep['n_out_contigs']} {total_bp/1e9:.3f} Gb "
          f"PLOIDI={rep['ploidy_written']:.3f} (plan {rep['ploidy_planned']}) "
          f"LOH={rep['loh_fraction']:.3f} kaynak[{bal}] CN-uyusmazlik={len(mism)}")
    if mism:
        print(f"  !! {len(mism)} segmentte yazilan kopya sayisi plandan farkli", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
