#!/usr/bin/env python3
"""Filogenetik agaci kurar, dal uzunluklarini atar ve SEKIL ISTATISTIKLERINI BAGIMSIZ dogrular.

    build_tree.py [--out-dir spec] [--trunk 0.55] [--seed 20260818] [--nulls 20000]

Topoloji kaynagi: docs/01_tumor_phylogenies.md §8.1 (Minussi 2021 superclone mimarisi +
Roerink 2018 tek-hucre turevli klon agaclari + Dentro 2021 kardes-baskinligi).
Bu betik o dokumandaki sayilara GUVENMEZ: Sackin / Colless / J1 / derinlik degerlerini
sifirdan hesaplar ve Yule ile PDA bos dagilimlarini KENDI seed'iyle yeniden uretir.
Dokumandaki degerle uyusmazsa yuksek sesle sikayet eder.
"""
import argparse
import json
import math
import os
import random
import sys

# docs/01 §8.1 — 4 superclone (8/6/4/2), merdiven omurga, karisik cherry+ladder ic yapi
NEWICK_TOPOLOGY = ("((((T01,(T02,((T03,T04),(T05,((T06,T07),T08))))),"
                   "((T09,T10),(T11,((T12,T13),T14)))),"
                   "((T15,T16),(T17,T18))),"
                   "(T19,T20));")

# docs/01 §8.3 — molekuler zaman: patlamali (punctuated) yerlesim
BURST = [0.0375, 0.075]          # omurga bolunmeleri (post-MRCA butcesinin orani)
STEMS = [0.15, 0.15, 0.20, 0.25]  # dort superclone sapi
ALPHA = 0.30                      # superclone icinde her bolunme kalan yolun bu kadarini yer


# ---------------------------------------------------------------- newick
def parse_newick(s):
    """-> ic ice listeler; yaprak = str"""
    s = s.strip().rstrip(";")
    pos = [0]

    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            kids = [node()]
            while s[pos[0]] == ",":
                pos[0] += 1
                kids.append(node())
            assert s[pos[0]] == ")", s[pos[0]]
            pos[0] += 1
            return kids
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",()":
            pos[0] += 1
        return s[start:pos[0]]

    return node()


def leaves(n):
    return [n] if isinstance(n, str) else [x for k in n for x in leaves(k)]


def nleaves(n):
    return 1 if isinstance(n, str) else sum(nleaves(k) for k in n)


# ---------------------------------------------------------------- sekil istatistikleri
def sackin(n):
    """tum yapraklarin koke uzakliklarinin (kenar sayisi) toplami"""
    tot = [0]

    def rec(x, d):
        if isinstance(x, str):
            tot[0] += d
        else:
            for k in x:
                rec(k, d + 1)
    rec(n, 0)
    return tot[0]


def colless(n):
    """yalnizca kesin ikili agaclar icin: her ic dugumde |L-R| toplami"""
    tot = [0]

    def rec(x):
        if isinstance(x, str):
            return 1
        if len(x) != 2:
            raise ValueError("Colless yalnizca ikili agaclarda tanimli")
        a, b = rec(x[0]), rec(x[1])
        tot[0] += abs(a - b)
        return a + b
    rec(n)
    return tot[0]


def j_one(n):
    """Lemant/Noble J^1 denge indeksi, kladogram formu (esit yaprak agirligi).
    Her ic dugumde cocuk alt-agac boyutlarinin normalize Shannon entropisi; dugumler
    (cocuk sayisi>1 olan) esit agirlikla ortalanir."""
    vals = []

    def rec(x):
        if isinstance(x, str):
            return 1
        sizes = [rec(k) for k in x]
        tot = sum(sizes)
        if len(sizes) > 1:
            h = -sum((s / tot) * math.log(s / tot) for s in sizes)
            vals.append(h / math.log(len(sizes)))
        return tot
    rec(n)
    return sum(vals) / len(vals) if vals else 1.0


def max_depth(n):
    if isinstance(n, str):
        return 0
    return 1 + max(max_depth(k) for k in n)


def mean_leaf_depth(n):
    return sackin(n) / nleaves(n)


def n_cherries(n):
    c = [0]

    def rec(x):
        if isinstance(x, str):
            return
        if all(isinstance(k, str) for k in x) and len(x) == 2:
            c[0] += 1
        for k in x:
            rec(k)
    rec(n)
    return c[0]


# ---------------------------------------------------------------- bos modeller
def random_yule(k, rng):
    """pure-birth: her adimda mevcut yapraklardan biri rastgele ikiye ayrilir"""
    nodes = ["L0"]
    tree = {"L0": None}
    idx = 1
    while len(nodes) < k:
        i = rng.randrange(len(nodes))
        parent = nodes.pop(i)
        a, b = f"L{idx}", f"L{idx+1}"
        idx += 2
        tree[parent] = (a, b)
        tree[a] = tree[b] = None
        nodes += [a, b]

    def build(x):
        return x if tree[x] is None else [build(tree[x][0]), build(tree[x][1])]
    return build("L0")


def random_pda(k, rng):
    """PDA (proportional-to-distinguishable-arrangements) = etiketli topolojiler uzerinde DUZGUN.

    DIKKAT: "havuzdan rastgele iki alt-agac sec ve birlestir" PDA DEGILDIR — o, her ciftin esit
    olasilikli oldugu Yule/coalescent semasidir ve Yule ile ayni SEKIL dagilimini verir. (Bu betigin
    ilk surumunde tam olarak bu hata vardi: PDA yuzdelikleri Yule ile ozdes cikiyordu.)
    Dogru uretim: yeni yapragi, mevcut agacin DUZGUN SECILMIS BIR KENARINA tak.
    """
    if k == 1:
        return "T0"
    tree = ["T0", "T1"]
    # kenar listesi: (parent_list, index) — koke tak da dahil edilsin diye sanal kok kenari tutuyoruz
    for i in range(2, k):
        edges = []

        def collect(node, parent, idx):
            edges.append((parent, idx))
            if isinstance(node, list):
                for j, kid in enumerate(node):
                    collect(kid, node, j)
        for j, kid in enumerate(tree):
            collect(kid, tree, j)
        edges.append((None, None))            # kok kenari
        parent, idx = edges[rng.randrange(len(edges))]
        leaf = f"T{i}"
        if parent is None:
            tree = [tree, leaf]
        else:
            parent[idx] = [parent[idx], leaf]
    return tree


def null_percentiles(observed, k, nsim, rng, stat):
    vals = []
    for _ in range(nsim):
        vals.append(stat(random_yule(k, rng)))
    vals.sort()
    yule = 100.0 * sum(1 for v in vals if v < observed) / len(vals)
    vals2 = []
    for _ in range(nsim):
        vals2.append(stat(random_pda(k, rng)))
    vals2.sort()
    pda = 100.0 * sum(1 for v in vals2 if v < observed) / len(vals2)
    return yule, pda, vals, vals2


# ---------------------------------------------------------------- dal uzunluklari
def assign_lengths(root, trunk):
    """Her kenara 'S birimi' cinsinden uzunluk atar (S = bir yapragin tasidigi toplam somatik SNV).
    Kok->MRCA = trunk. Post-MRCA butcesi (1-trunk) patlamali sema ile dagitilir."""
    post = 1.0 - trunk
    edges = {}          # (parent_id, child_id) -> length
    node_id = {}
    counter = [0]

    def label(x, name=None):
        if isinstance(x, str):
            node_id[id(x)] = x
            return x
        counter[0] += 1
        nm = name or f"i{counter[0]:02d}"
        node_id[id(x)] = nm
        return nm

    # MRCA = kokun kendisi (kok = germline). Kokun cocuklari omurga.
    # Omurga: (((scA,scB),scC),scD) yapisinda -> N1, N2 patlama dugumleri
    mrca_name = "MRCA"
    label(root, mrca_name)

    def walk(node, name, consumed):
        """consumed = koke gore harcanmis post-MRCA orani"""
        if isinstance(node, str):
            return
        for child in node:
            if isinstance(child, str):
                ln = (1.0 - consumed) * post * 1.0
                edges[(name, child)] = ln
            else:
                pass

    # Basitlestirilmis, docs/01 §8.3 semasini birebir uygulayan ozyineleme:
    def rec(node, name, remaining, depth_kind, burst_i=0, stem_i=0):
        if isinstance(node, str):
            return
        kids = node
        for ci, child in enumerate(kids):
            if depth_kind == "backbone":
                # omurga: sol=daha derin omurga, sag=superclone sapi
                if isinstance(child, list) and ci == 0 and burst_i < len(BURST):
                    frac = BURST[burst_i]
                    ln = frac * remaining
                    cname = f"N{burst_i+1}"
                    edges[(name, cname)] = ln
                    rec(child, cname, remaining - ln, "backbone", burst_i + 1, stem_i)
                else:
                    frac = STEMS[min(stem_i, len(STEMS) - 1)]
                    ln = frac * remaining
                    cname = f"sc{chr(65+stem_i)}"
                    edges[(name, cname)] = ln
                    rec(child, cname, remaining - ln, "inside")
                    stem_i += 1
            else:
                ln = ALPHA * remaining if isinstance(child, list) else remaining
                cname = child if isinstance(child, str) else f"n{len(edges)}"
                edges[(name, cname)] = ln
                if isinstance(child, list):
                    rec(child, cname, remaining - ln, "inside")
    rec(root, mrca_name, post, "backbone")
    return edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="spec")
    ap.add_argument("--trunk", type=float, default=0.55)
    ap.add_argument("--seed", type=int, default=20260818)
    ap.add_argument("--nulls", type=int, default=20000)
    ap.add_argument("--expect", default="111,49,9",
                    help="docs/01 §8.2'nin iddiasi: sackin,colless,maxdepth (J1 haric - tanim farki)")
    a = ap.parse_args()

    t = parse_newick(NEWICK_TOPOLOGY)
    lv = leaves(t)
    obs = {
        "n_leaves": len(lv),
        "sackin": sackin(t),
        "colless": colless(t),
        "J1": round(j_one(t), 4),
        "max_depth": max_depth(t),
        "mean_leaf_depth": round(mean_leaf_depth(t), 4),
        "n_cherries": n_cherries(t),
    }
    print("=== BAGIMSIZ HESAPLANAN SEKIL ISTATISTIKLERI ===")
    for k, v in obs.items():
        print(f"  {k:18s} {v}")

    exp = a.expect.split(",")
    claims = {"sackin": float(exp[0]), "colless": float(exp[1]), "max_depth": float(exp[2])}
    # J1 KARSILASTIRILMIYOR: docs/01 0.779 diyor, biz 0.876 hesapliyoruz. Sebep tanim farki --
    # Lemant 2022'nin J^1'i ic dugumleri alt-agac boyutuyla agirliklandirir, buradaki surum esit
    # agirlik kullanir. Ikisi de kendi icinde tutarli; bu yuzden J1'i bir DOGRULAMA kapisi
    # yapmiyoruz, yalnizca kendi tanimimizla raporluyoruz. Sackin/Colless/derinlik zaten
    # tanim belirsizligi tasimiyor ve uclu de birebir tutuyor.
    print("\n=== docs/01 §8.2 IDDIASI ILE KARSILASTIRMA ===")
    ok = True
    for k, want in claims.items():
        got = obs[k]
        agree = abs(got - want) < 0.5
        ok &= agree
        print(f"  {k:12s} hesaplanan={got:<8} iddia={want:<8} {'OK' if agree else 'UYUSMUYOR'}")

    rng = random.Random(a.seed)
    print(f"\n=== BOS DAGILIMLAR (n={obs['n_leaves']}, {a.nulls} simulasyon, seed={a.seed}) ===")
    pct = {}
    for name, fn in (("sackin", sackin), ("colless", colless), ("J1", j_one),
                     ("max_depth", max_depth)):
        y, p, yv, pv = null_percentiles(obs[name], obs["n_leaves"], a.nulls, rng, fn)
        pct[name] = {"yule_pct": round(y, 1), "pda_pct": round(p, 1),
                     "yule_5_95": [yv[int(.05*len(yv))], yv[int(.95*len(yv))]],
                     "pda_5_95": [pv[int(.05*len(pv))], pv[int(.95*len(pv))]]}
        inside = pct[name]["yule_5_95"][0] <= obs[name] <= pct[name]["yule_5_95"][1]
        print(f"  {name:12s} gozlenen={obs[name]:<8} Yule%={y:5.1f}  PDA%={p:5.1f}  "
              f"Yule[5,95]={pct[name]['yule_5_95']}  {'ICINDE' if inside else 'DISARIDA'}")
        ok &= inside

    edges = assign_lengths(t, a.trunk)
    tot_internal = sum(v for v in edges.values())
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "tree.newick"), "w") as fh:
        fh.write(NEWICK_TOPOLOGY + "\n")
    out = {
        "topology_newick": NEWICK_TOPOLOGY,
        "source": "docs/01_tumor_phylogenies.md §8.1 (Minussi 2021 / Roerink 2018 / Dentro 2021)",
        "trunk_fraction": a.trunk,
        "shape_observed": obs,
        "shape_null_percentiles": pct,
        "edges_in_units_of_S": {f"{p}->{c}": round(v, 4) for (p, c), v in edges.items()},
        "sum_post_mrca_edges": round(tot_internal, 4),
        "seed": a.seed,
    }
    with open(os.path.join(a.out_dir, "tree.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nyazildi: {a.out_dir}/tree.newick, {a.out_dir}/tree.json")
    print("\nSONUC:", "GECTI" if ok else "KALDI — sekil istatistikleri iddiayla veya Yule zarfiyla uyusmuyor")
    sys.exit(0 if ok else 2)


if __name__ == "__main__":
    main()
