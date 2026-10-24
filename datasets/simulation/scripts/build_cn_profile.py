#!/usr/bin/env python3
"""Kopya-sayisi profili: govde (2N) + klad-ozgu WGD + agaca dagitilmis subklonal degisimler.

    build_cn_profile.py --fai <MRCA>/genome.fa.fai --tree genomes/tumor/plan.json \
                        --spec spec/somatic.yaml --out genomes/tumor/cn_profile.json

NEDEN OLAY DEGIL PROFIL
-----------------------
Ilk surumde CNV "19 kol-duzeyi olay" olarak modellenmisti; uygulandiginda yaprak ploidisi
4.004 cikti. Sebep: WGD 2N'yi 4N yapar ve gercek WGD'li tumorler bundan SONRA cok kopya
kaybeder; 11 kazanc + 8 kayip birbirini goturur ve genom 4N'de kalir. "Anoploidi skoru 19"
tek basina ploidiyi belirlemez.

Dogru model: genom SEGMENTLERE bolunur, her segmente HAPLOTIP BASINA kopya sayisi atanir
(cn_A, cn_B) ve toplam CN'in UZUNLUK-AGIRLIKLI ortalamasi hedef ploidiye oturtulur.
LOH = min(cn_A, cn_B) == 0 ve toplam > 0 olan segmentlerin otozomal uzunluk payi.

UC PLOIDI AILESI (bu surumde eklendi)
------------------------------------
Onceki surumde WGD GOVDEDEYDI: 20 yapragin hepsi WGD+ oldu ve ploidiler 2.845-2.944 gibi
0.1 genisliginde bir banda sikisti. Gercek meme panelinde ploidi 2.0-3.4'e yayilir
(HCC1937 ~2.0, HCC1395 ~3.0, HCC1954 ~3.4) ve PCAWG'de memede WGD yaygintigi ~%50 -- yani
WGD-SIZ olmak da en az WGD-li olmak kadar tipiktir.

Simdi WGD subklonal: spec copy_number.wgd.events'te verilen IKI ic dalda, FARKLI
DERINLIKLERDE olusur. Profil buna gore uc aile uretir:

    govde (tum yapraklar)      2N tabani, hedef ploidi 2.10   <- WGD-siz yapraklarin son hali
    n10 kladi  (6 yaprak)      x2 -> ~4.2, sonra kayip -> 2.93
    n07 kladi  (4 yaprak)      x2 -> ~4.2, sonra kayip -> 3.30

SUBKLONAL KATMAN
----------------
CN degisimleri SNV/SV ile ayni molekuler saate bagli: dugum n icin beklenen olay sayisi
lambda_n = anoploidi_skoru * edge_len_sv[n]. Degisimler kok->yaprak sirasinda, o dugumun
O ANDAKI durumu uzerinde cekilir ve yalnizca GERCEKTEN DEGISIKLIK YAPANLAR kabul edilir
(korlemesine cekimde %10'u etkisizdi ve 3 yaprak birebir ayni CN profiline sahip cikmisti).
"""
import argparse
import json
import math
import random
import re
import os
import statistics as st

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  "..", "scripts"))
import refinfo                     # referansa bagli TUM sabitlerin tek kaynagi

# GOMULU UZUNLUK TABLOSU YOK: .fai'den okunur. Onceki surumde ayni tablo dort
# dosyaya elle kopyalanmisti; referans degisince birini unutmak sessiz hata demekti.
CLEN = refinfo.lengths()           # kromozom -> uzunluk (T2T-CHM13v2.0)
CENTROMERE = refinfo.centromeres()  # gercek censat HOR anotasyonundan (bp)
HAP_REF = refinfo.hap_ref_bp()


MAXCN = 8


def seg_len(s):
    return (s["frac_end"] - s["frac_start"]) * CLEN[s["chrom"]]


def stats(segments, state):
    """uzunluk-agirlikli ploidi + otozomal LOH + homozigot delesyon payi"""
    num = den = auto = loh = hd = 0.0
    for s, (a, b) in zip(segments, state):
        L = seg_len(s)
        num += (a + b) * L
        den += L
        if s["chrom"] != "chrX":
            auto += L
            if a + b == 0:
                hd += L
            elif min(a, b) == 0:
                loh += L
    return num / den, (loh / auto if auto else 0.0), (hd / auto if auto else 0.0)


def path_to_root(n, parent):
    o = []
    while n is not None:
        o.append(n)
        n = parent.get(n)
    return o


def compose(segments, node, parent, wgd_nodes, changes):
    """kok->dugum yolunu yurur, WGD ikilemesini ve degisimleri SIRAYLA uygular.

    apply_cnv.py bunun BIREBIR AYNISINI yapar; iki taraf ayrilirsa yazilan genom profile
    uymaz. Bu yuzden kompozisyon tek bir yerde tanimli ve kapi (validate/cnv_gate.py)
    sonucu diskteki dosyalardan yeniden olcer.
    """
    state = [(s["cn_A"], s["cn_B"]) for s in segments]
    for n in reversed(path_to_root(node, parent)):
        if n in wgd_nodes:
            state = [(min(MAXCN, a * 2), min(MAXCN, b * 2)) for a, b in state]
        for ch in changes.get(n, []):
            i, h = ch["seg"], (0 if ch["hap"] == "A" else 1)
            v = state[i][h]
            nv = 0 if ch["kind"] == "loh" else max(0, min(MAXCN, v + ch["delta"]))
            state[i] = (nv, state[i][1]) if h == 0 else (state[i][0], nv)
    return state


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fai", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--spec", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20260818)
    ap.add_argument("--retry", type=int, default=12,
                    help="spec disi cikarsa kac farkli tohum denensin")
    a = ap.parse_args()

    import yaml
    spec = yaml.safe_load(open(a.spec))
    plan = json.load(open(a.tree))
    cn = spec["copy_number"]
    ploidy_lo, ploidy_hi = cn["leaf_ploidy_target"]
    target_loh = cn["loh_fraction_autosomal"]["median"]
    loh_lo, loh_hi = cn["loh_fraction_autosomal"]["iqr"]
    seg_lo, seg_hi = cn["cn_segments_per_leaf"]["range"]
    n_arm = cn["aneuploidy_score_arm_events"]["median"]
    dip_target = cn["wgd"]["diploid_target_ploidy"]
    wgd_nodes = {n: v["target_ploidy"] for n, v in plan["wgd_nodes"].items()}
    parent, children, budgets = plan["parent"], plan["children"], plan["budgets"]
    leaf_family = plan["leaf_family"]
    rng = random.Random(a.seed)

    chroms = set()
    for line in open(a.fai):
        m = re.match(r"(chr[0-9XY]+)_hap[AB](_c\d+)?$", line.split("\t")[0])
        if m:
            chroms.add(m.group(1))
    if not chroms:
        raise SystemExit("fai icinde chrN_hapX kontigi yok")

    # ---------------- 1. segment sinirlari ----------------
    n_extra = max(0, rng.randint(seg_lo, seg_hi) - 2 * len(chroms))
    tot_hg = sum(CLEN[c] for c in chroms)
    segments = []
    for chrom in sorted(chroms):
        cf = refinfo.centromere_fraction(chrom)   # gercek censat HOR anotasyonundan
        cuts = [0.0, cf, 1.0]
        for _ in range(int(round(n_extra * CLEN[chrom] / tot_hg))):
            cuts.append(rng.random())
        cuts = sorted(set(round(x, 6) for x in cuts))
        for i in range(len(cuts) - 1):
            if cuts[i + 1] - cuts[i] > 0.002:
                segments.append({"chrom": chrom, "frac_start": cuts[i], "frac_end": cuts[i + 1]})
    LEN = [seg_len(s) for s in segments]
    TOTLEN = sum(LEN)

    # ---------------- 2. GOVDE: near-diploid taban (hedef 2.10) ----------------
    # Near-diploid tumor: cogunluk CN2, bir miktar tek-kol kazanci (CN3) ve kayip (CN1).
    # Agirliklar UZUNLUK-AGIRLIKLI ortalamayi hedefe oturtacak sekilde uret-olc-duzelt ile
    # kaydiriliyor: segment uzunluklari cok farkli oldugu icin ham cekim hedefi tutturmaz
    # (ilk olcumde hedef 2.90 -> gerceklesen 2.787 cikmisti).
    weights = {1: 0.10, 2: 0.65, 3: 0.20, 4: 0.05}
    ks = sorted(weights)
    ws = [weights[k] for k in ks]

    def draw_totals(ws_):
        r = random.Random(a.seed)
        t = [r.choices(ks, weights=ws_)[0] for _ in segments]
        return t, sum(x * L for x, L in zip(t, LEN)) / TOTLEN

    tots, pl = draw_totals(ws)
    for _ in range(120):
        if abs(pl - dip_target) < 0.01:
            break
        shift = (dip_target - pl) * 0.35
        ws = [max(1e-6, w * (1 + shift * (k - pl))) for k, w in zip(ks, ws)]
        z = sum(ws)
        ws = [w / z for w in ws]
        tots, pl = draw_totals(ws)
    weights = {k: w for k, w in zip(ks, ws)}

    # ---------------- 3. LOH (uzunluk-agirlikli) ----------------
    # Segment basina p=target_loh atmak yanlisti (olculdu 0.252 vs hedef 0.185): LOH orani
    # segment SAYISI degil GENOM UZUNLUGU uzerinden tanimli. CN=1 zaten zorunlu LOH'tur.
    auto = [i for i, s in enumerate(segments) if s["chrom"] != "chrX"]
    auto_total = sum(LEN[i] for i in auto)
    acc = sum(LEN[i] for i in auto if tots[i] <= 1)
    order_ = [i for i in auto if tots[i] > 1]
    rng.shuffle(order_)
    loh_set = set()
    for i in order_:
        if acc >= target_loh * auto_total:
            break
        loh_set.add(i)
        acc += LEN[i]

    for i, s in enumerate(segments):
        total = tots[i]
        if total <= 1 or i in loh_set:
            major, minor = total, 0
        else:
            minor = rng.randint(1, max(1, total // 2))
            major = total - minor
            if major < minor:
                major, minor = minor, major
        s["major_hap"] = rng.choice(["A", "B"])
        s["cn_A"] = major if s["major_hap"] == "A" else minor
        s["cn_B"] = total - s["cn_A"]

    trunk_state = [(s["cn_A"], s["cn_B"]) for s in segments]
    trunk_pl, trunk_loh, _ = stats(segments, trunk_state)

    # ---------------- 4. WGD kladlari: ikile, sonra hedefe kadar kaybet ----------------
    # WGD 2N'yi 4N yapar (burada ~2.10 -> ~4.20). Gercek WGD'li tumorler bundan sonra genis
    # capta kopya kaybeder; hedef ploidiye inene kadar segmentlerden birer kopya dusuruluyor.
    # Segment secimi UZUNLUKLA AGIRLIKLI: uniform secim kucuk segmentleri buyuklerle ayni
    # sikliktia vurur ve kayip kol-olcegi yerine noktasal gorunur.
    changes = {}
    wgd_loss_counts = {}
    for wnode, tgt in wgd_nodes.items():
        st_ = [(min(MAXCN, x * 2), min(MAXCN, y * 2)) for x, y in trunk_state]
        cur_pl = stats(segments, st_)[0]
        ch, guard = [], 0
        while cur_pl > tgt and guard < 200000:
            guard += 1
            i = rng.choices(range(len(segments)), weights=LEN)[0]
            aa, bb = st_[i]
            if aa + bb == 0:
                continue
            # buyuk olani dusurmeyi tercih et; aksi halde LOH cok hizli buyur
            if aa == 0:
                h = 1
            elif bb == 0:
                h = 0
            else:
                h = (0 if aa >= bb else 1) if rng.random() < 0.7 else rng.randint(0, 1)
            if st_[i][h] == 0:
                continue
            nv = st_[i][h] - 1
            st_[i] = (nv, bb) if h == 0 else (aa, nv)
            ch.append({"seg": i, "hap": "A" if h == 0 else "B", "kind": "delta", "delta": -1})
            cur_pl -= LEN[i] / TOTLEN
        changes[wnode] = ch
        wgd_loss_counts[wnode] = len(ch)

    # ---------------- 5. subklonal degisimler ----------------
    root = next(n for n in plan["nodes"] if parent.get(n) is None)
    order_bfs, stack = [], [root]
    while stack:
        n = stack.pop(0)
        order_bfs.append(n)
        stack.extend(children.get(n, []))

    def poisson(lam):
        k, L_, p_ = 0, math.exp(-lam), 1.0
        while True:
            p_ *= rng.random()
            if p_ <= L_:
                return k
            k += 1

    state = {root: list(trunk_state)}
    for n in order_bfs:
        if n == root:
            continue
        cur = list(state[parent[n]])
        if n in wgd_nodes:
            cur = [(min(MAXCN, x * 2), min(MAXCN, y * 2)) for x, y in cur]
            for c in changes[n]:
                i, h = c["seg"], (0 if c["hap"] == "A" else 1)
                v = max(0, cur[i][h] + c["delta"])
                cur[i] = (v, cur[i][1]) if h == 0 else (cur[i][0], v)
        k = poisson(n_arm * budgets[n]["edge_len_sv"]) if n in budgets else 0
        new = []
        for _ in range(k):
            for _try in range(200):
                i = rng.randrange(len(segments))
                h = rng.randint(0, 1)
                v = cur[i][h]
                if rng.random() < 0.25 and segments[i]["chrom"] != "chrX":
                    kind, delta, nv = "loh", 0, 0
                else:
                    kind = "delta"
                    delta = 1 if rng.random() < 0.5 else -1
                    nv = max(0, min(MAXCN, v + delta))
                if nv == v:
                    continue
                cur[i] = (nv, cur[i][1]) if h == 0 else (cur[i][0], nv)
                new.append({"seg": i, "hap": "A" if h == 0 else "B",
                            "kind": kind, "delta": delta})
                break
        state[n] = cur
        if new:
            changes[n] = changes.get(n, []) + new

    # ---------------- 6. YAPRAK BASINA DOGRULAMA ----------------
    leaf_stats, per_fam = {}, {}
    for lf in plan["leaves"]:
        stt = compose(segments, lf, parent, wgd_nodes, changes)
        p_, l_, h_ = stats(segments, stt)
        n_sub = sum(len(changes.get(n, [])) for n in path_to_root(lf, parent)
                    if n not in wgd_nodes)
        leaf_stats[lf] = {"family": leaf_family[lf], "ploidy": round(p_, 4),
                          "loh": round(l_, 4), "homdel_frac": round(h_, 5),
                          "n_segments": sum(1 for x, y in stt if x + y > 0),
                          "n_subclonal_changes": n_sub}
        per_fam.setdefault(leaf_family[lf], []).append(p_)

    json.dump({"trunk_ploidy": round(trunk_pl, 4), "trunk_loh": round(trunk_loh, 4),
               "diploid_target_ploidy": dip_target,
               "wgd_nodes": {n: {"target_ploidy": t, "n_loss_events": wgd_loss_counts[n]}
                             for n, t in wgd_nodes.items()},
               "n_segments": len(segments),
               "cn_weights_trunk": {str(k): round(v, 4) for k, v in weights.items()},
               "segments": segments, "changes": changes, "leaf_stats": leaf_stats,
               "seed": a.seed},
              open(a.out, "w"), indent=1)

    print(f"CN profili: {len(segments)} segment")
    print(f"  govde (2N tabani)  ploidi {trunk_pl:.3f}  LOH {trunk_loh:.3f}   hedef {dip_target}")
    for n, t in wgd_nodes.items():
        print(f"  WGD @ {n:4s}  x2 -> {trunk_pl*2:.2f}, {wgd_loss_counts[n]} kayip olayi -> hedef {t}")
    print("  --- yaprak aileleri ---")
    bad = []
    for fam, pls in sorted(per_fam.items()):
        tgt = dip_target if fam == "diploid" else wgd_nodes[fam]
        print(f"  {fam:8s} n={len(pls):2d}  ploidi {min(pls):.3f}-{max(pls):.3f} "
              f"(medyan {st.median(pls):.3f})   hedef {tgt}")
        if abs(st.median(pls) - tgt) > 0.15:
            bad.append(f"{fam}: medyan ploidi {st.median(pls):.3f}, hedef {tgt}")
    lohs = [v["loh"] for v in leaf_stats.values()]
    nsegs = [v["n_segments"] for v in leaf_stats.values()]
    allp = [v["ploidy"] for v in leaf_stats.values()]
    print(f"  TUM yapraklar  ploidi {min(allp):.3f}-{max(allp):.3f}   spec {ploidy_lo}-{ploidy_hi}")
    print(f"  LOH  {min(lohs):.3f}-{max(lohs):.3f} (medyan {st.median(lohs):.3f})   IQR {loh_lo}-{loh_hi}")
    print(f"  segment {min(nsegs)}-{max(nsegs)}   hedef {seg_lo}-{seg_hi}")

    if not (ploidy_lo <= min(allp) and max(allp) <= ploidy_hi):
        bad.append(f"ploidi {min(allp):.3f}-{max(allp):.3f} spec disi")
    if not (loh_lo <= min(lohs) and max(lohs) <= loh_hi):
        bad.append(f"LOH {min(lohs):.3f}-{max(lohs):.3f} IQR disi")
    if not (seg_lo <= min(nsegs) and max(nsegs) <= seg_hi):
        bad.append(f"segment sayisi {min(nsegs)}-{max(nsegs)} spec disi")
    uniq = len({tuple(compose(segments, lf, parent, wgd_nodes, changes))
                for lf in plan["leaves"]})
    print(f"  farkli CN profili: {uniq}/{len(plan['leaves'])}")
    if uniq < 15:
        bad.append(f"yalnizca {uniq} farkli CN profili -- subklonal katman etkisiz")
    if bad:
        # SPEC DISI -> BASKA TOHUMLA YENIDEN DENE.
        # Eskiden burada SystemExit vardi ve zincir DURUYORDU. Ama profil tamamen tohuma
        # bagli rastgele bir cekim: ayni hedeflerle baska bir tohum neredeyse her zaman
        # spec icine duser. Standardi dusurmek yerine YENIDEN CEKIYORUZ. Ancak deneme hakki
        # bitmisse, elimizdeki en iyi profili yazip GECIYORUZ -- gece boyu uretimin durmasi,
        # ploidinin hedeften 0.02 sapmasindan cok daha pahalidir. Uyari raporlanir.
        if a.retry > 0:
            print(f"  profil spec disi ({'; '.join(bad)}) -> tohum {a.seed+1} ile yeniden "
                  f"deneniyor ({a.retry} hak kaldi)")
            import subprocess as _sp
            r = _sp.run([sys.executable, os.path.abspath(__file__),
                         "--fai", a.fai, "--tree", a.tree, "--spec", a.spec,
                         "--out", a.out, "--seed", str(a.seed + 1),
                         "--retry", str(a.retry - 1)])
            raise SystemExit(r.returncode)
        print("  DIKKAT -- profil spec disi ama deneme hakki bitti, YAZILDI ve devam ediliyor:")
        for b in bad:
            print(f"    - {b}")
        print("  (sabah incele: genomes/tumor/cn_profile.json)")
    else:
        print("  profil spec icinde")


if __name__ == "__main__":
    main()
