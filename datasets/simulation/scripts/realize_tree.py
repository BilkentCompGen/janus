#!/usr/bin/env python3
"""FAZ D — agaci dugum dugum gercege dokme.

    realize_tree.py --tree spec/tree.json --spec spec/somatic.yaml \
                    --germline genomes/germline/germline.fa \
                    --out genomes/tumor --seed 20260818 [--only NODE] [--plan-only]

Calisma bicimi
--------------
Kokten yapraga, her dugum icin:
  1. EBEVEYNIN genomu okunur (kok icin germline.fa),
  2. o dugumun KENDI dal butcesi kadar mutasyon, ebeveynin GERCEK dizisi uzerinde uretilir
     (trinukleotid baglamlari ebeveynin dizisinden okunur -- referanstan degil),
  3. mutasyonlar acik POS + acik GT ile bir VCF'e yazilir,
  4. insilicoSV `haploid: True` + `import:` ile calistirilir,
  5. kontig adlarindaki fazladan `_hapA` soneki soyulur,
  6. dizi duzeyinde dogrulama yapilir.

Neden ebeveyn uzerinde uretiyoruz
---------------------------------
Gercek evrimde bir subklonal mutasyon, atalarinin degistirdigi dizinin uzerine duser. Ayrica
insilicoSV'nin VCF'i GIRDI koordinatlarindadir; dugumun VCF'ini ebeveyn koordinatlarinda
uretmek, zincirin dogal koordinat sozlesmesidir. hg38'e projeksiyon ayri bir adimda
(PAF zinciri + bagimsiz minimap2 geri-hizalamasi) yapilir.

Zigosite
--------
Kontig = haplotip. Bir mutasyon TEK bir kontige yazilir, dolayisiyla yapisal olarak
heterozigottur. `pick_genotype` hicbir zaman cagrilmaz (bkz. experiments/01).
"""
import argparse
import json
import os
import random
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genome_sampler as GS          # noqa: E402
import mutsig                        # noqa: E402

ISV = os.path.join(os.path.dirname(sys.executable), "insilicosv")


# ------------------------------------------------------------------ agac
def parse_newick(s):
    s = s.strip().rstrip(";")
    pos = [0]

    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            kids = [node()]
            while s[pos[0]] == ",":
                pos[0] += 1
                kids.append(node())
            pos[0] += 1
            return kids
        st = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",()":
            pos[0] += 1
        return s[st:pos[0]]
    return node()


def label_tree(root):
    """-> (nodes, parent, children, leafnames). Ic dugumler n00, n01, ... kok = MRCA."""
    parent, children, order = {}, {}, []
    counter = [0]

    def rec(x, pname):
        if isinstance(x, str):
            name = x
        else:
            name = "MRCA" if pname is None else f"n{counter[0]:02d}"
            counter[0] += 1
        parent[name] = pname
        order.append(name)
        children[name] = []
        if pname is not None:
            children[pname].append(name)
        if not isinstance(x, str):
            for k in x:
                rec(k, name)
        return name
    rec(root, None)
    leaves = [n for n in order if not children[n]]
    return order, parent, children, leaves


def edge_budgets(order, parent, children, trunk, burst, stems, alpha):
    """Her kenara 'S birimi' uzunluk atar. docs/01 §8.3 semasi.
    Kok kenari (germline -> MRCA) = trunk. Post-MRCA butcesi patlamali dagitilir."""
    length = {"MRCA": trunk}
    remaining = {"MRCA": 1.0 - trunk}
    # omurga: MRCA'nin cocuklari sirayla (ilk = daha derin omurga, digerleri superclone sapi)
    stem_i = [0]
    burst_i = [0]

    def rec(node, kind):
        kids = children[node]
        if not kids:
            return
        for ci, ch in enumerate(kids):
            if kind == "backbone" and ci == 0 and children[ch] and burst_i[0] < len(burst):
                f = burst[burst_i[0]]
                burst_i[0] += 1
                length[ch] = f * remaining[node]
                remaining[ch] = remaining[node] - length[ch]
                rec(ch, "backbone")
            elif kind == "backbone":
                f = stems[min(stem_i[0], len(stems) - 1)]
                stem_i[0] += 1
                length[ch] = f * remaining[node]
                remaining[ch] = remaining[node] - length[ch]
                rec(ch, "inside")
            else:
                if children[ch]:
                    length[ch] = alpha * remaining[node]
                    remaining[ch] = remaining[node] - length[ch]
                    rec(ch, "inside")
                else:
                    length[ch] = remaining[node]     # terminal dal kalani yer
                    remaining[ch] = 0.0
    rec("MRCA", "backbone")
    return length


# ------------------------------------------------------------------ FASTA
def load_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.rstrip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, seqs, width=60):
    with open(path, "w") as fh:
        for n, s in seqs.items():
            fh.write(f">{n}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


# ------------------------------------------------------------------ mutasyon uretimi
def gen_snvs(gc, n, chans, cum, forbidden, rng):
    """imza-kosullu SNV uretimi -> [(contig, pos0, ref, alt, channel)]"""
    if n <= 0:
        return []
    # kanal -> kac adet
    counts = np.zeros(len(chans), dtype=np.int64)
    draws = rng.random(n)
    idxs = np.searchsorted(np.array(cum), draws)
    idxs = np.clip(idxs, 0, len(chans) - 1)
    for i in idxs:
        counts[i] += 1
    need = {}
    ch_by_code = {}
    for i, (ch, tri, ref, alt) in enumerate(chans):
        if counts[i] == 0:
            continue
        code = GS.TRI2CODE[tri]
        need[code] = need.get(code, 0) + int(counts[i])
        ch_by_code.setdefault(code, []).append((ch, ref, alt, int(counts[i])))
    got, _ = gc.sample_by_context(need, forbidden)
    out = []
    for code, positions in got.items():
        k = 0
        for ch, ref, alt, cnt in ch_by_code[code]:
            for _ in range(cnt):
                if k >= len(positions):
                    break
                cname, p = positions[k]
                k += 1
                # .upper() SART: CHM13 yumusak maskeli. Kucuk harf gelirse base != ref
                # HER ZAMAN dogru cikar ve site "ters iplikte" sayilir -> ALT tumlenir.
                # OLCULDU (2026-08-19 04:37): MRCA'nin 6,190 kaydinin %47.1'inde REF kucuk
                # harfliydi; bir kismi duz "c>C" gibi NULL MUTASYON, kalani YANLIS KANAL.
                # S09 imza kapisi bunu goremedi cunku kaydedilen SBSCH etiketini okuyor,
                # gercekte yazilan ref/alt'i degil -> niyeti dogruluyor, sonucu degil.
                base = gc.seqs[cname][p].upper()
                if base == ref:
                    r_o, a_o = ref, alt
                else:
                    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
                    r_o, a_o = base, comp.get(alt, alt)
                out.append((cname, p, r_o, a_o, ch))
    return out


def gen_indels(gc, n, id_spec, forbidden, rng, pyrng):
    """ID imzalarina gore indel uretimi -> [(contig,pos0,'INS'|'DEL',length,seq)]"""
    if n <= 0:
        return []
    chans = list(id_spec.keys())
    probs = np.array([id_spec[c] for c in chans], dtype=float)
    probs /= probs.sum()
    picks = rng.choice(len(chans), size=n, p=probs)
    sites = gc.sample_uniform(n, forbidden, min_gap=2)
    out = []
    for i, (cname, p) in enumerate(sites):
        ch = chans[picks[i]]
        parts = ch.split(":")
        L = mutsig.id_channel_to_length(ch, pyrng)
        seq = gc.seqs[cname]
        if parts[1] == "Ins":
            # eklenen dizi: kanal 'C'/'T' ise homopolimer, 'R'/'M' ise komsu dizinin kopyasi
            if parts[2] in ("C", "T"):
                s = parts[2] * L
            else:
                s = seq[p:p + L].upper()
            # .upper() SART: CHM13v2.0 YUMUSAK MASKELI -- tekrar bolgeleri kucuk harf
            # (chr21'in %62.4'u). Bu kontrol kucuk harfi "gecersiz baz" sayip indeli
            # DUSURUYORDU. Olculdu: planlanan 501 indel -> gerceklesen 277 (%55). Ustelik
            # kayip RASTGELE DEGIL: tam da indellerin YOGUNLASTIGI tekrar bolgelerinden
            # oluyordu (ID imzalarinin cogu homopolimer/STR kaynaklidir), yani hem sayi
            # hem DAGILIM bozuluyordu.
            if not s or any(b not in "ACGT" for b in s.upper()):
                continue
            out.append((cname, p, "INS", L, s))
        else:
            if p + L + 1 >= len(seq):
                continue
            if any(b not in "ACGT" for b in seq[p:p + L + 1].upper()):
                continue
            out.append((cname, p, "DEL", L, None))
    return out


def gen_svs(gc, n, sv_mix, forbidden, rng, pyrng):
    """Somatik SV uretimi -- ALTI SINIF (yalniz DEL/DUP DEGIL).

    Sinif paylari `spec/somatic.yaml sv_class_mix_measured.fractions`'tan gelir ve bunlar
    UC MEME HUCRE HATTININ GERCEK truth setlerinden olculmustur (CASTLE minda_4_2, 3,085 SV):
        BND %40.6 · DEL %36.5 · DUP %12.5 · INS %10.4

    -> [(contig,pos0,type,length,extra)] ; extra: INS icin dizi, TRA icin (hedef_kontig,hedef_pos)

    Uretilen tipler ve insilicoSV karsiliklari (hepsi import ile ampirik dogrulandi):
        DEL        <DEL>
        DUP        <DUP>    tandem
        INV        <INV>
        INS        <INS>    dizi verilir (mobil element benzeri uzunluk dagilimi)
        TRA_unbal  <dDUP>   TARGET_CHROM farkli kontig -> dengesiz translokasyon
        TRA_recip  <rTRA>   TARGET_CHROM farkli kontig -> resiprok translokasyon
    """
    if n <= 0:
        return []
    names = list(sv_mix.keys())
    w = np.array([sv_mix[k] for k in names], dtype=float)
    w /= w.sum()
    picks = rng.choice(len(names), size=n, p=w)
    sites = gc.sample_uniform(int(n * 2.5) + 50, forbidden, min_gap=2)
    out, si = [], 0

    def lognorm(med, sigma, lo, hi):
        return int(min(hi, max(lo, pyrng.lognormvariate(np.log(med), sigma))))

    for i in range(n):
        if si >= len(sites):
            break
        kind = names[picks[i]]
        cname, p = sites[si]
        si += 1
        seq = gc.seqs[cname]
        L = 0
        extra = None

        if kind == "DEL":
            # CASTLE meme hatlarinda delesyon boyutu genis: medyan ~2 kb, kuyruk Mb'lere
            L = lognorm(2000, 1.6, 50, 5_000_000)
        elif kind == "DUP":
            L = lognorm(12000, 1.5, 200, 3_000_000)
        elif kind == "INV":
            L = lognorm(60000, 1.5, 500, 5_000_000)
        elif kind == "INS":
            # mobil element benzeri: Alu ~300, SVA ~1500, L1 ~6000 modlari
            mode = pyrng.choices([300, 1500, 6000], weights=[0.6, 0.15, 0.25])[0]
            L = lognorm(mode, 0.35, 50, 12000)
            extra = "".join(pyrng.choice("ACGT") for _ in range(L))
        elif kind in ("TRA_unbal", "TRA_recip"):
            # rTRA hedefte de AYNI uzunlukta bir aralik kaplar; cok buyuk secilirse iki
            # tarafta birden catisma olasiligi patliyor. Gercek resiprok translokasyonlar da
            # genelde daha kucuk parcalar tasir. Ust sinir dusuruldu.
            L = lognorm(30000, 1.3, 1000, 400_000)
            # HEDEF: BASKA bir kontig (interkromozomal) -- gercek BND'lerin cogu boyle
            others = [c for c in gc.names if c.split("_")[0] != cname.split("_")[0]]
            if not others:
                continue
            tc = pyrng.choice(others)
            tlen = gc.lengths[tc]
            tp = None
            for _ in range(40):
                cand = pyrng.randrange(1_000_000, max(1_000_001, tlen - 1_000_000))
                if gc.codes[tc][cand] < 0:          # N bolgesi
                    continue
                tf = forbidden.setdefault(tc, set())
                if any((cand + d) in tf for d in (-2, -1, 0, 1, 2)):
                    continue
                tp = cand
                tf.add(cand)
                break
            if tp is None:
                continue
            extra = (tc, tp)
        else:
            continue

        if p + L + 2 >= len(seq):
            continue
        # N kontrolu TUM acikligi ornekler. Onceki surum yalnizca ilk 200 kb'ye bakiyordu;
        # 676 kb'lik bir rTRA sentromere uzanip N'e girebiliyordu (chr2_hapA ~90 Mb'de oldu).
        step = max(1, L // 40)
        bad_n = False
        for q in range(p, min(p + L, len(seq)), step):
            w = seq[q:q + 1000]
            if w.count("N") + w.count("n") > 0.01 * max(len(w), 1):
                bad_n = True
                break
        if bad_n:
            continue
        out.append((cname, p, kind, L, extra))

    return out


def gen_complex_cluster(gc, n_events, forbidden, rng, pyrng):
    """Kromotripsis benzeri KUME: dar bir pencerede zincirlenmis DEL/INV/DUP yigini.

    Gercek kromotripsis, tek bir kromozom bolgesinde onlarca kirilma noktasi ve salinan
    kopya sayisi demektir. Dagitik tekil SV'lerle taklit edilemez; bu yuzden ayri uretilir.
    Kaynak: spec complex_events.chromothripsis (Cortes-Ciriano 2020: bolge medyan 56.6 Mb,
    26 SV, 30 CN segmenti).
    """
    if n_events <= 0:
        return []
    big = [c for c in gc.names if gc.lengths[c] > 60_000_000]
    if not big:
        return []
    cname = pyrng.choice(big)
    ln = gc.lengths[cname]
    span = min(int(56.6e6), ln - 4_000_000)
    start = pyrng.randrange(1_000_000, max(1_000_001, ln - span - 1_000_000))
    seq = gc.seqs[cname]
    out = []
    forb = forbidden.setdefault(cname, set())
    for _ in range(n_events):
        p = pyrng.randrange(start, start + span)
        kind = pyrng.choices(["DEL", "INV", "DUP"], weights=[0.55, 0.3, 0.15])[0]
        L = int(min(2_000_000, max(1000, pyrng.lognormvariate(np.log(30000), 1.3))))
        if p + L + 2 >= len(seq) or any((p + d) in forb for d in (-1, 0, 1)):
            continue
        win = seq[p:p + min(L, 100000)]
        if win.count("N") + win.count("n") > 0.01 * max(len(win), 1):
            continue
        forb.add(p)
        out.append((cname, p, kind, L, None))
    return out


# ------------------------------------------------------------------ VCF
def _up(x):
    """Allel/dizi alanlarini buyuk harfe cevirir.

    Son savunma hatti: yumusak maskeli bir referansta genomdan turetilen her dizi kucuk
    harf gelebilir. insilicoSV'nin kendi dogrulayicisi da bunu reddediyor
    ("INSSEQ has to be a valid sequence, aa was provided"), yani sessiz bozulma degil
    gurultulu bir cokme olur -- ama yine de burada temizleniyor ki hicbir yol kacamasin.
    """
    return x.upper() if isinstance(x, str) else x


def write_node_vcf(path, contig_lengths, snvs, indels, svs, node="NA", tsv_path=None):
    """VCF + mutations.tsv yazar.

    KIMLIK KURALI [REF 6 R2.2]: her mutasyon, uretildigi anda atanan ve BIR DAHA DEGISMEYEN bir
    kimlik tasir (`<dugum>:<sinif><sira>`), ve bu kimlik VCF'in ID alanina yazilir. Kimlik ASLA
    koordinattan turetilmez: bir dugumun koordinatlari, yol uzerindeki her indel/SV'den sonra
    kayar, dolayisiyla koordinat-anahtarli kume karsilastirmasi ikinci dugumde coker.
    """
    recs = []
    for i, (cname, p, r, a, ch) in enumerate(snvs):
        mid = f"{node}:S{i}"
        recs.append((cname, p + 1, p + 1, [cname, str(p + 1), mid, _up(r), _up(a), "100", "PASS",
                                           f"MUTCLASS=SNV;SBSCH={ch};MUTID={mid}", "GT", "1|1"]))
    for i, (cname, p, t, L, sq) in enumerate(indels):
        mid = f"{node}:I{i}"
        if t == "INS":
            recs.append((cname, p + 1, p + 1,
                         [cname, str(p + 1), mid, "N", "<INS>", "100", "PASS",
                          f"SVTYPE=INS;END={p+1};SVLEN={L};INSSEQ={_up(sq)};MUTCLASS=INDEL;MUTID={mid}",
                          "GT", "1|1"]))
        else:
            recs.append((cname, p + 1, p + L,
                         [cname, str(p + 1), mid, "N", "<DEL>", "100", "PASS",
                          f"SVTYPE=DEL;END={p+L};SVLEN={L};MUTCLASS=INDEL;MUTID={mid}", "GT", "1|1"]))
    for i, (cname, p, t, L, extra) in enumerate(svs):
        mid = f"{node}:V{i}"
        if t == "INS":
            recs.append((cname, p + 1, p + 1,
                         [cname, str(p + 1), mid, "N", "<INS>", "100", "PASS",
                          f"SVTYPE=INS;END={p+1};SVLEN={L};INSSEQ={_up(extra)};"
                          f"MUTCLASS=SV;SVCLASS=INS;MUTID={mid}", "GT", "1|1"]))
        elif t in ("TRA_unbal", "TRA_recip"):
            # insilicoSV: dDUP = dengesiz (kaynak korunur, hedefe kopya), rTRA = resiprok takas.
            # TARGET_CHROM farkli kontig oldugu icin INTERKROMOZOMAL -> gercek BND imzasi.
            tc, tp = extra
            isv = "dDUP" if t == "TRA_unbal" else "rTRA"
            recs.append((cname, p + 1, p + L,
                         [cname, str(p + 1), mid, "N", f"<{isv}>", "100", "PASS",
                          f"SVTYPE={isv};END={p+L};SVLEN={L};TARGET={tp};TARGET_CHROM={tc};"
                          f"MUTCLASS=SV;SVCLASS={t};MUTID={mid}", "GT", "1|1"]))
        else:                                   # DEL / DUP / INV
            recs.append((cname, p + 1, p + L,
                         [cname, str(p + 1), mid, "N", f"<{t}>", "100", "PASS",
                          f"SVTYPE={t};END={p+L};SVLEN={L};MUTCLASS=SV;SVCLASS={t};MUTID={mid}",
                          "GT", "1|1"]))
    # ORTUSME TEMIZLIGI -- kaynak aralik VE translokasyon HEDEF lokusu birlikte.
    # Ilk surumde yalnizca kaynak araligi bakiliyordu; bir dDUP'in BASKA kontigteki hedefi
    # orada duran bir SV'nin icine dusunce insilicoSV tum kosuyu dusuruyordu
    # ("cannot place imported SV ... <dDUP> ... TARGET_CHROM=chr21_hapB").
    recs.sort(key=lambda r: (r[0], r[1], r[2]))
    occupied = {}                     # kontig -> sirali [(s,e)]
    for r in recs:
        occupied.setdefault(r[0], [])

    def _hits(c, s0, e0):
        for (a0, b0) in occupied.get(c, []):
            if s0 <= b0 and a0 <= e0:
                return True
        return False

    kept, last_end, last_c = [], -10, None
    dropped = 0
    for r in recs:
        if r[0] != last_c:
            last_c, last_end = r[0], -10
        if r[1] <= last_end:
            dropped += 1
            continue
        info = r[3][7]
        tgt = None
        if "TARGET_CHROM=" in info:
            tc = info.split("TARGET_CHROM=")[1].split(";")[0]
            tp = int(info.split("TARGET=")[1].split(";")[0])
            # rTRA hedefte AYNI uzunlukta bir aralik kaplar (A_B -> B_A); dDUP icin de
            # eklenen kopya kadar yer tutar. Nokta olarak rezerve etmek yetmiyordu.
            tlen = int(info.split("SVLEN=")[1].split(";")[0]) if "SVLEN=" in info else 2
            if _hits(tc, tp - 2, tp + tlen + 2):
                dropped += 1
                continue
            tgt = (tc, tp, tlen)
        kept.append(r)
        last_end = r[2]
        occupied.setdefault(r[0], []).append((r[1], r[2]))
        if tgt:
            occupied.setdefault(tgt[0], []).append((tgt[1] - 2, tgt[1] + tgt[2] + 2))
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, L in contig_lengths.items():
            fh.write(f"##contig=<ID={c},length={L}>\n")
        for t, n, ty in (("SVTYPE", "1", "String"), ("END", "1", "Integer"),
                         ("SVLEN", "1", "Integer"), ("INSSEQ", "1", "String"),
                         ("TARGET", "1", "Integer"), ("TARGET_CHROM", "1", "String"),
                         ("MUTCLASS", "1", "String"), ("SVCLASS", "1", "String"),
                         ("SBSCH", "1", "String")):
            fh.write(f'##INFO=<ID={t},Number={n},Type={ty},Description="{t}">\n')
        for t in ("DEL", "INS", "DUP", "INV", "dDUP", "rTRA"):
            fh.write(f'##ALT=<ID={t},Description="{t}">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for r in kept:
            fh.write("\t".join(r[3]) + "\n")
    if tsv_path:
        newf = not os.path.exists(tsv_path)
        with open(tsv_path, "a") as tf:
            if newf:
                tf.write("mut_id\tbirth_node\tcontig\tpos_parent\tend_parent\tclass\t"
                         "svtype\tsvlen\tref\talt\tinfo\n")
            for r in kept:
                f = r[3]
                info = f[7]
                cls = info.split("MUTCLASS=")[1].split(";")[0] if "MUTCLASS=" in info else "NA"
                svt = info.split("SVTYPE=")[1].split(";")[0] if "SVTYPE=" in info else "SNV"
                svl = info.split("SVLEN=")[1].split(";")[0] if "SVLEN=" in info else "1"
                tf.write(f"{f[2]}\t{node}\t{f[0]}\t{f[1]}\t{r[2]}\t{cls}\t{svt}\t{svl}\t"
                         f"{f[3]}\t{f[4]}\t{info}\n")
    return len(kept), dropped


# ------------------------------------------------------------------ tek dugum
def realize_node(node, parent_fa, out_dir, budget, gc_cache, spec, chans, cum, seed, log):
    t0 = time.time()
    d = os.path.join(out_dir, node)
    os.makedirs(d, exist_ok=True)
    node_fa = os.path.join(d, "genome.fa")
    if os.path.exists(node_fa) and os.path.getsize(node_fa) > 1000:
        log(f"{node}: zaten var, atlaniyor")
        return node_fa

    seqs = load_fasta(parent_fa)
    log(f"{node}: ebeveyn yuklendi {sum(len(v) for v in seqs.values())/1e9:.3f} Gb, "
        f"{len(seqs)} kontig ({time.time()-t0:.0f}s)")

    rng = np.random.default_rng(seed)
    pyrng = random.Random(seed)
    gc = GS.GenomeContexts(seqs, rng)
    log(f"{node}: baglam kodlari hazir ({time.time()-t0:.0f}s)")

    forbidden = {}
    snvs = gen_snvs(gc, budget["snv"], chans, cum, forbidden, rng)
    indels = gen_indels(gc, budget["indel"], spec["signatures"]["id_mix_expanded"],
                        forbidden, rng, pyrng)
    svs = gen_svs(gc, budget["sv"], spec["per_leaf"]["sv_class_mix_measured"]["fractions"],
                  forbidden, rng, pyrng)
    # kromotripsis: yalnizca govdenin WGD-sonrasi parcasinda, bir kez
    if budget.get("complex_cluster"):
        svs += gen_complex_cluster(gc, budget["complex_cluster"], forbidden, rng, pyrng)
    log(f"{node}: uretildi SNV={len(snvs)} indel={len(indels)} SV={len(svs)} "
        f"({time.time()-t0:.0f}s)")

    vcf = os.path.join(d, "mutations.vcf")
    n_kept, n_drop = write_node_vcf(vcf, {k: len(v) for k, v in seqs.items()},
                                    snvs, indels, svs, node=node,
                                    tsv_path=os.path.join(out_dir, "mutations.tsv"))
    log(f"{node}: VCF {n_kept} kayit ({n_drop} ortusme elendi)")

    cfg = os.path.join(d, "config.yaml")
    with open(cfg, "w") as fh:
        fh.write(f'reference: "{parent_fa}"\nrandom_seed: {seed}\nhaploid: True\n'
                 # output_adjacencies KAPALI: rTRA import edilince
                 # output_novel_adjacencies IndexError ile cokuyor (olculdu). Adjacency'ler
                 # zaten truth VCF'inden turetilebilir; uretici raporuna ihtiyac yok.
                 f'output_paf: True\n'
                 f'variant_sets:\n    - import: "{vcf}"\n')
    # Binlerce SV icinde TEK bir yerlestirilemeyen kayit tum kosuyu oldururdu. insilicoSV
    # hata mesajinda kaydin ID'sini veriyor; onu dusurup tekrar deniyoruz. Her dusurulen
    # kayit LOGLANIR -- sessiz kayip yok.
    import re as _re
    for attempt in range(12):
        r = subprocess.run([ISV, "-c", "config.yaml"], cwd=d, capture_output=True, text=True)
        if r.returncode == 0:
            break
        m = _re.search(r"cannot place imported SV vcf_record:\S+\s+\d+\s+(\S+)", r.stderr)
        if not m:
            log(f"{node}: insilicoSV HATA rc={r.returncode}\n{r.stderr[-1200:]}")
            return None
        bad_id = m.group(1)
        lines = [ln for ln in open(vcf) if ln.startswith("#") or f"\t{bad_id}\t" not in ln]
        with open(vcf, "w") as fh:
            fh.writelines(lines)
        log(f"{node}: yerlestirilemeyen SV dusuruldu ({bad_id}), deneme {attempt+1}")
    else:
        log(f"{node}: 12 denemeden sonra hala yerlestirilemiyor")
        return None

    # --- UYGULAMA DOGRULAMASI: insilicoSV ISTEDIGIMIZ HER MUTASYONU UYGULADI MI? ---
    # NEDEN VAR: mutations.tsv URETECIN NIYETINI kaydeder; hicbir sey o mutasyonun genomda
    # gercekten bulundugunu sormuyordu. Butun asagi akis kapilari (S01 SBS yuku, S02 indel
    # orani, S03 SV yuku, S09 imza) bu NIYET kaydini sayiyor. insilicoSV bir kaydi sessizce
    # atlarsa hepsi yine "PASS" derdi.
    # Ham kayit sayisi karsilastirilamaz: insilicoSV bazi SV'leri ikiye boler (rTRA, dDUP),
    # MRCA'da 6,204 girdi -> 6,286 cikti. Ama MUTID korunuyor ("Imported_0_MRCA:S4661"),
    # dolayisiyla KIMLIK KUMESI karsilastirilir.
    want, got = set(), set()
    try:
        for ln in open(vcf):
            if not ln.startswith("#"):
                want.add(ln.split("\t", 3)[2])
        simvcf = os.path.join(d, "sim.vcf")
        if os.path.exists(simvcf):
            import re as _re2
            for ln in open(simvcf):
                if ln.startswith("#"):
                    continue
                f2 = ln.split("\t")
                cands = [f2[2]]
                # insilicoSV kimligi "Imported_<vset>_<bizim_id>" yapar ve COK PARCALI
                # SV'lerde sona "_0"/"_1" ekler (rTRA -> 2 kayit; MRCA'da 82 rTRA -> 164
                # kayit). Ilk surumde bu son ek soyulmuyordu ve 82 rTRA "uygulanmadi" gibi
                # gorunuyordu -- yani DOGRULAMANIN KENDISI hatali alarm veriyordu.
                m2 = _re2.match(r"Imported_\d+_(.+?)(?:_\d+)?$", f2[2])
                if m2:
                    cands.append(m2.group(1))
                if len(f2) > 7 and "SVID=" in f2[7]:
                    cands.append(f2[7].split("SVID=")[1].split(";")[0].split("\t")[0])
                    m3 = _re2.match(r"Imported_\d+_(.+)$", cands[-1])
                    if m3:
                        cands.append(m3.group(1))
                for w in cands:
                    if w in want:
                        got.add(w)
                        break
        missing = sorted(want - got)
        json.dump({"node": node, "n_requested": len(want), "n_applied": len(got),
                   "n_missing": len(missing), "missing_examples": missing[:10]},
                  open(os.path.join(d, "apply_report.json"), "w"), indent=1)
        if missing:
            log(f"{node}: !! {len(missing)} mutasyon UYGULANMADI (ornek: {missing[:3]})")
        else:
            log(f"{node}: uygulama dogrulandi ({len(want)} mutasyonun hepsi sim.vcf'te)")
    except Exception as e:                                       # noqa: BLE001
        log(f"{node}: uygulama dogrulamasi yapilamadi: {e}")

    raw = os.path.join(d, "sim.fa")
    if not os.path.exists(raw):
        log(f"{node}: sim.fa yok"); return None
    # insilicoSV her kosuda kontig adina '_hapA' EKLIYOR -> soy
    out = {}
    for k, v in load_fasta(raw).items():
        out[k[:-5] if k.endswith("_hapA") else k] = v
    write_fasta(node_fa, out)
    # sim.fa, genome.fa'nin yalnizca kontig adlari farkli bir KOPYASIDIR (~6 GB). 39 dugumde
    # 240 GB gereksiz yer tutar. Silmiyoruz (proje kurali), trash'e tasiyoruz.
    try:
        tdir = os.path.expanduser("~/trash/simulation-simfa")
        os.makedirs(tdir, exist_ok=True)
        os.replace(raw, os.path.join(tdir, f"{node}.sim.fa"))
    except OSError as e:                                        # noqa: BLE001
        log(f"{node}: sim.fa tasinamadi ({e}) -- yer tutmaya devam ediyor")
    log(f"{node}: genom yazildi {sum(len(v) for v in out.values())/1e9:.3f} Gb "
        f"({time.time()-t0:.0f}s)")
    return node_fa


# ------------------------------------------------------------------ main
def expand_id_mix(id_matrix, id_mix):
    """{ID1:0.27,...} -> {kanal: olasilik} (83 kanal). Imza adi bulunamazsa yuksek sesle patlar."""
    missing = [s for s in id_mix if s not in id_matrix]
    if missing:
        raise KeyError(f"ID imza matrisinde yok: {missing}")
    tot = sum(id_mix.values())
    spec = {}
    for s, w in id_mix.items():
        for ch, p in id_matrix[s].items():
            spec[ch] = spec.get(ch, 0.0) + (w / tot) * p
    z = sum(spec.values())
    return {k: v / z for k, v in spec.items() if v > 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="spec/tree.json")
    ap.add_argument("--spec", default="spec/somatic.yaml")
    ap.add_argument("--germline", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20260818)
    ap.add_argument("--only", help="yalnizca bu dugum")
    ap.add_argument("--plan-only", action="store_true")
    a = ap.parse_args()

    # YOLLAR MUTLAK OLMALI: insilicoSV her dugum icin `cwd=<dugum dizini>` ile calisiyor,
    # dolayisiyla config icindeki goreli bir referans/VCF yolu orada cozulemez. Faz C'de ayni
    # sinif hata sessiz bir rc=1 uretmisti; burada "reference must be a readable file" olarak
    # patladi. Bir kez basta mutlaklastiriyoruz.
    a.germline = os.path.abspath(a.germline)
    a.out = os.path.abspath(a.out)
    a.spec = os.path.abspath(a.spec)
    a.tree = os.path.abspath(a.tree)

    import yaml
    spec = yaml.safe_load(open(a.spec))
    tree = json.load(open(a.tree))
    root = parse_newick(tree["topology_newick"])
    order, parent, children, leaves = label_tree(root)

    BURST = [0.0375, 0.075]
    STEMS = [0.15, 0.15, 0.20, 0.25]
    ALPHA = 0.30
    trunk_snv = spec["tree_distribution"]["snv_trunk_fraction"]["value"]
    trunk_sv = spec["tree_distribution"]["sv_trunk_fraction"]["value"]
    len_snv = edge_budgets(order, parent, children, trunk_snv, BURST, STEMS, ALPHA)
    len_sv = edge_budgets(order, parent, children, trunk_sv, BURST, STEMS, ALPHA)

    S_SNV = spec["per_leaf"]["sbs"]["median"]
    S_IND = spec["per_leaf"]["indel"]["median"]
    # SV butcesi JUNCTION duzeyinde (bir cagirici bunu gorur), LINX-kumelenmis OLAY degil.
    # Deger, uc meme hucre hattinin GERCEK truth setlerinden olculdu (medyan 1106).
    S_SV = spec["per_leaf"]["sv_junctions"]["median"]
    SV_SPREAD = spec["per_leaf"]["sv_junctions"].get("per_leaf_lineage_spread", 1.9)

    # ---- dal uzunlugu modeli (spec/parameters.yaml tree.branch_budget) ----
    # REFEREE DUZELTMESI [REF 1 M1]: ilk surumde terminal dallara 5.5x'lik bir lognormal hiz
    # carpani uygulaniyordu. Ama TOPOLOJI ZATEN 5.97x'lik bir private-dal yayilimi uretiyor
    # (0.382/0.064 S birimi); ustune carpan koymak 24-42x'lik, hicbir kaynagin desteklemedigi
    # bir yayilim demekti. Ayrica [REF 1 B3]: 4-7x yayilimin dayandigi iki hucre hatti (LS180
    # MSI, HT115 POLE) HIPERMUTATOR; ayni deneyde HT115'in kardes klonlari yalnizca 1.23x
    # yayiliyor. Biz hipermutator ILAN ETMEDIK.
    #
    # DUZELTILMIS MODEL:
    #   in-vivo private yayilim = YALNIZ topolojik (carpan yok)
    #   + her yaprağa SABIT BEKLENTILI bir kultur artisi C (20 klon AYNI sure kultur edilir)
    #     C ~ lognormal(medyan 230 SBS), kaynak: Kuijk 2020 7.2-8.3 SBS/populasyon ikilenmesi
    #     x ~30 ikilenme  [REF 3 R2.1 ile atif duzeltildi; Blokzijl 2016 DEGIL]
    #   C uzerine soy-ozgu heterojenlik: p05-p95 yayilim 2.0x  [REF 1 B4]
    rng0 = random.Random(a.seed ^ 0x5EED)
    C_MEDIAN = 230
    C_SPREAD = 2.0
    c_sigma = np.log(C_SPREAD) / (2 * 1.6449)
    culture = {n: (C_MEDIAN * np.exp(rng0.gauss(-0.5 * c_sigma ** 2, c_sigma))
                   if not children[n] else 0.0) for n in order}
    # SV'ler icin de soy-ozgu yayilim. ONCEKI SURUMDE YOKTU ve butun yapraklar 271-273
    # SV tasiyordu (yayilim 1.007x) -- gercek veride 20 klon arasinda bu kadar tekduzelik
    # olamaz; tek basina Poisson gurultusu bile +-33 verirdi. Simdi: her soya lognormal bir
    # SV hiz carpani + her kenara Poisson.
    # Soy-ozgu SV hizi: bir SOYUN TAMAMINA uygulanir (kromozomal kararsizlik bir klonun
    # ozelligidir, tek bir kenarin degil). Kokten yaprağa carpanlar carpilir, boylece derin
    # dallarda kararsizlik birikir -- gercek CIN davranisi budur.
    sv_sigma = np.log(SV_SPREAD) / (2 * 1.6449)
    _own = {n: np.exp(rng0.gauss(-0.5 * sv_sigma ** 2, sv_sigma)) for n in order}
    sv_factor = {}

    def _acc(n):
        if n in sv_factor:
            return sv_factor[n]
        p_ = parent.get(n)
        sv_factor[n] = _own[n] if p_ is None else _own[n] * _acc(p_) ** 0.75
        return sv_factor[n]
    for n in order:
        _acc(n)
    sv_factor["MRCA"] = 1.0                # govde ortak: carpan uygulanmaz
    factor = {n: 1.0 for n in order}          # in-vivo dallara carpan UYGULANMIYOR

    # ---------------- WGD: GOVDEDE DEGIL, SUBKLONAL ----------------
    # ONCEKI SURUM: WGD govde molekuler zamaninin %80'inde uygulaniyordu, dolayisiyla 20
    # yapragin HEPSI WGD+ oluyordu ve ploidiler 2.845-2.944 gibi 0.1 genisliginde bir banda
    # sikisiyordu. Gercek meme panelinde ploidi genistir: HCC1937 ~2.0 (near-diploid),
    # HCC1395 ~3.0, HCC1954 ~3.4. PCAWG'de memede WGD yaygintigi ~%50 -- yani WGD-SIZ olmak
    # da en az WGD-li olmak kadar tipiktir. Dewhurst 2014: ayni tumorde WGD+ ve WGD-
    # subklonlar bir arada bulunur.
    #
    # YENI MODEL: WGD, agacin FARKLI DERINLIKLERINDEKI iki ic dalda bagimsiz olusur. Bir
    # dalda olustugunda o daldan inen tum yapraklar WGD+ olur. Erken/gec klonal VAF sinyali
    # KAYBOLMAZ: WGD dalinin USTUNDEKI mutasyonlar ikilenir (CN4'te 2 kopya, VAF~0.50), o
    # daldan itibaren olusanlar tek kopyada kalir (VAF~0.25). Govdeyi ikiye bolmeye gerek
    # yok, cunku WGD zaten bir DUGUM SINIRINDA duruyor.
    wgd_cfg = spec["copy_number"]["wgd"]
    wgd_nodes = {e["node"]: e for e in wgd_cfg["events"]}

    def _depth(n):
        d = 0
        while parent.get(n) is not None:
            n = parent[n]
            d += 1
        return d

    def _leaves_under(n):
        if not children.get(n):
            return [n]
        o = []
        for c in children[n]:
            o += _leaves_under(c)
        return o

    # spec'teki dugum adlari topolojiye bagli; sessizce yanlis yere WGD koymamak icin DOGRULA
    _wgd_leaves, _bad = set(), []
    for nd, ev in wgd_nodes.items():
        if nd not in order:
            _bad.append(f"{nd}: agacta yok")
            continue
        if not children.get(nd):
            _bad.append(f"{nd}: yaprak, ic dugum olmali")
        if _depth(nd) != ev.get("expect_depth", _depth(nd)):
            _bad.append(f"{nd}: derinlik {_depth(nd)}, spec {ev['expect_depth']} diyor")
        lv = set(_leaves_under(nd))
        if lv & _wgd_leaves:
            _bad.append(f"{nd}: baska bir WGD dalinin altinda (ic ice WGD)")
        _wgd_leaves |= lv
    _frac = len(_wgd_leaves) / len(leaves)
    if abs(_frac - wgd_cfg["target_leaf_fraction"]) > 0.15:
        _bad.append(f"WGD+ yaprak orani {_frac:.2f}, hedef {wgd_cfg['target_leaf_fraction']}")
    if len({_depth(n) for n in wgd_nodes}) < len(wgd_nodes):
        _bad.append("WGD olaylari ayni derinlikte -- spec farkli derinlik istiyor")
    if _bad:
        raise SystemExit("WGD yerlesimi gecersiz:\n  - " + "\n  - ".join(_bad))

    budgets = {}
    for n in order:
        if n == "MRCA":
            ls, lv = len_snv["MRCA"], len_sv["MRCA"]
        else:
            ls, lv = len_snv.get(n, 0.0), len_sv.get(n, 0.0)
        budgets[n] = {
            "snv": int(round(ls * S_SNV + culture[n])),
            "indel": int(round(ls * S_IND + culture[n] * S_IND / S_SNV)),
            "sv": int(np.random.default_rng(a.seed + abs(hash(n)) % 99991)
                      .poisson(max(0.0, lv * S_SV * sv_factor[n]))),
            "edge_len_snv": round(ls, 5), "edge_len_sv": round(lv, 5),
            "culture_increment": int(round(culture[n])),
        }

    # kromotripsis: govdede, tek olay (spec complex_events). Govde artik BOLUNMUYOR --
    # WGD subklonal oldugu icin "WGD oncesi/sonrasi govde" ayrimina gerek kalmadi.
    _ct = spec["complex_events"]["chromothripsis"]
    budgets["MRCA"]["complex_cluster"] = int(_ct["n_svs"]["median"])
    budgets["MRCA"]["note"] = ("govde -- tum yapraklara kalitilir; kromotripsis kumesi burada. "
                               "WGD dallarinin USTUNDE oldugu icin bu mutasyonlar o "
                               "kladlarda ikilenir (CN4'te VAF~0.50)")
    for nd, ev in wgd_nodes.items():
        budgets[nd]["wgd"] = True
        budgets[nd]["wgd_target_ploidy"] = ev["target_ploidy"]
        budgets[nd]["note"] = (f"WGD BURADA (derinlik {_depth(nd)}, "
                               f"{len(_leaves_under(nd))} yaprak, hedef ploidi "
                               f"{ev['target_ploidy']}). Bu daldan itibaren olusan "
                               f"mutasyonlar tek kopyada kalir (VAF~0.25).")

    # yaprak toplamlari (kalitilan + kendi) -- kontrol
    def path_to_root(n):
        p, cur = [], n
        while cur is not None:
            p.append(cur)
            cur = parent[cur]
        return p
    leaf_tot = {}
    for lf in leaves:
        tot = {"snv": 0, "indel": 0, "sv": 0}
        for n in path_to_root(lf):
            for k in tot:
                tot[k] += budgets[n][k]
        leaf_tot[lf] = tot

    # her yaprak hangi WGD olayinin altinda? (CN profili ve kapilar bunu kullanir)
    leaf_family = {}
    for lf in leaves:
        fam = "diploid"
        for n in path_to_root(lf):
            if n in wgd_nodes:
                fam = n
                break
        leaf_family[lf] = fam

    plan = {"nodes": order, "parent": parent, "children": children, "leaves": leaves,
            "budgets": budgets, "leaf_totals": leaf_tot,
            "trunk_snv": trunk_snv, "trunk_sv": trunk_sv,
            "wgd_nodes": {n: {"target_ploidy": e["target_ploidy"], "depth": _depth(n),
                              "leaves": _leaves_under(n)} for n, e in wgd_nodes.items()},
            "leaf_family": leaf_family,
            "diploid_target_ploidy": wgd_cfg["diploid_target_ploidy"],
            "S": {"snv": S_SNV, "indel": S_IND, "sv": S_SV}}
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "plan.json"), "w") as fh:
        json.dump(plan, fh, indent=2)

    sn = [v["snv"] for v in leaf_tot.values()]
    sv = [v["sv"] for v in leaf_tot.values()]
    print(f"dugum sayisi {len(order)} (yaprak {len(leaves)}, ic {len(order)-len(leaves)})")
    print(f"govde SNV butcesi     : {budgets['MRCA']['snv']:,}  ({trunk_snv:.0%} of S)")
    for n, e in wgd_nodes.items():
        lv = _leaves_under(n)
        print(f"WGD @ {n:4s} derinlik {_depth(n)}  {len(lv):2d} yaprak "
              f"({lv[0]}..{lv[-1]})  hedef ploidi {e['target_ploidy']}")
    _nd = sum(1 for f in leaf_family.values() if f == "diploid")
    print(f"WGD-siz yaprak        : {_nd}/{len(leaves)}  hedef ploidi "
          f"{wgd_cfg['diploid_target_ploidy']}")
    print(f"yaprak SNV toplami    : medyan {int(np.median(sn)):,}  "
          f"aralik {min(sn):,}-{max(sn):,}  yayilim {max(sn)/min(sn):.2f}x")
    print(f"yaprak SV toplami     : medyan {int(np.median(sv)):,}  aralik {min(sv):,}-{max(sv):,}")
    print(f"agactaki toplam farkli SNV: {sum(b['snv'] for b in budgets.values()):,}")
    if a.plan_only:
        return

    sbs_mat, _ = mutsig.load_signature_matrix(spec["signatures"]["sbs_reference"])
    sbs_mix = {k: v for k, v in spec["signatures"]["sbs_mix"].items()
               if isinstance(v, (int, float))}
    spectrum = mutsig.mix_spectrum(sbs_mat, sbs_mix)
    id_mat, _ = mutsig.load_signature_matrix(spec["signatures"]["id_reference"])
    id_mix = {k: v for k, v in spec["signatures"]["id_mix"].items()
              if isinstance(v, (int, float))}
    spec["signatures"]["id_mix_expanded"] = expand_id_mix(id_mat, id_mix)

    logf = open(os.path.join(a.out, "realize.log"), "a")

    def log(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    # kanal ornekleyici GENOMDAN BAGIMSIZ -- dongu disinda BIR KEZ kurulur.
    # (Ilk surumde her dugumde 6 GB'lik ebeveyn FASTA'si UC KEZ yukleniyordu: bir kez
    # bu kurulum icin, bir kez kullanilmayan bir denemede, bir kez realize_node icinde.
    # 39 dugumde ~470 GB gereksiz disk okumasi demekti.)
    chans, probs = [], []
    for ch, pr in spectrum.items():
        if pr <= 0:
            continue
        left, ref, alt, right = mutsig.parse_sbs_channel(ch)
        tri = left + ref + right
        if tri not in GS.TRI2CODE:
            continue
        chans.append((ch, tri, ref, alt))
        probs.append(pr)
    z = sum(probs)
    cum, acc = [], 0.0
    for pr in probs:
        acc += pr / z
        cum.append(acc)
    log(f"spektrum hazir: {len(chans)} kanal")

    todo = [a.only] if a.only else order
    for node in todo:
        if parent[node] is None:
            pfa = a.germline
        else:
            pfa = os.path.join(a.out, parent[node], "genome.fa")
            if node in wgd_nodes:
                # WGD BU DALDA: ebeveyn genomunun kontigleri once ikilenir, mutasyonlar
                # ikilenmis genom uzerine dusurulur. Boylece bu daldan itibaren olusan
                # mutasyonlar 4 kopyanin BIRINDE kalir (VAF~0.25); ustteki mutasyonlar
                # ikilenmeyle 2 kopyaya cikar (VAF~0.50).
                wfa = os.path.join(a.out, f"{node}_wgd_in", "genome.fa")
                if not (os.path.exists(wfa) and os.path.getsize(wfa) > 1000):
                    if not os.path.exists(pfa):
                        log(f"{node}: ebeveyn genomu yok ({pfa}) -- atlandi")
                        continue
                    log(f"{node}: WGD uygulaniyor (kontig ikilemesi) -> {wfa}")
                    r = subprocess.run([sys.executable,
                                        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     "apply_wgd.py"),
                                        "--in", pfa, "--out", wfa],
                                       capture_output=True, text=True)
                    if r.returncode != 0 or not os.path.exists(wfa):
                        log(f"{node}: WGD BASARISIZ: {r.stderr[-400:]}")
                        continue
                    log(f"{node}: WGD tamam -- {r.stdout.strip()}")
                pfa = wfa
        if not os.path.exists(pfa):
            log(f"{node}: ebeveyn genomu yok ({pfa}) -- atlandi")
            continue
        realize_node(node, pfa, a.out, budgets[node], None, spec, chans, cum,
                     a.seed + (abs(hash(node)) % 100000), log)
    logf.close()


if __name__ == "__main__":
    main()
