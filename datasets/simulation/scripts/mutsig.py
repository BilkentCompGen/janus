#!/usr/bin/env python3
"""Mutasyonel imza tabanli somatik SNV/indel ornekleyici.

Ne yapar
--------
COSMIC v3.4 SBS (96 kanal) ve ID (83 kanal) matrislerinden, meme-HRD karisimina karsilik gelen
HEDEF spektrumu kurar; sonra verilen bir dizi uzerinde o spektruma uyan gercek konumlar secer.

Neden boyle
-----------
Somatik mutasyonlari genoma DUZGUN dagitmak (uniform rastgele) gercek disidir: gercek tumorlerde
mutasyonlarin trinukleotid baglami imzalarla belirlenir (ornegin SBS1 neredeyse yalnizca CpG'de
C>T yapar). Duzgun dagitim, SigProfiler ile geri ayristirildiginda hicbir COSMIC imzasina
benzemez ve "spektrum gercek disi" diye geri doner. Bu yuzden:

  1. hedef 96-kanal olasilik vektoru = Sum_i w_i * SBS_i   (w = spec/somatic.yaml sbs_mix)
  2. dizideki 32 pirimidin-merkezli trinukleotid baglaminin GERCEK frekanslari sayilir
  3. konum secimi baglam-kosullu yapilir: once hedeften bir kanal cekilir, sonra o baglama
     sahip konumlar arasindan duzgun secim yapilir

Boylece uretilen mutasyon kumesi, kendi 96-kanal profili cikarildiginda hedefe kosinus ~1.0
ile oturur -- ve bu, uretecin kendi raporu degil, bagimsiz olarak olculebilir bir sey.
"""
import bisect
import random
from collections import defaultdict

COMP = str.maketrans("ACGT", "TGCA")
PYR = ("C", "T")


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_signature_matrix(path):
    """COSMIC TSV -> {sig_name: {channel: prob}}, kanal adlari dosyadaki gibi."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        sigs = header[1:]
        out = {s: {} for s in sigs}
        channels = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            ch = f[0]
            channels.append(ch)
            for s, v in zip(sigs, f[1:]):
                out[s][ch] = float(v)
    return out, channels


def mix_spectrum(matrix, weights):
    """{sig: w} -> normalize edilmis {channel: prob}. Bilinmeyen imza adinda hata verir."""
    missing = [s for s in weights if s not in matrix]
    if missing:
        raise KeyError(f"imza matrisinde yok: {missing}")
    tot_w = sum(weights.values())
    spec = defaultdict(float)
    for s, w in weights.items():
        for ch, p in matrix[s].items():
            spec[ch] += (w / tot_w) * p
    z = sum(spec.values())
    return {k: v / z for k, v in spec.items()}


def parse_sbs_channel(ch):
    """'A[C>A]G' -> (5' baz, ref, alt, 3' baz)"""
    left, rest = ch.split("[")
    mut, right = rest.split("]")
    ref, alt = mut.split(">")
    return left, ref, alt, right


class ContextIndex:
    """Bir dizideki 32 pirimidin-merkezli trinukleotid baglaminin konumlarini indeksler.

    Bellek: konumlari 1 Mb'lik kovalarda SAYAR, konum listesini tutmaz. Ornekleme sirasinda
    once kova secilir (sayilarla agirlikli), sonra kova taranip o baglamdaki konumlar bulunur.
    3 Gb'lik bir genom icin liste tutmak ~24 GB olurdu; sayaclar ~100 kB.
    """

    BIN = 1 << 20

    def __init__(self, seq, bin_size=None):
        self.seq = seq
        self.BIN = bin_size or ContextIndex.BIN
        nbins = (len(seq) + self.BIN - 1) // self.BIN
        self.counts = defaultdict(lambda: [0] * nbins)
        self.total = defaultdict(int)
        self.nbins = nbins
        for i in range(1, len(seq) - 1):
            c = seq[i]
            if c == "C" or c == "T":
                tri = seq[i - 1:i + 2]
            elif c == "G" or c == "A":
                tri = revcomp(seq[i - 1:i + 2])
            else:
                continue
            if len(tri) != 3 or any(b not in "ACGT" for b in tri):
                continue
            self.counts[tri][i // self.BIN] += 1
            self.total[tri] += 1
        self._cum = {}
        for tri, arr in self.counts.items():
            cum, s = [], 0
            for v in arr:
                s += v
                cum.append(s)
            self._cum[tri] = cum

    def context_frequencies(self):
        z = sum(self.total.values())
        return {k: v / z for k, v in self.total.items()} if z else {}

    def sample_position(self, tri, rng, used):
        """Verilen trinukleotid baglamina sahip bir konum (0-tabanli merkez) dondur.
        `used` kumesindekiler atlanir. Bulunamazsa None."""
        cum = self._cum.get(tri)
        if not cum or cum[-1] == 0:
            return None
        for _ in range(60):
            r = rng.random() * cum[-1]
            b = bisect.bisect_right(cum, r)
            if b >= self.nbins:
                b = self.nbins - 1
            lo, hi = b * self.BIN, min(len(self.seq) - 1, (b + 1) * self.BIN)
            hits = []
            s = self.seq
            for i in range(max(1, lo), hi):
                c = s[i]
                if c == "C" or c == "T":
                    t = s[i - 1:i + 2]
                elif c == "G" or c == "A":
                    t = revcomp(s[i - 1:i + 2])
                else:
                    continue
                if t == tri and i not in used:
                    hits.append(i)
            if hits:
                return rng.choice(hits)
        return None


def build_sbs_sampler(spectrum, ctx_index):
    """Hedef 96-kanal spektrumu + gercek baglam frekanslari -> ornekleyici.

    ONEMLI: hedef spektrum, mutasyonlarin kanal dagilimidir. Bir kanali secmek demek,
    o kanalin trinukleotid baglamina sahip bir konum secmek demektir. Baglamin genomdaki
    bolluguna gore DUZELTME YAPILMAZ -- cunku COSMIC profilleri zaten gozlenen mutasyon
    dagilimidir, genom kompozisyonuna gore normalize edilmemistir. (COSMIC dokumantasyonu
    bunlari 'genome-normalized' degil, gozlenen sayimlar olarak verir.)
    """
    chans, probs = [], []
    for ch, p in spectrum.items():
        if p <= 0:
            continue
        left, ref, alt, right = parse_sbs_channel(ch)
        tri = left + ref + right
        if tri not in ctx_index.total or ctx_index.total[tri] == 0:
            continue
        chans.append((ch, tri, ref, alt))
        probs.append(p)
    z = sum(probs)
    cum, s = [], 0.0
    for p in probs:
        s += p / z
        cum.append(s)
    return chans, cum


def sample_snvs(n, chans, cum, ctx_index, seq, rng, used):
    """-> [(pos0, ref, alt, channel)] ; pos0 = 0-tabanli merkez konumu."""
    out = []
    tries = 0
    while len(out) < n and tries < n * 50:
        tries += 1
        r = rng.random()
        i = bisect.bisect_right(cum, r)
        if i >= len(chans):
            i = len(chans) - 1
        ch, tri, ref, alt = chans[i]
        pos = ctx_index.sample_position(tri, rng, used)
        if pos is None:
            continue
        base = seq[pos]
        # kanal pirimidin merkezlidir; genomda purin ise ters ipliktedir
        if base == ref:
            r_out, a_out = ref, alt
        elif base == ref.translate(COMP):
            r_out, a_out = base, alt.translate(COMP)
        else:
            continue
        used.add(pos)
        out.append((pos, r_out, a_out, ch))
    return out


# ------------------------------------------------------------------ indel
def parse_id_channel(ch):
    """COSMIC ID kanal adi, ornegin '1:Del:C:0' / '5:Ins:R:2' / '2:Del:M:1'
    -> (uzunluk_kovasi, Del|Ins, alt_tip, tekrar/mikrohomoloji_kovasi)"""
    p = ch.split(":")
    return p[0], p[1], p[2], p[3]


def id_channel_to_length(ch, rng):
    """Kanal adindan somut bir indel uzunlugu uret.
    COSMIC ID kanallarinda uzunluk kovasi '1'..'5+' seklindedir; '5+' icin gercekci bir kuyruk
    gerekir. Kuyruk ussel; ortalama 8 bp (docs/03: ID6 '>=5 bp delesyon' ile tutarli)."""
    size = ch.split(":")[0]
    if size == "5+":
        return 5 + int(rng.expovariate(1 / 3.0))
    try:
        return int(size)
    except ValueError:
        return 1
