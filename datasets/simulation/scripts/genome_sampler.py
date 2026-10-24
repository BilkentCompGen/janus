#!/usr/bin/env python3
"""Genom olcekli, imza-kosullu konum ornekleyici (numpy).

Neden ayri bir modul
--------------------
Saf-Python bir trinukleotid indeksi 6.2 Gb'lik diploid bir genomda dugum basina ~1 saat surerdi
(olculdu: ~1.5 M baz/s). Burada dizi bir kez uint8 diziye cevrilir, 32 pirimidin-merkezli baglam
kodu VEKTOREL hesaplanir, ve konumlar REDDETME ORNEKLEMESI ile toplu cekilir. Dugum basina
saniyeler surer.

Sozlesme
--------
* Konumlar 0-tabanlidir ve kontig-icidir.
* Baglam kodu: pirimidin merkezli 3-mer -> 0..31. Merkez purin ise ters tamamlayan alinir.
* N iceren veya kenar 3-mer'ler -> kod -1 (asla secilmez).
"""
import numpy as np

BASES = "ACGT"
_B2I = np.full(256, -1, dtype=np.int8)
for i, b in enumerate(BASES):
    _B2I[ord(b)] = i
    _B2I[ord(b.lower())] = i
_COMP_IDX = np.array([3, 2, 1, 0], dtype=np.int8)          # A<->T, C<->G

# pirimidin merkezli 3-mer listesi (32 tane), COSMIC sirasiyla uyumlu olacak sekilde uretilir
TRIS = [a + c + b for c in "CT" for a in BASES for b in BASES]
TRI2CODE = {t: i for i, t in enumerate(TRIS)}
CODE2TRI = {i: t for t, i in TRI2CODE.items()}


def encode(seq_str):
    """FASTA dizisi -> int8 dizi (A0 C1 G2 T3, diger -1)."""
    arr = np.frombuffer(seq_str.encode("ascii", "ignore"), dtype=np.uint8)
    return _B2I[arr]


def context_codes(idx):
    """int8 baz dizisi -> her konum icin pirimidin-merkezli 3-mer kodu (0..31) veya -1.

    Merkez baz pirimidin (C=1,T=3) ise 3-mer dogrudan; purin (A=0,G=2) ise 3-mer'in
    ters tamamlayani alinir. Bu, COSMIC'in pirimidin-merkezli sozlesmesidir.
    """
    n = idx.size
    code = np.full(n, -1, dtype=np.int8)
    if n < 3:
        return code
    left, mid, right = idx[:-2], idx[1:-1], idx[2:]
    valid = (left >= 0) & (mid >= 0) & (right >= 0)

    is_pyr = (mid == 1) | (mid == 3)                        # C veya T
    # ileri yon: (mid_pyr_index, left, right)
    fwd_c = np.where(mid == 1, 0, 1)                        # C->0, T->1
    fwd = fwd_c * 16 + left * 4 + right
    # ters yon (merkez purin): tamamla ve sagi/solu takas et
    rmid = _COMP_IDX[np.clip(mid, 0, 3)]
    rl = _COMP_IDX[np.clip(right, 0, 3)]
    rr = _COMP_IDX[np.clip(left, 0, 3)]
    rev_c = np.where(rmid == 1, 0, 1)
    rev = rev_c * 16 + rl * 4 + rr

    res = np.where(is_pyr, fwd, rev).astype(np.int8)
    res[~valid] = -1
    code[1:-1] = res
    return code


class GenomeContexts:
    """Bir genomun (kontig -> dizi) baglam kodlarini tutar ve konum ornekler."""

    def __init__(self, contigs, rng):
        """contigs: {name: str}. Bellek: her kontig icin 1 byte/baz kod dizisi."""
        self.names = list(contigs)
        self.seqs = contigs
        self.codes = {}
        self.lengths = {}
        for n, s in contigs.items():
            idx = encode(s)
            self.codes[n] = context_codes(idx)
            self.lengths[n] = idx.size
        self.total_len = sum(self.lengths.values())
        self._names_arr = np.array(self.names, dtype=object)
        self._w = np.array([self.lengths[n] for n in self.names], dtype=np.float64)
        self._w /= self._w.sum()
        self.rng = rng

    def genome_context_freq(self, sample=5_000_000):
        """Baglam frekanslarini (32 kod) tum genomdan sayar."""
        tot = np.zeros(32, dtype=np.int64)
        for n in self.names:
            c = self.codes[n]
            tot += np.bincount(c[c >= 0], minlength=32)
        return tot / tot.sum()

    def sample_by_context(self, need, forbidden, batch_mult=40, max_rounds=40):
        """need: {context_code: adet}. -> {code: [(contig, pos0), ...]}

        Reddetme ornekleme: genom boyunca duzgun konum cek, kodunu oku, ihtiyaci olan kodlara
        dagit. Nadir baglamlar icin turlar tekrarlanir. `forbidden`: {contig: set(pos)} --
        sonsuz-siteler kurali (ayni konum iki kez mutasyona ugramaz).
        """
        out = {c: [] for c in need}
        remaining = dict(need)
        rng = self.rng
        for _ in range(max_rounds):
            tot_need = sum(v for v in remaining.values() if v > 0)
            if tot_need == 0:
                break
            k = max(10_000, tot_need * batch_mult)
            # kontig sec (uzunlukla agirlikli), sonra kontig icinde konum
            ci = rng.choice(len(self.names), size=k, p=self._w)

            # ---- TUM kontiglerden aday topla, SONRA kuresel olarak karistir ----
            # HATA (2026-08-18, MRCA dugumunde yakalandi): onceki surum kontigleri
            # np.unique(ci) sirasiyla dolasip her kontigte "kalan ihtiyaci" doyuruyordu.
            # Ilk islenen kontig (chr1_hapA) parti buyuk oldugu icin ihtiyacin neredeyse
            # tamamini tek basina karsiliyordu: 5,644 mutasyonun 3,842'si chr1_hapA'ya,
            # 630'u chr1_hapB'ye dustu ve sayilar kontig indeksiyle monoton azaldi.
            # Beklenen ~230 idi. Duzeltme: adaylari once TOPLA, kuresel karistir, sonra ata --
            # boylece atama kontig islenme sirasindan bagimsiz ve uzunlukla orantili olur.
            cand = {}
            for cidx in np.unique(ci):
                name = self.names[cidx]
                m = int((ci == cidx).sum())
                pos = rng.integers(1, self.lengths[name] - 1, size=m)
                codes = self.codes[name][pos]
                for code in remaining:
                    if remaining[code] <= 0:
                        continue
                    sel = pos[codes == code]
                    if sel.size:
                        cand.setdefault(code, []).append((name, sel))
            for code, chunks in cand.items():
                if remaining[code] <= 0:
                    continue
                flat = [(nm, int(p)) for nm, arr in chunks for p in arr]
                rng.shuffle(flat)
                for nm, p in flat:
                    if remaining[code] <= 0:
                        break
                    forb = forbidden.setdefault(nm, set())
                    if p in forb:
                        continue
                    forb.add(p)
                    out[code].append((nm, p))
                    remaining[code] -= 1
        return out, remaining

    def sample_uniform(self, k, forbidden, min_gap=1):
        """Baglamdan bagimsiz, duzgun k konum (indel/SV yerlesimi icin)."""
        res = []
        rng = self.rng
        tries = 0
        while len(res) < k and tries < k * 200 + 1000:
            tries += 1
            cidx = int(rng.choice(len(self.names), p=self._w))
            name = self.names[cidx]
            p = int(rng.integers(1000, max(1001, self.lengths[name] - 1000)))
            if self.codes[name][p] < 0:          # N bolgesi
                continue
            forb = forbidden.setdefault(name, set())
            if any((p + d) in forb for d in range(-min_gap, min_gap + 1)):
                continue
            forb.add(p)
            res.append((name, p))
        return res
