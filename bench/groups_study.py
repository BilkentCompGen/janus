#!/usr/bin/env python3
import sys, numpy as np
from collections import defaultdict
MM = int(sys.argv[2]) if len(sys.argv) > 2 else 3
S = defaultdict(list)
for l in open(sys.argv[1]):
    s, r, sc, d = l.split()
    S[int(s)].append(int(d))
D = {s: np.sort(np.array(v)) for s, v in S.items()}
n = np.array([len(v) for v in D.values()])
print(f"sites with at least one measured tumor read: {len(D)}; reads per site: median {np.median(n):.0f}, "
      f"mean {n.mean():.1f}, sites with >= {MM}: {(n >= MM).sum()}")
big = {s: v for s, v in D.items() if len(v) >= MM}
med = np.array([np.median(v) for v in big.values()])
spread = np.array([v[-1] - v[0] for v in big.values()])
mad = np.array([np.median(np.abs(v - np.median(v))) for v in big.values()])
print(f"\nanchor distance (median per site, bp): quartiles {np.percentile(med, [25, 50, 75]).round(0)}, 90% {np.percentile(med, 90):.0f}")
print(f"spread max-min per site (bp): quartiles {np.percentile(spread, [25, 50, 75]).round(1)}, 90% {np.percentile(spread, 90):.0f}, 99% {np.percentile(spread, 99):.0f}")
print(f"MAD per site (bp): quartiles {np.percentile(mad, [25, 50, 75]).round(1)}, 90% {np.percentile(mad, 90):.1f}, 99% {np.percentile(mad, 99):.1f}")
rel = mad / np.maximum(med, 1)
print(f"MAD relative to the distance: median {np.median(rel) * 100:.3f} %, 90% {np.percentile(rel, 90) * 100:.3f} %, 99% {np.percentile(rel, 99) * 100:.2f} %")
for lo, hi in ((0, 1000), (1000, 3000), (3000, 6000), (6000, 10000), (10000, 1e9)):
    k = (med >= lo) & (med < hi)
    if k.sum(): print(f"  distance {lo:>5}-{hi if hi < 1e9 else 'inf':>5}: sites {k.sum():6d}, MAD median {np.median(mad[k]):.1f} bp, 90% {np.percentile(mad[k], 90):.1f}, 99% {np.percentile(mad[k], 99):.1f}")
gaps = np.concatenate([np.diff(v) for v in big.values()])
print(f"\ngaps between neighbouring distances within a site: n {len(gaps)}, zero {np.mean(gaps == 0) * 100:.1f} %, <=2 {np.mean(gaps <= 2) * 100:.1f} %, "
      f"<=5 {np.mean(gaps <= 5) * 100:.1f} %, <=10 {np.mean(gaps <= 10) * 100:.1f} %, <=20 {np.mean(gaps <= 20) * 100:.2f} %, >50 {np.mean(gaps > 50) * 100:.2f} %, >200 {np.mean(gaps > 200) * 100:.2f} %")
hist = [(0, 1), (1, 3), (3, 6), (6, 11), (11, 21), (21, 31), (31, 51), (51, 101), (101, 201), (201, 501), (501, 2001), (2001, 10**9)]
print("gap histogram:", ", ".join(f"{a}-{b - 1 if b < 10**9 else 'inf'}: {((gaps >= a) & (gaps < b)).sum()}" for a, b in hist))

def split(v, T):
    g, c = [], 1
    for x in np.diff(v):
        if x <= T: c += 1
        else: g.append(c); c = 1
    g.append(c)
    return g

print(f"\nsites giving >= 2 groups of >= {MM} reads (what adds pieces), by absolute threshold:")
for T in (2, 5, 10, 20, 30, 50, 100, 200, 500, 1000):
    k = sum(1 for v in big.values() if sum(1 for x in split(v, T) if x >= MM) >= 2)
    print(f"  {T:>5} bp: {k:6d} sites ({k / len(big) * 100:.2f} %)")
print("by relative threshold (fraction of the site's median distance):")
for f in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05):
    k = sum(1 for v in big.values() if sum(1 for x in split(v, max(1, f * np.median(v))) if x >= MM) >= 2)
    print(f"  {f * 100:5.1f} %: {k:6d} sites")
sep = []
for v in big.values():
    gs, cur = [], [v[0]]
    for a, b in zip(v[:-1], v[1:]):
        if b - a <= 20: cur.append(b)
        else: gs.append(cur); cur = [b]
    gs.append(cur)
    gs = [g for g in gs if len(g) >= MM]
    if len(gs) >= 2:
        gs.sort(key=len, reverse=True)
        sep.append(abs(np.median(gs[0]) - np.median(gs[1])))
if sep:
    sep = np.array(sep)
    print(f"\nsites split at 20 bp: {len(sep)}; distance between the two largest groups' medians (bp): quartiles {np.percentile(sep, [25, 50, 75]).round(0)}, "
          f"<=50: {np.mean(sep <= 50) * 100:.0f} %, <=100: {np.mean(sep <= 100) * 100:.0f} %, >500: {np.mean(sep > 500) * 100:.0f} %")
