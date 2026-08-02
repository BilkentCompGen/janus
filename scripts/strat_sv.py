#!/usr/bin/env python3
"""DEL and DUP, each split into unique and repetitive reference, on two builds side by side.

Same metric as scan_strat.py -- a candidate is attributed to the truth record it shares the most
cores with, precision is "that record is a >=500 bp event of the target type", recall is over the
truth records of that type. The only thing added here is the split: every candidate and every truth
record gets a region class from the copy-number map of ITS OWN build (how many places in the
reference emit that signature), and the table is cut into

    unique     the signature occurs once in the reference
    problem    it occurs 2+ times -- collapsed sequence, where a dosage test is diluted by 1/C

so calls and truth in a stratum are scored against each other and nothing else.

Detectors are the mirror image of one another. DEL: runs of tcov/ncov < LO. DUP: runs of
tcov/ncov > HI with both counts positive. With REP finite, a coverage-derived repeat node may be
traversed without carrying the run; REP is left off by default because measurement says the
transparency merges unrelated runs through the repeat network and costs far more precision than the
bridging buys -- and, decisively, both columns of a comparison must use the same rule.

usage: strat_sv.py <labelA> <graphA> <gtA> <cnA> <labelB> <graphB> <gtB> <cnB> [REP=inf]
"""
import collections
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
from scan_dosage import components, load_graph

BIG = 500
LO = 0.1
HI = 1.5
MINC = 3
REP = float("inf")
XS = (3, 5, 6, 8, 10, 12, 15, 20)


def load_cn(path):
    n = int(np.fromfile(path, dtype=np.uint64, count=1)[0])
    raw = np.fromfile(path, dtype=np.uint8, offset=8, count=n * 12).reshape(n, 12)
    return raw[:, :8].copy().view(np.uint64).ravel(), raw[:, 8:].copy().view(np.uint32).ravel()


def cn_of(sigs, cs, cc):
    q = np.asarray(sigs, dtype=np.uint64)
    i = np.clip(np.searchsorted(cs, q), 0, len(cs) - 1)
    return np.where(cs[i] == q, cc[i], 0)


def stratum(v):
    """v = copy numbers observed at a locus; 'unique' unless the median says 2+ copies."""
    v = v[v > 0]
    return "unique" if not len(v) or np.median(v) < 2 else "problem"


def build(label, gpath, gt, cnpath):
    cs, cc = load_cn(cnpath)
    sig, ncov, tcov, src, dst = load_graph(gpath)
    dose = np.where(ncov > 0, tcov.astype(np.float64) / np.maximum(ncov, 1), np.inf)
    backbone = float(np.median(ncov[(ncov > 0) & (tcov > 0)]))
    repeat = (ncov.astype(np.float64) / backbone) > REP
    masks = {"DEL": (ncov > 0) & (dose < LO),
             "DUP": (ncov > 0) & (tcov > 0) & (dose > HI)}
    print(f"# {label}: K={len(sig):,} backbone={backbone:.0f} repeat={int(repeat.sum()):,} "
          f"depleted={int(masks['DEL'].sum()):,} amplified={int(masks['DUP'].sum()):,} "
          f"refsigs={len(cs):,}", flush=True)

    # truth: every core any record touches, plus the per-record class from its interior
    idx = collections.defaultdict(list)
    recs, truth = {}, collections.defaultdict(lambda: collections.defaultdict(set))
    for line in open(gt):
        r = json.loads(line)
        rid, t, sl = r["id"], r["svtype"], abs(r["svlen"])
        recs[rid] = (t, sl)
        for k in ("over", "over_rc", "inside", "inside_rc"):
            for c in r.get(k) or []:
                idx[c[2]].append(rid)
        if t in ("DEL", "DUP") and sl >= BIG:
            ins = sorted({c[2] for c in (r.get("inside") or [])} |
                         {c[2] for c in (r.get("inside_rc") or [])})
            truth[t][stratum(cn_of(ins, cs, cc)) if ins else "unique"].add(rid)
    keys = np.array(sorted(idx), dtype=np.uint64)

    out = {}
    for t, mask in masks.items():
        raw = components(mask | repeat, src, dst, MINC)
        cands = []
        for c in raw:
            a = np.asarray(c, dtype=np.int64)
            d = a[mask[a]]
            if len(d) >= MINC:
                cands.append(d)
        rows = []
        for c in cands:
            s = sig[c]
            i = np.clip(np.searchsorted(keys, s), 0, len(keys) - 1)
            hits = collections.Counter()
            for k in s[keys[i] == s]:
                hits.update(idx[int(k)])
            rid = hits.most_common(1)[0][0] if hits else None
            good = rid is not None and recs[rid][0] == t and recs[rid][1] >= BIG
            rows.append((len(c), stratum(cn_of(s, cs, cc)), rid if good else None))
        out[t] = rows
        print(f"#   {t}: {len(cands):,} candidates", flush=True)
    return out, truth


def table(t, strat, A, B, la, lb):
    (ra, ta), (rb, tb) = A, B
    print(f"\n### {t} >= {BIG} bp  --  {strat} region"
          f"   (truth: {la} {len(ta[t][strat])}, {lb} {len(tb[t][strat])})")
    print(f"{'X':>4} | {la+' calls':>12}{'prec':>8}{'recall':>9} | "
          f"{lb+' calls':>12}{'prec':>8}{'recall':>9}")
    for X in XS:
        cells = []
        for rows, truth, in ((ra[t], ta), (rb[t], tb)):
            sub = [r for r in rows if r[0] >= X and r[1] == strat]
            good = [r for r in sub if r[2]]
            den = truth[t][strat]
            hit = {r[2] for r in good} & den          # a call may sit on a record of the other stratum
            cells.append(f"{len(sub):>12,}{(len(good)/len(sub) if sub else 0):>8.1%}"
                         f"{(len(hit)/len(den) if den else 0):>9.1%}")
        print(f"{X:>4} | {cells[0]} | {cells[1]}")


def main():
    global REP
    la, ga, gta, cna, lb, gb, gtb, cnb = sys.argv[1:9]
    if len(sys.argv) > 9:
        REP = float(sys.argv[9])
    A = build(la, ga, gta, cna)
    B = build(lb, gb, gtb, cnb)
    for t in ("DEL", "DUP"):
        for s in ("unique", "problem"):
            table(t, s, A, B, la, lb)


if __name__ == "__main__":
    main()
