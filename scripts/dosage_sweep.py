#!/usr/bin/env python3
"""Can DEL and DUP be called from COVERAGE alone?

A deletion's interior should lose its tumour reads; a duplication's should carry twice as many. Both
claims are about the cores INSIDE the variant, measured against the local backbone (its own anchors) --
never against a global depth, because with the blacklist gone a repeat core can sit at 39,000x.

For every truth record we compute
    dosage = median(tcov inside) / median(ncov inside)      the tumour/normal ratio of the segment
    local  = median(ncov inside) / median(ncov anchors)      a sanity check that the segment is normal-depth
and then sweep a threshold to see whether DEL, DUP and everything else separate.

usage: dosage_sweep.py <graph_v3.bin> <gt.jsonl>
"""
import collections
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
from gt_vs_graph import load_graph, lookup


def med(x):
    return float(np.median(x)) if len(x) else float("nan")


def main():
    gpath, gt = sys.argv[1], sys.argv[2]
    sig, ncov, tcov = load_graph(gpath)
    print(f"graph {gpath}: K={len(sig):,}")

    rows = collections.defaultdict(list)
    for line in open(gt):
        r = json.loads(line)
        ins = [c[2] for k in ("inside", "inside_rc") for c in r.get(k) or []]
        anc = [c[2] for k in ("pre", "post", "pre_rc", "post_rc") for c in r[k]]
        if len(ins) < 3 or not anc:
            continue
        _, i_n, i_t = lookup(sig, ncov, tcov, ins)
        _, a_n, _ = lookup(sig, ncov, tcov, anc)
        base = med(a_n)
        if base <= 0:
            continue
        rows[r["svtype"]].append((med(i_t) / max(med(i_n), 1e-9), med(i_n) / base, len(ins), abs(r["svlen"])))

    print(f"\n=== dosage = median(tcov inside) / median(ncov inside) ===")
    print(f"{'type':<12}{'n':>7}{'5%':>8}{'25%':>8}{'50%':>8}{'75%':>8}{'95%':>8}   {'local depth 50%':>16}")
    order = sorted(rows, key=lambda t: -len(rows[t]))
    for t in order:
        d = np.array([x[0] for x in rows[t]])
        l = np.array([x[1] for x in rows[t]])
        q = np.percentile(d, [5, 25, 50, 75, 95])
        print(f"{t:<12}{len(d):>7}" + "".join(f"{v:>8.2f}" for v in q) + f"   {np.median(l):>16.2f}")

    print(f"\n=== sweep: DEL = dosage < LO, DUP = dosage > HI ===")
    print(f"{'LO':>6}{'DEL kept':>10}{'others<LO':>11}     {'HI':>6}{'DUP kept':>10}{'others>HI':>11}")
    dels = np.array([x[0] for x in rows.get("DEL", [])])
    dups = np.array([x[0] for x in rows.get("DUP", [])])
    other = np.array([x[0] for t in rows if t not in ("DEL", "DUP") for x in rows[t]])
    for lo, hi in ((0.05, 1.3), (0.1, 1.4), (0.2, 1.5), (0.3, 1.6), (0.5, 1.8), (0.7, 2.0)):
        dk = float((dels < lo).mean()) if len(dels) else 0
        ol = float((other < lo).mean()) if len(other) else 0
        uk = float((dups > hi).mean()) if len(dups) else 0
        oh = float((other > hi).mean()) if len(other) else 0
        print(f"{lo:>6}{dk:>10.1%}{ol:>11.3%}     {hi:>6}{uk:>10.1%}{oh:>11.3%}")

    print(f"\n=== dosage vs size (DEL / DUP) ===")
    for t in ("DEL", "DUP"):
        if t not in rows:
            continue
        arr = sorted(rows[t], key=lambda x: x[3])
        for lo, hi in ((0, 100), (100, 500), (500, 2000), (2000, 10000), (10000, 10**9)):
            sub = [x[0] for x in arr if lo <= x[3] < hi]
            if sub:
                print(f"  {t} {lo:>7}-{hi:<10} n={len(sub):>5}  dosage median {np.median(sub):>6.2f}")


if __name__ == "__main__":
    main()
