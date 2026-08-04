#!/usr/bin/env python3
"""Junction-level scoring against a breakend-pair truth set (COLO829 / Valle-Inclan).

That truth set writes every somatic SV as a pair of BND records joined by MATEID, which is the right
representation -- an SV is a novel adjacency -- but it defeats any comparator that matches on SVTYPE,
because every record types as BND while callers emit DEL/DUP/INV/INS. Here type is not asked for at
all: a call matches a truth junction when BOTH of its breakpoints land within TOL bp of the truth's
two breakpoints, on the right contigs, in either pairing order.

usage: colo_bp_eval.py <truth.vcf> <caller1.vcf> [caller2.vcf ...] [tolerance]
"""
import collections
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from severus_bp_eval import load_calls, match

BNDALT = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def load_truth_bnd(path):
    """Breakend pairs -> one junction each: (chrom1, pos1, chrom2, pos2, id, svtype)."""
    rec = {}
    for line in open(path):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        rec[f[2]] = f
    out, seen = [], set()
    for i, r in rec.items():
        m = re.search(r"MATEID=([^;]*)", r[7])
        if not m or m.group(1) not in rec:
            continue
        k = tuple(sorted((i, m.group(1))))
        if k in seen:
            continue
        seen.add(k)
        a, b = rec[k[0]], rec[k[1]]
        t = re.search(r"SVTYPE=([^;]*)", a[7])
        out.append((a[0], int(a[1]), b[0], int(b[1]), k[0],
                    t.group(1) if t else "?", a[0] != b[0]))
    return out


def main():
    args = sys.argv[1:]
    tol = 500
    if args[-1].isdigit():
        tol = int(args[-1]); args = args[:-1]
    truth = load_truth_bnd(args[0])
    print(f"truth: {len(truth)} junction (breakend cifti), tolerans +/-{tol} bp\n")
    print(f"  {'arm':<28s} {'TP':>4s} {'FN':>4s} {'FP':>5s} {'precision':>10s} {'recall':>8s} {'F1':>7s}")
    detail = {}
    for vcf in args[1:]:
        calls, cins = load_calls(vcf)
        hitT, hitC = set(), set()
        for i, a in enumerate(truth):
            for j, c in enumerate(calls):
                if match(c, a[:4] + (a[4],), tol):
                    hitT.add(i); hitC.add(j)
        tp, fn, fp = len(hitT), len(truth) - len(hitT), len(calls) - len(hitC)
        p = tp / (tp + fp) if tp + fp else 0.0
        r = tp / len(truth) if truth else 0.0
        f1 = 2 * p * r / (p + r) if p + r else 0.0
        name = os.path.basename(os.path.dirname(os.path.dirname(vcf)))
        print(f"  {name:<28s} {tp:>4d} {fn:>4d} {fp:>5d} {p:>10.3f} {r:>8.3f} {f1:>7.3f}")
        detail[name] = (hitT, len(calls), len(cins))

    print("\n  bulunma, truth tipine gore:")
    by = collections.defaultdict(list)
    for i, a in enumerate(truth):
        by[a[5] + (" (kromozomlar arasi)" if a[6] else "")].append(i)
    hdr = "  ".join(f"{n:>14s}" for n in detail)
    print(f"  {'tip':<26s} {'n':>4s}  {hdr}")
    for t in sorted(by):
        idx = by[t]
        cells = "  ".join(f"{sum(1 for i in idx if i in detail[n][0]):>6d} "
                          f"{100*sum(1 for i in idx if i in detail[n][0])/len(idx):>6.0f}%" for n in detail)
        print(f"  {t:<26s} {len(idx):>4d}  {cells}")
    for n, (h, nc, ni) in detail.items():
        print(f"  {n}: {nc} junction cagrisi + {ni} insertion cagrisi")


if __name__ == "__main__":
    main()
