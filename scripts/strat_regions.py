#!/usr/bin/env python3
"""Stratify the DEL detector by how repetitive the reference is where the deletion sits.

The aggregate precision/recall of a coverage detector is a mixture, and the mixing weight is the
genome itself: in unique sequence a deleted copy takes the depth to zero, and in a C-copy collapse it
moves it by 1/C, which the ratio test cannot see. So the only honest way to compare two builds is
per region class, on each build's OWN reference map -- the classes are not the same sequence twice,
because the parse decides what counts as repetitive.

Region class of a core = how many places in the reference produce its signature:
    unique(1) | low(2-9) | repeat(10-99) | high(100+)
A deletion's class is the median over its interior cores; a call's class is the median over its nodes.

Reports, per build:
  A  the reference map itself -- what fraction of cores and of the genome falls in each class
  B  how the truth deletions distribute over the classes
  C  recall per class at a range of thresholds
  D  precision per class (a call is right if it contains a core from inside a DEL >= 500 bp)

usage: strat_regions.py <label> <graph.bin> <gt.jsonl> <copynum.bin> [<label2> <graph2> <gt2> <cn2> ...]
"""
import collections
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
from scan_dosage import components, load_graph

BIG = 500
LO = 0.1
CLASSES = ["unique", "low(2-9)", "repeat(10-99)", "high(100+)"]


def klass(cn):
    if cn >= 100:
        return "high(100+)"
    if cn >= 10:
        return "repeat(10-99)"
    if cn >= 2:
        return "low(2-9)"
    return "unique"


def load_cn(path):
    n = int(np.fromfile(path, dtype=np.uint64, count=1)[0])
    raw = np.fromfile(path, dtype=np.uint8, offset=8, count=n * 12).reshape(n, 12)
    return raw[:, :8].copy().view(np.uint64).ravel(), raw[:, 8:].copy().view(np.uint32).ravel()


def cn_of(sigs, cs, cc):
    """Reference copy number of each signature; 0 when the signature is not in the reference."""
    q = np.asarray(sigs, dtype=np.uint64)
    i = np.clip(np.searchsorted(cs, q), 0, len(cs) - 1)
    return np.where(cs[i] == q, cc[i], 0)


def run(label, gpath, gt, cnpath):
    cs, cc = load_cn(cnpath)
    sig, ncov, tcov, src, dst = load_graph(gpath)
    dose = np.where(ncov > 0, tcov.astype(np.float64) / np.maximum(ncov, 1), np.inf)
    depleted = (ncov > 0) & (dose < LO)
    cands = components(depleted, src, dst, 3)

    print(f"\n{'='*96}\n### {label}   graph K={len(sig):,}  depleted={int(depleted.sum()):,}  "
          f"candidates={len(cands):,}  reference sigs={len(cs):,}\n{'='*96}")

    # --- A: the reference map ---
    ka = np.where(cc >= 100, 3, np.where(cc >= 10, 2, np.where(cc >= 2, 1, 0)))
    print(f"\nA. reference core signatures by class")
    print(f"{'class':>16}{'signatures':>14}{'share':>9}{'occurrences':>15}{'occ share':>11}")
    tot_occ = cc.sum()
    for i, name in enumerate(CLASSES):
        m = ka == i
        print(f"{name:>16}{int(m.sum()):>14,}{m.mean():>9.2%}{int(cc[m].sum()):>15,}"
              f"{cc[m].sum()/tot_occ:>11.2%}")

    # --- truth ---
    recs = {}
    for line in open(gt):
        r = json.loads(line)
        if r["svtype"] != "DEL" or abs(r["svlen"]) < BIG:
            continue
        cores = [c[2] for c in (r.get("inside") or [])] + [c[2] for c in (r.get("inside_rc") or [])]
        if not cores:
            continue
        v = cn_of(sorted(set(cores)), cs, cc)
        v = v[v > 0]
        recs[r["id"]] = dict(cores=set(cores), cn=int(np.median(v)) if len(v) else 0,
                             svlen=abs(r["svlen"]))
    for r in recs.values():
        r["class"] = klass(r["cn"])

    print(f"\nB. truth DEL >= {BIG}bp by class ({len(recs)} records)")
    print(f"{'class':>16}{'n':>7}{'median size':>14}")
    byclass = collections.defaultdict(list)
    for rid, r in recs.items():
        byclass[r["class"]].append(rid)
    for name in CLASSES:
        ids = byclass.get(name, [])
        if ids:
            print(f"{name:>16}{len(ids):>7}{int(np.median([recs[i]['svlen'] for i in ids])):>14,}")

    # --- candidates: which truth records they contain, and their own class ---
    idx = {}
    for rid, r in recs.items():
        for s in r["cores"]:
            idx[s] = rid
    keys = np.array(sorted(idx), dtype=np.uint64)
    cand_info = []
    for c in cands:
        node = np.asarray(c, dtype=np.int64)
        s = sig[node]
        i = np.clip(np.searchsorted(keys, s), 0, len(keys) - 1)
        hits = {idx[int(k)] for k in s[keys[i] == s]}
        v = cn_of(s, cs, cc)
        v = v[v > 0]
        cand_info.append((len(c), hits, klass(int(np.median(v))) if len(v) else "unique"))

    print(f"\nC. recall per class  (of the deletions in that class, how many a call contains)")
    hdr = "".join(f"{c:>16}" for c in CLASSES)
    print(f"{'X':>4}{hdr}{'TOTAL':>12}")
    for X in (3, 5, 6, 8, 10, 12, 15, 20):
        found = set()
        for ln, hits, _ in cand_info:
            if ln >= X:
                found |= hits
        row = ""
        for name in CLASSES:
            ids = byclass.get(name, [])
            row += f"{(sum(1 for i in ids if i in found)/len(ids)) if ids else 0:>15.1%} " if ids else f"{'-':>16}"
        print(f"{X:>4}{row}{len(found)/max(len(recs),1):>11.1%} ")

    print(f"\nD. precision per class  (of the calls whose nodes are in that class, how many contain a deletion)")
    print(f"{'X':>4}{hdr}{'TOTAL':>12}")
    for X in (3, 5, 6, 8, 10, 12, 15, 20):
        sub = [(ln, h, k) for ln, h, k in cand_info if ln >= X]
        row = ""
        for name in CLASSES:
            s2 = [x for x in sub if x[2] == name]
            row += f"{(sum(1 for x in s2 if x[1])/len(s2)) if s2 else 0:>14.1%}{'':>1}" + f"{len(s2):>1}" if False else (
                f"{(sum(1 for x in s2 if x[1])/len(s2)) if s2 else 0:>15.1%} " if s2 else f"{'-':>16}")
        tot = sum(1 for x in sub if x[1]) / len(sub) if sub else 0
        print(f"{X:>4}{row}{tot:>11.1%} ")

    print(f"\n   call volume per class")
    print(f"{'X':>4}{hdr}{'TOTAL':>12}")
    for X in (3, 6, 10, 15):
        sub = [(ln, h, k) for ln, h, k in cand_info if ln >= X]
        row = "".join(f"{sum(1 for x in sub if x[2]==name):>15,} " for name in CLASSES)
        print(f"{X:>4}{row}{len(sub):>11,} ")


def main():
    a = sys.argv[1:]
    for i in range(0, len(a), 4):
        run(a[i], a[i + 1], a[i + 2], a[i + 3])


if __name__ == "__main__":
    main()
