#!/usr/bin/env python3
"""Which deletions does one graph find that the other does not, and why.

Runs the same detector on two (graph, truth) pairs -- normally HPC vs no-HPC -- takes the set of
DEL >= 500 bp each one detects at the loosest threshold, and reports the difference. For every
deletion in the difference it then asks the questions that decide whether a coverage detector can
see it at all:

  present   are its interior cores in the graph? (absent = the parse never produced them)
  depleted  of those present, how many actually satisfy tcov/ncov < 0.1?
  ncov      their median normal depth -- far above backbone means a collapsed repeat, where
            deleting one copy moves the ratio by 1/C and the test cannot fire
  region    the reference copy number the truth record was annotated with

usage: lost_dels.py <graphA.bin> <gtA.jsonl> <labelA> <graphB.bin> <gtB.jsonl> <labelB>
"""
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
from scan_dosage import components, load_graph

BIG = 500
LO = 0.1


def truth_dels(gt):
    """id -> record, keeping the interior cores of every DEL >= BIG."""
    out = {}
    for line in open(gt):
        r = json.loads(line)
        if r["svtype"] != "DEL" or abs(r["svlen"]) < BIG:
            continue
        r["_cores"] = [c[2] for c in (r.get("inside") or [])] + [c[2] for c in (r.get("inside_rc") or [])]
        out[r["id"]] = r
    return out


def detected(gpath, gt):
    """The set of DEL ids some candidate run touches, plus the graph arrays for diagnosis."""
    sig, ncov, tcov, src, dst = load_graph(gpath)
    dose = np.where(ncov > 0, tcov.astype(np.float64) / np.maximum(ncov, 1), np.inf)
    depleted = (ncov > 0) & (dose < LO)
    cands = components(depleted, src, dst, 3)
    recs = truth_dels(gt)
    idx = {}
    for rid, r in recs.items():
        for s in r["_cores"]:
            idx[s] = rid
    keys = np.array(sorted(idx), dtype=np.uint64)
    found = set()
    for c in cands:
        s = sig[np.asarray(c, dtype=np.int64)]
        i = np.clip(np.searchsorted(keys, s), 0, len(keys) - 1)
        found |= {idx[int(k)] for k in s[keys[i] == s]}
    backbone = float(np.median(ncov[(ncov > 0) & (tcov > 0)]))
    return found, recs, (sig, ncov, tcov, depleted, backbone)


def diagnose(rid, recs, g):
    sig, ncov, tcov, depleted, backbone = g
    r = recs[rid]
    q = np.array(sorted(set(r["_cores"])), dtype=np.uint64)
    i = np.clip(np.searchsorted(sig, q), 0, len(sig) - 1)
    hit = sig[i] == q
    j = i[hit]
    if not len(j):
        return dict(cores=len(q), present=0.0, dep=0.0, ncov=0, dose=float("nan"))
    d = np.where(ncov[j] > 0, tcov[j] / np.maximum(ncov[j], 1), np.inf)
    return dict(cores=len(q), present=hit.mean(), dep=depleted[j].mean(),
                ncov=int(np.median(ncov[j])), dose=float(np.median(d)), backbone=backbone)


def main():
    ga, ta, la, gb, tb, lb = sys.argv[1:7]
    fa, ra, aa = detected(ga, ta)
    fb, rb, bb = detected(gb, tb)
    print(f"{la}: {len(fa)}/{len(ra)} detected   |   {lb}: {len(fb)}/{len(rb)} detected\n")

    for tag, ids, recs, g, other in ((f"found by {la}, MISSED by {lb}", fa - fb, ra, aa, bb),
                                     (f"found by {lb}, MISSED by {la}", fb - fa, rb, bb, aa)):
        print(f"=== {tag}: {len(ids)} ===")
        if not ids:
            continue
        print(f"{'id':>10}{'chrom':>7}{'pos':>12}{'svlen':>9}{'region':>15}"
              f"{'cores':>7}{'present':>9}{'depleted':>10}{'ncov':>7}{'dose':>7}"
              f"{'| other: dep':>13}{'ncov':>7}")
        for rid in sorted(ids, key=lambda x: -abs(recs[x]["svlen"])):
            r = recs[rid]
            a = diagnose(rid, recs, g)
            o = diagnose(rid, recs, other) if rid in recs else None
            print(f"{rid:>10}{str(r['chrom']):>7}{r['pos']:>12,}{abs(r['svlen']):>9,}"
                  f"{str(r.get('region', '?')):>15}{a['cores']:>7}{a['present']:>9.1%}"
                  f"{a['dep']:>10.1%}{a['ncov']:>7}{a['dose']:>7.2f}"
                  f"{o['dep']:>13.1%}{o['ncov']:>7}")
        print(f"(backbone depth: {la}={aa[4]:.0f}  {lb}={bb[4]:.0f})\n")


if __name__ == "__main__":
    main()
