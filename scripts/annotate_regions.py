#!/usr/bin/env python3
"""Tag every truth variant with the copy number of the reference sequence it sits in.

tools/ref_copynum counted, for every LCP core of CHM13, how many places in the reference produce that
signature. This joins that count back onto the truth records, because copy number is what decides
whether a coverage-based caller can see a variant at all: deleting one copy of a C-copy sequence moves
its depth by 1/C.

Adds to each record of outputs/svp/gt/*.jsonl:
    cn_inside   median copy number of the cores INSIDE the variant   (0 if it has none)
    cn_over     median copy number of the cores it OVERLAPS
    cn_flank    median copy number of its pre/post anchors
    cn_max      the largest of the three medians -- the number that decides difficulty
    region      unique | low(2-9) | repeat(10-99) | high(100+)

usage: annotate_regions.py <copynum.bin> <gt.jsonl> [out.jsonl]
"""
import collections
import json
import sys

import numpy as np


def load_copynum(path):
    n = int(np.fromfile(path, dtype=np.uint64, count=1)[0])
    raw = np.fromfile(path, dtype=np.uint8, offset=8, count=n * 12).reshape(n, 12)
    sig = raw[:, :8].copy().view(np.uint64).ravel()
    cnt = raw[:, 8:].copy().view(np.uint32).ravel()
    assert np.all(sig[:-1] < sig[1:]), "copynum.bin is not sorted by sig"
    return sig, cnt


def med_cn(sig, cnt, cores):
    if not cores:
        return 0
    q = np.array([c[2] for c in cores], dtype=np.uint64)
    i = np.clip(np.searchsorted(sig, q), 0, len(sig) - 1)
    v = np.where(sig[i] == q, cnt[i], 0)
    v = v[v > 0]
    return int(np.median(v)) if len(v) else 0


def region_of(cn):
    if cn >= 100:
        return "high(100+)"
    if cn >= 10:
        return "repeat(10-99)"
    if cn >= 2:
        return "low(2-9)"
    return "unique"


def main():
    cpath, gt = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else gt
    sig, cnt = load_copynum(cpath)
    print(f"copynum: {len(sig):,} signatures, max copies {cnt.max():,}")

    rows = []
    stat = collections.defaultdict(collections.Counter)
    for line in open(gt):
        r = json.loads(line)
        ins = (r.get("inside") or []) + (r.get("inside_rc") or [])
        ovr = (r.get("over") or []) + (r.get("over_rc") or [])
        fl = r["pre"] + r["post"] + r["pre_rc"] + r["post_rc"]
        r["cn_inside"] = med_cn(sig, cnt, ins)
        r["cn_over"] = med_cn(sig, cnt, ovr)
        r["cn_flank"] = med_cn(sig, cnt, fl)
        r["cn_max"] = max(r["cn_inside"], r["cn_over"], r["cn_flank"])
        r["region"] = region_of(r["cn_max"])
        stat[r["svtype"]][r["region"]] += 1
        rows.append(json.dumps(r))

    with open(out + ".tmp", "w") as f:
        f.write("\n".join(rows) + "\n")
    import os
    os.replace(out + ".tmp", out)

    cols = ["unique", "low(2-9)", "repeat(10-99)", "high(100+)"]
    print(f"\n{'svtype':<12}{'n':>7}" + "".join(f"{c:>16}" for c in cols))
    for t, c in sorted(stat.items(), key=lambda kv: -sum(kv[1].values())):
        n = sum(c.values())
        print(f"{t:<12}{n:>7}" + "".join(f"{c[k]/n:>15.1%} " for k in cols))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
