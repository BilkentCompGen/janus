#!/usr/bin/env python3
"""Join the ground-truth cores (outputs/svp/gt/*.jsonl) to a v3 graph, by signature.

The graph stores no sequence and no positions -- a node is a 64-bit sig plus a normal and a tumour
count. gt_cores gives every truth variant as sigs. So the join needs nothing else, and it answers the
question a caller has to answer first: for each variant, are its ANCHOR cores in the graph (and shared
between the two samples, as normal sequence should be), and are its JUNCTION cores in the graph as
tumour-only?

usage: gt_vs_graph.py <graph_v3.bin> <gt.jsonl> [max_records]
"""
import collections
import json
import sys

import numpy as np

COLOUR = {0: "absent", 1: "normal-only", 2: "tumour-only", 3: "shared"}


def load_graph(path):
    K = int(np.fromfile(path, dtype=np.uint64, count=1)[0])
    sig = np.fromfile(path, dtype=np.uint64, count=K, offset=8)
    ncov = np.fromfile(path, dtype=np.uint32, count=K, offset=8 + K * 8)
    tcov = np.fromfile(path, dtype=np.uint32, count=K, offset=8 + K * 12)
    assert np.all(sig[:-1] < sig[1:]), "sig[] is not sorted+unique -- searchsorted would be wrong"
    return sig, ncov, tcov


def lookup(sig, ncov, tcov, q):
    """-> (colour, ncov, tcov) arrays for query sigs; colour 0 = not a node at all."""
    q = np.asarray(q, dtype=np.uint64)
    i = np.searchsorted(sig, q)
    i_clip = np.clip(i, 0, len(sig) - 1)
    hit = sig[i_clip] == q
    n = np.where(hit, ncov[i_clip], 0)
    t = np.where(hit, tcov[i_clip], 0)
    col = np.where(~hit, 0, np.where((n > 0) & (t > 0), 3, np.where(n > 0, 1, 2)))
    return col, n, t


def sweep(sig, ncov, tcov, recs):
    """Where does a ncov/tcov threshold actually separate anchors from junctions?

    An anchor is normal sequence: it should look balanced (n ~ t). A junction is tumour-private: it
    should look extreme (n << t). A threshold is only meaningful if those two populations do not
    overlap, so print both distributions before picking one, then sweep.
    """
    anc_n, anc_t, jun_n, jun_t = [], [], [], []
    for r in recs:
        a = [c[2] for k in ("pre", "post", "pre_rc", "post_rc") for c in r[k]]
        j = [c[2] for k in ("alt", "alt_rc") if r.get(k) for c in r[k] if c[3] == 1]
        for q, dn, dt in ((a, anc_n, anc_t), (j, jun_n, jun_t)):
            if not q:
                continue
            _, n, t = lookup(sig, ncov, tcov, q)
            dn.append(n)
            dt.append(t)
    anc_n, anc_t = np.concatenate(anc_n), np.concatenate(anc_t)
    jun_n, jun_t = np.concatenate(jun_n), np.concatenate(jun_t)

    print(f"\n=== coverage distributions (anchors n={len(anc_n):,}, junctions n={len(jun_n):,}) ===")
    print(f"{'population':<12}{'':>6}" + "".join(f"{p:>8}%" for p in (1, 5, 25, 50, 75, 95, 99)))
    for name, a, b in (("anchor", anc_n, anc_t), ("junction", jun_n, jun_t)):
        for lbl, v in ((f"{name} ncov", a), (f"{name} tcov", b)):
            qs = np.percentile(v, [1, 5, 25, 50, 75, 95, 99])
            print(f"{lbl:<18}" + "".join(f"{q:>9.0f}" for q in qs))

    # ratio r = ncov/tcov, with the tcov==0 case called out separately (it is not a ratio at all)
    print(f"\n=== ncov/tcov ===")
    for name, a, b in (("anchor", anc_n, anc_t), ("junction", jun_n, jun_t)):
        z = int((b == 0).sum())
        m = b > 0
        r = a[m] / b[m]
        qs = np.percentile(r, [1, 5, 25, 50, 75, 95, 99])
        print(f"{name:<10} tcov==0: {z/len(b):>7.2%}   ratio pct: " + " ".join(f"{q:.2f}" for q in qs))

    print(f"\n=== sweep: 'tumour-private' = ncov/tcov < T (tcov==0 excluded as no-ratio) ===")
    print(f"{'T':>6}{'junction kept':>16}{'anchor leaked':>16}")
    for T in (0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
        jm = jun_t > 0
        am = anc_t > 0
        jk = float(((jun_n[jm] / jun_t[jm]) < T).mean())
        al = float(((anc_n[am] / anc_t[am]) < T).mean())
        print(f"{T:>6}{jk:>16.2%}{al:>16.4%}")


def main():
    gpath, gt = sys.argv[1], sys.argv[2]
    cap = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    sig, ncov, tcov = load_graph(gpath)
    print(f"graph {gpath}: K={len(sig):,}")

    recs = []
    for line in open(gt):
        recs.append(json.loads(line))
        if cap and len(recs) >= cap:
            break

    # anchors = pre+post (both families); junctions = the novel cores of alt (both families)
    stat = collections.defaultdict(lambda: collections.Counter())
    for r in recs:
        t = r["svtype"]
        anc = [c[2] for k in ("pre", "post", "pre_rc", "post_rc") for c in r[k]]
        jun = [c[2] for k in ("alt", "alt_rc") if r.get(k) for c in r[k] if c[3] == 1]
        for name, q in (("anchor", anc), ("junction", jun)):
            if not q:
                continue
            col, _, _ = lookup(sig, ncov, tcov, q)
            s = stat[(t, name)]
            s["n"] += len(col)
            for c in (0, 1, 2, 3):
                s[COLOUR[c]] += int((col == c).sum())
        stat[(t, "recs")]["n"] += 1

    types = sorted({k[0] for k in stat}, key=lambda t: -stat[(t, "recs")]["n"])
    for what in ("anchor", "junction"):
        print(f"\n=== {what} cores ===")
        print(f"{'type':<12}{'records':>8}{'cores':>10}{'absent':>9}{'normal':>9}{'tumour':>9}{'shared':>9}")
        for t in types:
            s = stat.get((t, what))
            if not s or not s["n"]:
                continue
            n = s["n"]
            print(f"{t:<12}{stat[(t,'recs')]['n']:>8}{n:>10}"
                  f"{s['absent']/n:>9.2%}{s['normal-only']/n:>9.2%}{s['tumour-only']/n:>9.2%}{s['shared']/n:>9.2%}")

    sweep(sig, ncov, tcov, recs)


if __name__ == "__main__":
    main()
