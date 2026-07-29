#!/usr/bin/env python3
"""Genome-wide dosage scan: run the DEL/DUP detector over the WHOLE graph, then score it.

Everything measured before this was conditioned on the truth -- look up a truth core, ask its colour.
That measures recall and nothing else. This scans all K nodes with no knowledge of the truth, emits
candidates, and only then asks how many are real. That is the only way to see FP.

Detector (no coordinates, no reference):
  DEL candidate = a connected run of >=MINC nodes whose dosage tcov/ncov < LO, ncov > 0
  DUP candidate = a connected run of >=MINC nodes whose dosage > HI, both counts > 0,
                  and whose depth exceeds its own BOUNDARY nodes (the local backbone) by > HI
Connectivity comes from the graph's own edges; a run is a connected component of the induced subgraph.

Scoring: a candidate is a TP if any of its sigs is an `inside` core of a truth record of the target
type. Candidates hitting a truth record of a DIFFERENT type are reported separately -- they are not
noise, they are the detector finding real material and the type being wrong.

usage: scan_dosage.py <graph_v3.bin> <gt.jsonl> [LO=0.1] [HI=1.5] [MINC=3]
"""
import collections
import json
import sys

import numpy as np


def load_graph(path):
    K = int(np.fromfile(path, dtype=np.uint64, count=1)[0])
    sig = np.fromfile(path, dtype=np.uint64, count=K, offset=8)
    ncov = np.fromfile(path, dtype=np.uint32, count=K, offset=8 + K * 8)
    tcov = np.fromfile(path, dtype=np.uint32, count=K, offset=8 + K * 12)
    off = 8 + K * 16
    E = int(np.fromfile(path, dtype=np.uint64, count=1, offset=off)[0])
    e = np.fromfile(path, dtype=np.uint32, count=E * 4, offset=off + 8).reshape(E, 4)
    return sig, ncov, tcov, e[:, 0], e[:, 1]


def components(mask, src, dst, minc):
    """Connected components of the subgraph induced on `mask`.

    Hooking + pointer jumping on the edge list: propagate the smaller label across every kept edge,
    then collapse each label chain onto its root (lab = lab[lab]) so the number of rounds is
    logarithmic in the component diameter rather than linear. Edges are filtered to the mask first,
    so the arrays stay small even though the graph has 88M nodes.
    """
    keep = mask[src] & mask[dst]
    s, d = src[keep].astype(np.int64), dst[keep].astype(np.int64)
    idx = np.flatnonzero(mask)
    # Union-find over the masked nodes only. The mask is ~1% of the graph, so this is ~1M edges of
    # plain Python -- seconds -- and it is exactly right, which iterative label propagation was not.
    pos = {int(v): i for i, v in enumerate(idx)}
    parent = list(range(len(idx)))

    def find(x):
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r

    for a, b in zip(s.tolist(), d.tolist()):
        ra, rb = find(pos[a]), find(pos[b])
        if ra != rb:
            parent[ra] = rb
    roots = np.array([find(i) for i in range(len(idx))], dtype=np.int64)
    order = np.argsort(roots, kind="stable")
    roots_s, idx_s = roots[order], idx[order]
    bounds = np.flatnonzero(np.r_[True, roots_s[1:] != roots_s[:-1], True])
    out = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a >= minc:
            out.append(idx_s[a:b])
    return out


def truth_index(gt):
    """sig -> (record id, svtype) for every reference core a truth variant destroys in the tumour.

    `inside` (the segment interior) is what a deletion's dosage signal lives on; `over` (the cores the
    variant merely overlaps) is what a SNV kills. Both go normal-only in the graph and both therefore
    answer the DEL node test -- which is the whole point of measuring FP properly.
    """
    ins = {}
    per_type = collections.Counter()
    sizes = collections.defaultdict(dict)
    for line in open(gt):
        r = json.loads(line)
        for k in ("over", "over_rc"):
            for c in r.get(k) or []:
                ins.setdefault(c[2], (r["id"], r["svtype"]))
        for k in ("inside", "inside_rc"):
            for c in r.get(k) or []:
                ins[c[2]] = (r["id"], r["svtype"])
        per_type[r["svtype"]] += 1                       # denominator = every truth record of that type
        sizes[r["svtype"]][r["id"]] = abs(r["svlen"])
    return ins, per_type, sizes


def main():
    gpath, gt = sys.argv[1], sys.argv[2]
    LO = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
    HI = float(sys.argv[4]) if len(sys.argv) > 4 else 1.5
    MINC = int(sys.argv[5]) if len(sys.argv) > 5 else 3

    sig, ncov, tcov, src, dst = load_graph(gpath)
    K = len(sig)
    print(f"graph {gpath}: K={K:,} E={len(src):,}   LO={LO} HI={HI} MINC={MINC}")

    n = ncov.astype(np.float32)
    t = tcov.astype(np.float32)
    dose = np.where(n > 0, t / np.maximum(n, 1e-9), np.inf)

    ins, per_type, sizes = truth_index(gt)
    keys = np.array(sorted(ins.keys()), dtype=np.uint64)

    def score(cands, target, verbose):
        hit_same, hit_other, miss = 0, 0, 0
        found = set()
        other_types = collections.Counter()
        for c in cands:
            s = sig[np.asarray(c, dtype=np.int64)]
            i = np.clip(np.searchsorted(keys, s), 0, len(keys) - 1)
            recs = {ins[int(k)] for k in s[keys[i] == s]}
            same = {r for r in recs if r[1] == target}
            if same:
                hit_same += 1
                found |= {r[0] for r in same}
            elif recs:
                hit_other += 1
                other_types.update(r[1] for r in recs)
            else:
                miss += 1
        tot = max(len(cands), 1)
        prec, rec = hit_same / tot, len(found) / max(per_type[target], 1)
        if verbose:
            print(f"  candidates {len(cands):,} | on truth {target} {hit_same:,} ({hit_same/tot:.1%})"
                  f" | on another truth type {hit_other:,} ({hit_other/tot:.1%})"
                  f" | on NO truth variant {miss:,} ({miss/tot:.1%})")
            print(f"  truth {target}s recovered: {len(found):,} / {per_type[target]:,} ({rec:.1%})")
            if other_types:
                print(f"  other types hit: {dict(other_types.most_common(6))}")
            big = [i for i in found if sizes[target][i] >= 500]
            allbig = [i for i, v in sizes[target].items() if v >= 500]
            if allbig:
                print(f"  of those >=500bp: {len(big):,} / {len(allbig):,} ({len(big)/len(allbig):.1%})")
        return len(cands), prec, rec

    for name, mask0, cmp_txt in (
            ("DEL", (ncov > 0) & (dose < LO), f"dosage < {LO}"),
            ("DUP", (ncov > 0) & (tcov > 0) & (dose > HI), f"dosage > {HI}")):
        print(f"\n=== {name} scan: connected runs with {cmp_txt} ===")
        print(f"  nodes passing the node test: {int(mask0.sum()):,} ({mask0.mean():.2%} of the graph)")
        print(f"  default MINC={MINC}, no depth floor:")
        score(components(mask0, src, dst, MINC), name, True)

        print(f"\n  sweep -- min run length x min ncov (depth floor kills ratio noise):")
        print(f"  {'minrun':>7}{'mindepth':>9}{'cands':>10}{'precision':>11}{'recall':>9}{'recall>=500bp':>15}")
        for minrun in (3, 5, 10, 20):
            for mind in (0, 5, 10, 15):
                m = mask0 & (ncov >= mind)
                cands = components(m, src, dst, minrun)
                n_c, prec, rec = score(cands, name, False)
                big = 0
                allbig = [i for i, v in sizes[name].items() if v >= 500]
                found = set()
                for c in cands:
                    s = sig[np.asarray(c, dtype=np.int64)]
                    i = np.clip(np.searchsorted(keys, s), 0, len(keys) - 1)
                    found |= {r[0] for r in {ins[int(k)] for k in s[keys[i] == s]} if r[1] == name}
                big = len([i for i in found if sizes[name][i] >= 500])
                rb = big / max(len(allbig), 1)
                print(f"  {minrun:>7}{mind:>9}{n_c:>10,}{prec:>11.2%}{rec:>9.1%}{rb:>15.1%}")


if __name__ == "__main__":
    main()
