#!/usr/bin/env python3
"""Two things the stratified table could not tell apart, measured separately.

(1) STRAND. Every read enters the graph both ways, so a locus has a forward-family core and a mirror
    core with a different signature, and a deletion should therefore produce TWO runs. scan_strat's
    mirror check went from 2 candidates per truth record to 1.0 after homopolymer normalisation. This
    counts, per truth deletion, how many candidates touch its `inside` cores and how many touch its
    `inside_rc` cores -- if one family stopped producing runs, it shows up here as a zero.

(2) ATTRIBUTION. scan_strat gives a candidate to whichever truth record it shares the most cores with.
    With 118k SNVs against 3.5k deletions, a deletion's run can be handed to an SNV that happens to sit
    in it. This asks the question directly instead: of the candidates that winner-take-all calls SMV,
    how many actually contain a core from inside a deletion >= 500 bp?

usage: check_runs.py <graph_v3.bin> <gt.jsonl>
"""
import collections
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
from scan_dosage import components, load_graph

BIG = 500


def main():
    gpath, gt = sys.argv[1], sys.argv[2]
    sig, ncov, tcov, src, dst = load_graph(gpath)
    dose = np.where(ncov > 0, tcov.astype(np.float64) / np.maximum(ncov, 1), np.inf)
    depleted = (ncov > 0) & (dose < 0.1)
    print(f"graph {gpath}: K={len(sig):,} | depleted {int(depleted.sum()):,}")
    cands = components(depleted, src, dst, 3)
    print(f"candidates (>=3 nodes): {len(cands):,}\n")

    # --- truth, keeping the two strand families apart ---
    fw, rc, smv = {}, {}, set()          # sig -> record id, per family
    size = {}
    for line in open(gt):
        r = json.loads(line)
        rid, t, sl = r["id"], r["svtype"], abs(r["svlen"])
        if t == "DEL" and sl >= BIG:
            size[rid] = sl
            for c in (r.get("inside") or []):
                fw[c[2]] = rid
            for c in (r.get("inside_rc") or []):
                rc[c[2]] = rid
        elif t in ("SNV", "INDEL") or sl < 50:
            for k in ("over", "over_rc", "inside", "inside_rc"):
                for c in (r.get(k) or []):
                    smv.add(c[2])
    print(f"truth DEL>={BIG}: {len(size):,} records | fw cores {len(fw):,} | rc cores {len(rc):,}")

    kf = np.array(sorted(fw), dtype=np.uint64)
    kr = np.array(sorted(rc), dtype=np.uint64)
    ks = np.array(sorted(smv), dtype=np.uint64)

    def look(s, keys):
        if len(keys) == 0:
            return np.zeros(len(s), bool)
        i = np.clip(np.searchsorted(keys, s), 0, len(keys) - 1)
        return keys[i] == s

    # --- (1) per truth deletion: how many candidates touch each family ---
    hit_fw = collections.Counter()
    hit_rc = collections.Counter()
    # --- (2) candidates that contain a deletion core, by run length ---
    rows = []
    for c in cands:
        s = sig[np.asarray(c, dtype=np.int64)]
        mf, mr, ms = look(s, kf), look(s, kr), look(s, ks)
        ids_f = {fw[int(x)] for x in s[mf]}
        ids_r = {rc[int(x)] for x in s[mr]}
        for rid in ids_f:
            hit_fw[rid] += 1
        for rid in ids_r:
            hit_rc[rid] += 1
        rows.append((len(c), int(mf.sum() + mr.sum()), int(ms.sum()), bool(ids_f or ids_r)))

    both = sum(1 for rid in size if hit_fw[rid] and hit_rc[rid])
    only_f = sum(1 for rid in size if hit_fw[rid] and not hit_rc[rid])
    only_r = sum(1 for rid in size if hit_rc[rid] and not hit_fw[rid])
    none = sum(1 for rid in size if not hit_fw[rid] and not hit_rc[rid])
    n = max(len(size), 1)
    print("\n=== (1) which strand family produced a run, per truth DEL >= 500bp ===")
    print(f"  both families : {both:>5} ({both/n:6.1%})")
    print(f"  forward only  : {only_f:>5} ({only_f/n:6.1%})")
    print(f"  mirror only   : {only_r:>5} ({only_r/n:6.1%})")
    print(f"  neither       : {none:>5} ({none/n:6.1%})")
    tot = [hit_fw[rid] + hit_rc[rid] for rid in size if hit_fw[rid] or hit_rc[rid]]
    if tot:
        print(f"  candidates per detected deletion: median {np.median(tot):.1f}, mean {np.mean(tot):.2f}")

    print("\n=== (2) what candidates of each run length actually CONTAIN ===")
    print(f"{'run':>9}{'cands':>9}{'has DEL core':>14}{'has SMV core':>14}{'DEL-only':>10}{'SMV-only':>10}{'both':>8}")
    for lo, hi in ((3, 4), (5, 6), (7, 9), (10, 14), (15, 29), (30, 59), (60, 10**9)):
        sub = [r for r in rows if lo <= r[0] <= hi]
        if not sub:
            continue
        nd = sum(1 for r in sub if r[1])
        ns = sum(1 for r in sub if r[2])
        do = sum(1 for r in sub if r[1] and not r[2])
        so = sum(1 for r in sub if r[2] and not r[1])
        bo = sum(1 for r in sub if r[1] and r[2])
        m = len(sub)
        lab = f"{lo}-{hi if hi < 10**9 else ''}"
        print(f"{lab:>9}{m:>9}{nd/m:>14.1%}{ns/m:>14.1%}{do/m:>10.1%}{so/m:>10.1%}{bo/m:>8.1%}")


    # --- threshold curve: ONE definition throughout ---------------------------------------------
    # call = a run of >= X nodes. precision = the call contains a core from inside a DEL >= 500bp.
    # recall = the fraction of the truth DEL >= 500bp that at least one call contains a core of.
    # Both come from the same pass over the same candidates, so the columns are comparable.
    print("\n=== threshold curve (call = run >= X nodes) ===")
    print(f"{'X':>4}{'calls':>10}{'precision':>12}{'recall':>10}{'deletions found':>18}")
    for X in (3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30):
        hf, hr = collections.Counter(), collections.Counter()
        ncall = ngood = 0
        for c in cands:
            if len(c) < X:
                continue
            ncall += 1
            s_ = sig[np.asarray(c, dtype=np.int64)]
            mf, mr = look(s_, kf), look(s_, kr)
            ids = {fw[int(x)] for x in s_[mf]} | {rc[int(x)] for x in s_[mr]}
            if ids:
                ngood += 1
                for rid in ids:
                    hf[rid] += 1
        found = len(hf)
        print(f"{X:>4}{ncall:>10,}{(ngood/ncall if ncall else 0):>12.1%}"
              f"{found/n:>10.1%}{f'{found}/{len(size)}':>18}")


if __name__ == "__main__":
    main()
