#!/usr/bin/env python3
"""What is each DEL candidate actually made of, stratified by its size?

The genome-wide scan says the dosage detector is 5% precise. That number is useless on its own: a
candidate's size decides what it can possibly be. A SNV kills 2-3 reference cores, so it can only ever
produce a SHORT run; a 500 bp deletion removes ~6 cores and a 10 kb one ~130. So the question is not
"how precise is the detector" but "at what run length does it stop finding point mutations".

Every candidate is attributed to the truth record it shares the most cores with, and candidates are
binned by run length (in nodes) and by the implied bp. No truth is used to build the candidates.

usage: scan_strat.py <graph_v3.bin> <gt.jsonl> [LO=0.1] [MINC=3] [REP=1.75]
"""
import collections
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
from scan_dosage import components, load_graph

SMALL = 50          # a truth variant below this is a SMV (SNV / small indel) by the new vocabulary
BIG = 500           # the size the user cares about


def truth_index(gt):
    """sig -> list of (record id, svtype, svlen). A core can be destroyed by more than one variant."""
    idx = collections.defaultdict(list)
    recs = {}
    for line in open(gt):
        r = json.loads(line)
        rid, t, sl = r["id"], r["svtype"], abs(r["svlen"])
        recs[rid] = (t, sl)
        for k in ("over", "over_rc", "inside", "inside_rc"):
            for c in r.get(k) or []:
                idx[c[2]].append(rid)
    return idx, recs


def klass(t, sl):
    if t in ("SNV", "INDEL") or sl < SMALL:
        return "SMV(<50bp)"
    if t == "DEL":
        return f"DEL>={BIG}" if sl >= BIG else f"DEL 50-{BIG}"
    if t in ("INS", "DUP", "INV"):
        return t
    return "CPX"


def main():
    gpath, gt = sys.argv[1], sys.argv[2]
    LO = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
    MINC = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    REP = float(sys.argv[5]) if len(sys.argv) > 5 else 1.75

    sig, ncov, tcov, src, dst = load_graph(gpath)
    dose = np.where(ncov > 0, tcov.astype(np.float64) / np.maximum(ncov, 1), np.inf)
    depleted = (ncov > 0) & (dose < LO)
    print(f"graph {gpath}: K={len(sig):,}  LO={LO} MINC={MINC} REP={REP}")

    # Coverage-derived repeat mask -- the reference-free replacement for the blacklist. A core that
    # occurs at C places in the genome collects ~C x the depth, in BOTH samples: that is not a tumour
    # excess, it is copy number, and it is exactly what breaks a deletion's run. Such a core is made
    # TRANSPARENT: it neither starts a call nor ends one, the run passes straight through it.
    backbone = float(np.median(ncov[(ncov > 0) & (tcov > 0)]))
    copies = (ncov.astype(np.float64) + tcov) / (2.0 * backbone)
    repeat = copies > REP
    print(f"backbone depth (median ncov of shared nodes): {backbone:.0f}")
    print(f"repeat nodes (>{REP}x backbone in n+t): {int(repeat.sum()):,} ({repeat.mean():.2%})")
    print(f"tumour-depleted nodes: {int(depleted.sum()):,} ({depleted.mean():.2%})")

    # a run may traverse repeats, but must be carried by at least MINC genuinely depleted nodes
    raw = components(depleted | repeat, src, dst, MINC)
    cands = []
    for c in raw:
        d = np.asarray(c)[depleted[np.asarray(c, dtype=np.int64)]]
        if len(d) >= MINC:
            cands.append(d)
    print(f"candidates: {len(cands):,}\n")

    idx, recs = truth_index(gt)
    keys = np.array(sorted(idx.keys()), dtype=np.uint64)

    # attribute every candidate to the truth record it shares the most cores with
    rows = []
    for c in cands:
        s = sig[np.asarray(c, dtype=np.int64)]
        i = np.clip(np.searchsorted(keys, s), 0, len(keys) - 1)
        hits = collections.Counter()
        for k in s[keys[i] == s]:
            hits.update(idx[int(k)])
        if hits:
            rid, n = hits.most_common(1)[0]
            t, sl = recs[rid]
            rows.append((len(c), klass(t, sl), rid, sl))
        else:
            rows.append((len(c), "NOTHING", None, 0))

    bins = [(3, 4), (5, 6), (7, 9), (10, 14), (15, 29), (30, 59), (60, 10**9)]
    cols = [f"DEL>={BIG}", f"DEL 50-{BIG}", "SMV(<50bp)", "INS", "DUP", "INV", "CPX", "NOTHING"]
    print(f"=== what a candidate of each run length is made of ===")
    print(f"{'run (nodes)':<13}{'cands':>9}" + "".join(f"{c:>13}" for c in cols))
    for lo, hi in bins:
        sub = [r for r in rows if lo <= r[0] <= hi]
        if not sub:
            continue
        cnt = collections.Counter(r[1] for r in sub)
        n = len(sub)
        print(f"{f'{lo}-{hi if hi<10**9 else ""}':<13}{n:>9}" + "".join(f"{cnt[c]/n:>13.1%}" for c in cols))

    # the call rule the user asked about: call DEL when the run is at least X nodes
    allbig = {rid for rid, (t, sl) in recs.items() if t == "DEL" and sl >= BIG}
    print(f"\n=== call DEL when run >= X nodes; scored against the {len(allbig):,} truth DELs >= {BIG} bp ===")
    # Two precisions, because under the DEL/INS/DUP/INV/CPX/SMV vocabulary they answer different
    # questions. "as DEL" = the call is a deletion. "as any SV" also credits calls landing on a complex
    # event's deleted block: the somatic loss is real and correctly located, only the event-level label
    # is CPX rather than DEL. A call on an SMV or on nothing is wrong either way.
    print(f"{'X':>4}{'calls':>9}{'prec as DEL':>13}{'prec as any SV':>16}{'recall':>9}")
    for X in (3, 4, 5, 6, 8, 10, 12, 15, 20, 30):
        sub = [r for r in rows if r[0] >= X]
        if not sub:
            continue
        good = [r for r in sub if r[1] == f"DEL>={BIG}"]
        anysv = [r for r in sub if r[1] not in ("SMV(<50bp)", "NOTHING")]
        found = {r[2] for r in good}
        print(f"{X:>4}{len(sub):>9}{len(good)/len(sub):>13.1%}{len(anysv)/len(sub):>16.1%}"
              f"{len(found)/max(len(allbig),1):>9.1%}")

    # the mirror-twin tax: every locus enters the graph twice, so a real event yields two candidates
    dbl = collections.Counter(r[2] for r in rows if r[2])
    if dbl:
        per = np.array(list(dbl.values()))
        print(f"\nmirror check: a matched truth record is hit by {np.median(per):.1f} candidates (median) --"
              f" both strand families produce their own run, so half of every count above is the twin.")


if __name__ == "__main__":
    main()
