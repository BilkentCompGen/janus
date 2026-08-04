#!/usr/bin/env python3
"""How much of the clean junction signal does the graph actually carry?

The clean signal is built by tools/junctions off the two FASTAs: for every novel adjacency a truth
variant creates, the +/-500 bp of TUMOUR sequence around the join, cut from out.fa and parsed with the
graph's own extract(). This script asks, per variant type, how much of that survives into graph.bin.

Four questions, in order, because each one can only fail after the previous one passed:

  SIGNAL    does the junction produce any core the NORMAL GENOME cannot?  A core is novel when its
            signature is absent from the whole-reference signature map. Nothing below this line can be
            recovered by any caller -- a junction with no novel core is invisible to a core graph, and
            that is a property of the sequence, not of the pipeline.
  IN GRAPH  is at least one of those novel cores a node at all? Failing here is read sampling plus the
            min_mult floor.
  SOMATIC   is a novel core a node with ncov == 0, i.e. derived colour 2? A novel core that picked up
            normal reads is not somatic to the graph, whatever the truth says.
  ANCHORED  is the novel run joined by a graph EDGE to a real (normal-carrying) core on that side? This
            is what turns "something is here" into "it starts HERE and ends THERE", and it is where
            start and end become separable for every type -- including a deletion, whose two
            breakpoints are one join and therefore one window.

Both strand families are built (the graph inserts every read forward and reverse-complemented, so a
locus has a mirror twin with different signatures and no core in common). A junction counts as covered
if EITHER family carries it: one locus is one event, not two.

Variants strictly inside somebody else's deletion are not in the tumour genome at all. They are
reported apart and kept out of every denominator -- the difference between impossible and missed is
the measurement.

False junctions: every tumour-only node in the graph is looked up in the TUMOUR genome's own signature
map, so "read error" and "real somatic sequence the truth list does not put at a junction" stop being
the same bucket.

usage: junction_cov.py <label> <graph.bin> <junctions.jsonl> <ref_copynum.bin> <tumour_copynum.bin>
"""
import collections
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
from scan_dosage import components, load_graph

TYPES_FIRST = ("SNV", "DEL", "INS", "DUP", "INV")   # then everything complex, alphabetically


def load_sigs(path):
    """The signature column of a copynum map: { u64 N; {u64 sig, u32 count}[N] }, sorted by sig."""
    n = int(np.fromfile(path, dtype=np.uint64, count=1)[0])
    raw = np.fromfile(path, dtype=np.uint8, offset=8, count=n * 12).reshape(n, 12)
    return raw[:, :8].copy().view(np.uint64).ravel()


def member(query, sorted_pool):
    """Which of `query` are in `sorted_pool`, and where (-1 when absent)."""
    if not len(query):
        return np.zeros(0, bool), np.zeros(0, np.int64)
    i = np.searchsorted(sorted_pool, query)
    ok = (i < len(sorted_pool)) & (sorted_pool[np.minimum(i, len(sorted_pool) - 1)] == query)
    return ok, np.where(ok, i, -1)


def flatten(path):
    """One row per core of every located junction, plus the per-junction metadata."""
    meta, sig, jid, fam, pos = [], [], [], [], []
    for line in open(path):
        d = json.loads(line)
        j = len(meta)
        meta.append((d["svtype"], d["svid"], d["side"], d["svlen"], d["located"], d["absent"],
                     min(d["pos1"], d["pos2"]), max(d["pos1"], d["pos2"])))
        if not d["located"]:
            continue
        for f, key in ((0, "cores"), (1, "cores_rc")):
            for k, (st, en, s) in enumerate(d[key]):
                sig.append(s); jid.append(j); fam.append(f); pos.append(k)
    return (meta, np.array(sig, dtype=np.uint64), np.array(jid, dtype=np.int64),
            np.array(fam, dtype=np.int8), np.array(pos, dtype=np.int32))


def edge_codes(src, dst):
    """Both directions of every edge as one sorted u64 key, so membership is a binary search."""
    a = (src.astype(np.uint64) << np.uint64(32)) | dst.astype(np.uint64)
    b = (dst.astype(np.uint64) << np.uint64(32)) | src.astype(np.uint64)
    return np.sort(np.concatenate([a, b]))


def main():
    label, gpath, jpath, refpath, tumpath = sys.argv[1:6]
    gsig, ncov, tcov, esrc, edst = load_graph(gpath)
    refsig, tumsig = load_sigs(refpath), load_sigs(tumpath)
    ecode = edge_codes(esrc, edst)
    meta, sig, jid, fam, pos = flatten(jpath)

    # ---- per-core facts ----
    in_ref, _ = member(sig, refsig)
    in_g, gi = member(sig, gsig)
    node = np.where(in_g, gi, 0)
    novel = ~in_ref                                  # the normal genome cannot produce this core
    somatic = in_g & (ncov[node] == 0)               # derived colour 2
    real = in_g & (ncov[node] > 0)                   # carries normal reads: usable as an anchor

    # ---- per (junction, family) run of novel cores ----
    key = jid * 2 + fam
    order = np.lexsort((pos, key))
    key = key[order]
    nv, sm, rl, ig, nd = novel[order], somatic[order], real[order], in_g[order], node[order]
    bounds = np.flatnonzero(np.r_[True, key[1:] != key[:-1], True])

    nj = len(meta)
    sig_ok = np.zeros(nj, bool); graph_ok = np.zeros(nj, bool); som_ok = np.zeros(nj, bool)
    left_ok = np.zeros(nj, bool); right_ok = np.zeros(nj, bool)
    qa, qb, qtag = [], [], []                        # boundary pairs to test for an edge
    for b in range(len(bounds) - 1):
        lo, hi = bounds[b], bounds[b + 1]
        j = int(key[lo] // 2)
        idx = np.flatnonzero(nv[lo:hi])
        if not len(idx):
            continue
        sig_ok[j] = True
        graph_ok[j] |= bool(ig[lo:hi][idx].any())
        som_ok[j] |= bool(sm[lo:hi][idx].any())
        i0, i1 = int(idx[0]), int(idx[-1])
        if i0 > 0 and rl[lo + i0 - 1] and sm[lo + i0]:
            qa.append(nd[lo + i0 - 1]); qb.append(nd[lo + i0]); qtag.append((j, 0))
        if i1 + 1 < hi - lo and rl[lo + i1 + 1] and sm[lo + i1]:
            qa.append(nd[lo + i1]); qb.append(nd[lo + i1 + 1]); qtag.append((j, 1))

    if qa:
        q = (np.array(qa, np.uint64) << np.uint64(32)) | np.array(qb, np.uint64)
        i = np.searchsorted(ecode, q)
        hit = (i < len(ecode)) & (ecode[np.minimum(i, len(ecode) - 1)] == q)
        for (j, side), h in zip(qtag, hit.tolist()):
            if h:
                (left_ok if side == 0 else right_ok)[j] = True

    report(label, meta, sig_ok, graph_ok, som_ok, left_ok, right_ok)
    false_junctions(label, gsig, ncov, tcov, refsig, tumsig, sig[novel], esrc, edst)


def row(name, n, *cols):
    cells = "".join(f"{c:>8d} {100 * c / n if n else 0:>5.1f}%" for c in cols)
    return f"  {name:<26s} {n:>8d}{cells}"


def order_types(seen):
    rest = sorted(t for t in seen if t not in TYPES_FIRST)
    return [t for t in TYPES_FIRST if t in seen] + rest


def report(label, meta, sig_ok, graph_ok, som_ok, left_ok, right_ok):
    by = collections.defaultdict(list)
    for j, m in enumerate(meta):
        by[m[0]].append(j)

    print(f"\n=== {label} :: JUNCTION LEVEL (one row per novel adjacency) ===")
    print(f"  {'type':<26s} {'located':>8s} {'SIGNAL':>14s} {'IN GRAPH':>14s} {'SOMATIC':>14s}"
          f" {'ANCH BOTH':>14s} {'ANCH LEFT':>14s} {'ANCH RIGHT':>14s} {'ANCH NONE':>14s}")
    for t in order_types(by):
        js = [j for j in by[t] if meta[j][4]]
        n = len(js)
        if not n:
            continue
        a = np.array(js)
        both = int((left_ok[a] & right_ok[a]).sum())
        lo = int((left_ok[a] & ~right_ok[a]).sum())
        ro = int((~left_ok[a] & right_ok[a]).sum())
        none = int((~left_ok[a] & ~right_ok[a]).sum())
        print(row(t, n, int(sig_ok[a].sum()), int(graph_ok[a].sum()), int(som_ok[a].sum()), both, lo, ro, none))

    absent = collections.Counter(meta[j][0] for j in range(len(meta)) if meta[j][5])
    lost = collections.Counter(meta[j][0] for j in range(len(meta)) if not meta[j][4] and not meta[j][5])
    print(f"  not in the tumour genome at all (inside another deletion): {dict(absent)}")
    print(f"  could not be located in out.fa: {dict(lost)}")

    # ---- variant level: does the START junction and the END junction each survive? ----
    # The simulator labels an adjacency by the symbol it touches, so a deletion or a tandem
    # duplication is `both` (its two breakpoints ARE one join) and an insertion or inversion gets a
    # `start` and an `end`. A multi-block complex event can have neither symbol on any of its
    # adjacencies -- there the leftmost and rightmost reference anchors define the two ends instead.
    var = collections.defaultdict(list)
    for j, m in enumerate(meta):
        if m[4]:
            var[m[1]].append(j)

    print(f"\n=== {label} :: VARIANT LEVEL (does each end of the variant have a usable junction?) ===")
    print("  a junction counts when its novel cores are in the graph as tumour-only nodes (SOMATIC)")
    print(f"  {'type':<26s} {'variants':>8s} {'jn/var':>7s} {'BOTH':>14s} {'START ONLY':>14s}"
          f" {'END ONLY':>14s} {'NEITHER':>14s} {'ALL JN':>14s}")
    tab = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0, 0])
    for svid, js in var.items():
        t = meta[js[0]][0]
        st = [j for j in js if meta[j][2] in ("start", "both")]
        en = [j for j in js if meta[j][2] in ("end", "both")]
        if not st:
            lo = min(meta[j][6] for j in js); st = [j for j in js if meta[j][6] == lo]
        if not en:
            hi = max(meta[j][7] for j in js); en = [j for j in js if meta[j][7] == hi]
        a, b = any(som_ok[j] for j in st), any(som_ok[j] for j in en)
        c = tab[t]
        c[0] += 1
        c[6] += len(js)
        c[1 if (a and b) else 2 if a else 3 if b else 4] += 1
        c[5] += int(all(som_ok[j] for j in js))
    for t in order_types(tab):
        n, bo, so, eo, no, alljn, tot = tab[t]
        print(f"  {t:<26s} {n:>8d} {tot / n:>7.1f}"
              + "".join(f"{v:>8d} {100 * v / n:>5.1f}%" for v in (bo, so, eo, no, alljn)))


def false_junctions(label, gsig, ncov, tcov, refsig, tumsig, junction_novel, esrc, edst):
    """Every tumour-only node, split by what -- if anything -- in the two genomes can explain it.

    A call is a connected RUN of such nodes, so each class is also counted in components: one false
    node is noise, a hundred joined into a run is a false junction call waiting to happen."""
    m = (ncov == 0) & (tcov > 0)
    s = gsig[m]
    n = len(s)
    in_tum, _ = member(s, tumsig)
    in_ref, _ = member(s, refsig)
    in_jn, _ = member(s, np.unique(junction_novel))
    ids = np.flatnonzero(m)
    print(f"\n=== {label} :: TUMOUR-ONLY NODES (ncov == 0, tcov > 0) ===")
    print(f"  {'total':<46s} {n:>10d}")
    print(f"  {'':<46s} {'nodes':>10s} {'share':>8s} {'runs':>9s} {'max run':>8s} {'tcov med':>9s}")
    for name, mask in (
        ("real: at a truth junction (+/-500 bp of a join)", in_jn),
        ("real: tumour sequence away from any junction", in_tum & ~in_ref & ~in_jn),
        ("both genomes have it -- normal-read dropout", in_ref & ~in_jn),
        ("NEITHER genome has it -- false", ~in_tum & ~in_ref & ~in_jn),
    ):
        c = int(mask.sum())
        runs = components(np.isin(np.arange(len(gsig)), ids[mask]), esrc, edst, 1) if c else []
        big = max((len(r) for r in runs), default=0)
        print(f"  {name:<46s} {c:>10d} {100 * c / n if n else 0:>7.2f}% {len(runs):>9d} {big:>8d} "
              f"{int(np.median(tcov[ids[mask]])) if c else 0:>9d}")


if __name__ == "__main__":
    main()
