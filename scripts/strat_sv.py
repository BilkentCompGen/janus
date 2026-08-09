#!/usr/bin/env python3
"""DEL and DUP, unique versus repetitive reference, two builds, one rule.

Every parameter is fixed before any score is computed, and each is fixed by a measurement that does
not look at the final precision/recall -- otherwise the comparison is just tuning.

  node test      DEL: tcov/ncov < 0.1     DUP: tcov/ncov > 1.5, both counts > 0
                 LO is not tuned -- a deleted core's dose is bimodal, so anything from 0.05 to 0.30
                 gives the same TPR (0.967 / 0.974) and the same FPR to three digits. HI is a real
                 trade: 1.3 doubles the FPR, 1.75 drops node TPR to 0.69, and the run rule needs
                 CONSECUTIVE passing nodes so a gap costs multiplicatively. 1.5 is where they meet.
  depth floor    ncov >= 0.5 x backbone, backbone = median ncov over nodes seen in both colours.
                 The largest floor still costing under 1% of node-level TPR in all four
                 (build x type) cells; at 0.75x TPR collapses. A gate on the NORMAL sample alone.
  run rule       a call is a connected component of passing nodes, and its threshold is stated as
                 the smallest event it is allowed to claim, in bp. X consecutive cores span
                 (X-1)*stride + corelen, and the builds share neither constant -- measured on the
                 truth interiors, 67.3 / 151 bp uncompressed against 97.8 / 217 bp compressed, so
                 the same X is 1.45x more genome under HPC and the shortest call it can physically
                 make is 413 bp against 286 bp. X is the CEILING of the inverse, so a row never
                 admits a call whose own footprint is under the row's bp -- rounding to nearest
                 silently judged the two arms at footprints up to 51 bp apart.
  scoring        containment, and class- and size-consistent. A call of class C at row B counts as
                 correct only if its node set meets an interior core of a truth event that is of the
                 target type, of class C, and at least B bp. Recall at the same cell is over the
                 truth events of that type and class with svlen >= B, found by the calls of that
                 class. Precision and recall therefore range over the same population; mixing them
                 is what made an earlier version report 100% precision and 21.7% recall for
                 repetitive sequence when every one of those calls was sitting on unique sequence.
  region class   unique = the signature occurs once in the reference, problem = twice or more, taken
                 as the median over a locus's signatures. A truth event's class is assigned ONCE,
                 from the uncompressed arm's map, and both columns use that same partition -- the
                 class is a property of the genome, and letting each arm draw its own line puts 3 of
                 371 events in different tables and makes the columns score different populations.
                 A call's class is per-arm, because a call is made of that arm's own nodes.

Note both builds enter every read twice, as sequenced and reverse-complemented, so a locus yields
two mirror candidates: the call counts below are fragments, not loci, and the median event is
covered by exactly 2 of them in both arms.

usage: strat_sv.py <labelA> <graphA> <gtA> <cnA> <labelB> <graphB> <gtB> <cnB>
       (A is the arm whose reference map defines the truth partition)
"""
import collections
import json
import math
import sys

import numpy as np

sys.path.insert(0, "scripts")
from scan_dosage import components, load_graph

BIG = 500                    # an event smaller than this is not what the dosage detector is for
LO, HI = 0.1, 1.5            # node tests
FLOOR = 0.5                  # x backbone, on ncov
BPS = (500, 750, 1000, 1500, 2000, 3000, 5000)   # smallest event a call is allowed to claim
ALT = 0.25                   # sensitivity rule: problem if this share of cores are multi-copy


def load_cn(path):
    n = int(np.fromfile(path, dtype=np.uint64, count=1)[0])
    raw = np.fromfile(path, dtype=np.uint8, offset=8, count=n * 12).reshape(n, 12)
    return raw[:, :8].copy().view(np.uint64).ravel(), raw[:, 8:].copy().view(np.uint32).ravel()


def cn_of(sigs, cs, cc):
    q = np.asarray(sigs, dtype=np.uint64)
    if not len(q):
        return np.zeros(0, np.uint32)
    i = np.clip(np.searchsorted(cs, q), 0, len(cs) - 1)
    return np.where(cs[i] == q, cc[i], 0)


def shape(cn):
    """(median copy number, share of copies that are multi-copy) over the signatures the reference
    actually knows. A locus none of whose signatures are in the map is unique by default: a
    signature the reference never produces is novel or error sequence, not a collapsed repeat."""
    cn = cn[cn > 0]
    if not len(cn):
        return 1.0, 0.0
    return float(np.median(cn)), float((cn >= 2).mean())


def klass(sh, alt=False):
    return "problem" if (sh[1] > ALT if alt else sh[0] >= 2) else "unique"


def build(label, gpath, gt, cnpath, shared=None):
    sig, ncov, tcov, src, dst = load_graph(gpath)
    cs, cc = load_cn(cnpath)
    backbone = float(np.median(ncov[(ncov > 0) & (tcov > 0)]))
    floor = max(int(math.ceil(FLOOR * backbone)), 1)
    dose = np.where(ncov > 0, tcov.astype(np.float64) / np.maximum(ncov, 1), np.inf)
    deep = ncov >= floor
    masks = {"DEL": deep & (dose < LO),
             "DUP": deep & (tcov > 0) & (dose > HI)}

    # ---- truth: interior cores per event, its shape on THIS map, and this build's core geometry --
    ev, cores, shp, stride_num, stride_den, clen = {}, {}, {}, 0, 0, []
    for line in open(gt):
        r = json.loads(line)
        t, sl = r["svtype"], abs(r["svlen"])
        if t not in ("DEL", "DUP") or sl < BIG:
            continue
        fw = r.get("inside") or []
        c = {x[2] for x in fw} | {x[2] for x in (r.get("inside_rc") or [])}
        if not c:
            continue
        ev[r["id"]] = (t, sl)
        cores[r["id"]] = c
        shp[r["id"]] = shape(cn_of(sorted(c), cs, cc))
        clen.append(np.array([x[1] - x[0] for x in fw], np.int64))
        if sl >= 5000:                      # stride from the big events, where end effects vanish
            stride_num += sl
            stride_den += len(fw)
    stride = stride_num / stride_den
    corelen = float(np.median(np.concatenate(clen)))
    shp = shared if shared is not None else shp        # one partition, shared by both arms

    # every core of every truth record, so a wrong call can be told apart from a blind one
    owner, target, small, midsv = collections.defaultdict(set), collections.defaultdict(set), set(), set()
    for line in open(gt):
        r = json.loads(line)
        pt = r["svtype"] in ("SNV", "INDEL") or abs(r["svlen"]) < 50
        for k in ("over", "over_rc", "inside", "inside_rc"):
            for c in r.get(k) or []:
                owner[c[2]].add(r["id"])
                (small if pt else midsv).add(c[2])
    for rid, c in cores.items():
        for s in c:
            target[s].add(rid)
    okeys = np.array(sorted(owner), dtype=np.uint64)
    tkeys = np.array(sorted(target), dtype=np.uint64)
    skeys = np.array(sorted(small), dtype=np.uint64)
    mkeys = np.array(sorted(midsv - small), dtype=np.uint64)

    print(f"# {label}: K={len(sig):,} backbone={backbone:.0f} floor=ncov>={floor} "
          f"stride={stride:.1f} corelen={corelen:.0f} -> shortest claimable event "
          f"{2*stride+corelen:.0f}bp", flush=True)

    def has(s, keys):
        if not len(keys):
            return False
        j = np.clip(np.searchsorted(keys, s), 0, len(keys) - 1)
        return bool((keys[j] == s).any())

    out, passfrac, contig = {}, {}, {}
    for t, mask in masks.items():
        # how spatially sticky is this node test? P(a neighbour of a passing node also passes).
        # This, not the depth floor, is what decides how long a false run gets.
        e_any = mask[src] | mask[dst]
        contig[t] = float((mask[src] & mask[dst]).sum() / max(e_any.sum(), 1))
        cands = components(mask, src, dst, 3)
        rows = []
        for c in cands:
            s = sig[np.asarray(c, dtype=np.int64)]
            i = np.clip(np.searchsorted(tkeys, s), 0, len(tkeys) - 1)
            hit = set()
            for k in s[tkeys[i] == s]:
                hit |= {r for r in target[int(k)] if ev[r][0] == t}
            rows.append(dict(n=len(c), sh=shape(cn_of(s, cs, cc)), hit=hit,
                             any=has(s, okeys), small=has(s, skeys), mid=has(s, mkeys)))
        out[t] = rows
        for rid, (tt, sl) in ev.items():                # can the node test see it at all?
            if tt != t:
                continue
            q = np.array(sorted(cores[rid]), dtype=np.uint64)
            i = np.clip(np.searchsorted(sig, q), 0, len(sig) - 1)
            j = i[sig[i] == q]
            passfrac[rid] = float(mask[j].mean()) if len(j) else 0.0
        print(f"#   {t}: {int(mask.sum()):,} nodes pass ({mask.mean():.2%}), "
              f"P(neighbour also passes)={contig[t]:.3f} -> {len(cands):,} candidates", flush=True)

    return dict(rows=out, ev=ev, shp=shp, stride=stride, corelen=corelen, label=label,
                passfrac=passfrac, contig=contig, K=len(sig), floor=floor,
                below_floor=int((ncov < floor).sum()))


def nodes_for(b, B):
    """The smallest node count whose footprint (X-1)*stride + corelen reaches B bp. Ceiling, not
    rounding: a row must never admit a call that cannot span the event it claims."""
    return max(int(math.ceil(1 + (B - b["corelen"]) / b["stride"])), 3)


def cell(b, t, strat, B, alt=False):
    """Everything about one cell, all of it over the same population: calls of class `strat`,
    truth events of type `t`, class `strat`, at least B bp."""
    X = nodes_for(b, B)
    den = {r for r, (tt, sl) in b["ev"].items()
           if tt == t and sl >= B and klass(b["shp"][r], alt) == strat}
    mine = [r for r in b["rows"][t] if r["n"] >= X and klass(r["sh"], alt) == strat]
    good = [r for r in mine if r["hit"] & den]
    found = set().union(*[r["hit"] for r in good]) & den if good else set()
    return dict(X=X, calls=len(mine), prec=(len(good) / len(mine) if mine else None),
                found=len(found), elig=len(den),
                # what the wrong ones are, as an exclusive partition
                loose=(sum(1 for r in mine if r["hit"]) / len(mine) if mine else None),
                pt=(sum(1 for r in mine if not r["hit"] and r["small"]) / len(mine) if mine else None),
                other=(sum(1 for r in mine if not r["hit"] and not r["small"] and r["any"]) / len(mine) if mine else None),
                blank=(sum(1 for r in mine if not r["any"]) / len(mine) if mine else None),
                perloc=(len(good) / len(found) if found else None),
                permed=(float(np.median(list(collections.Counter(
                    e for r in good for e in (r["hit"] & den)).values()))) if found else None))


def pc(v):
    return "      -" if v is None else f"{v:>7.1%}"


def table(t, strat, A, B_, alt=False):
    print(f"\n### {t} -- {strat} region"
          f"{'   [ALT class rule: >' + str(int(ALT*100)) + '% of cores multi-copy]' if alt else ''}")
    print(f"{'min event':>10} | {A['label']:>34} | {B_['label']:>34}")
    print(f"{'bp':>10} | {'X':>4}{'calls':>9}{'prec':>8}{'recall':>8}{'found/elig':>13} | "
          f"{'X':>4}{'calls':>9}{'prec':>8}{'recall':>8}{'found/elig':>13}")
    for B in BPS:
        line = f"{B:>10} |"
        for b in (A, B_):
            c = cell(b, t, strat, B, alt)
            rec = "      -" if not c["elig"] else f"{c['found']/c['elig']:>7.1%}"
            frac = "{}/{}".format(c["found"], c["elig"])
            line += f" {c['X']:>4}{c['calls']:>9,}{pc(c['prec'])}{rec}{frac:>13} |"
        print(line)


def main():
    a = sys.argv[1:]
    A = build(*a[0:4])
    B_ = build(*a[4:8], shared=A["shp"])
    A["shp"] = A["shp"]

    n_pro = {t: sum(1 for r, (tt, sl) in A["ev"].items() if tt == t and klass(A["shp"][r]) == "problem")
             for t in ("DEL", "DUP")}
    n_uni = {t: sum(1 for r, (tt, sl) in A["ev"].items() if tt == t and klass(A["shp"][r]) == "unique")
             for t in ("DEL", "DUP")}
    print(f"\nthe bp axis: X is the smallest run that spans the row, so the two arms are held to")
    print(f"footprints that differ only by node granularity -- printed here so it is auditable")
    print(f"{'row bp':>8} | " + " | ".join(f"{b['label']:>22}" for b in (A, B_)))
    print(f"{'':>8} | " + " | ".join(f"{'X':>6}{'actual footprint':>16}" for _ in (A, B_)))
    for B in BPS:
        line = f"{B:>8} |"
        for b in (A, B_):
            X = nodes_for(b, B)
            fp = "{:.0f} bp".format((X - 1) * b["stride"] + b["corelen"])
            line += f" {X:>6}{fp:>16} |"
        print(line)

    print(f"\nshared truth partition (from {A['label']}): "
          f"DEL {n_uni['DEL']} unique / {n_pro['DEL']} problem, "
          f"DUP {n_uni['DUP']} unique / {n_pro['DUP']} problem")

    for t in ("DEL", "DUP"):
        for s in ("unique", "problem"):
            table(t, s, A, B_)

    print("\n\n" + "=" * 118)
    print("SUPPORTING MEASUREMENTS -- what the four tables are made of")
    print("=" * 118)

    print("\nS1. what the calls of a cell sit on -- an exclusive partition of the same call set")
    print(f"{'type':>5}{'region':>9}{'bp':>7} | " + " | ".join(f"{b['label']:>39}" for b in (A, B_)))
    print(f"{'':>21} | " + " | ".join(
        f"{'on target':>9}{'target,wrong size/class':>11}{'point var':>10}{'other SV':>9}{'nothing':>8}"
        for _ in (A, B_)))
    for t in ("DEL", "DUP"):
        for s in ("unique", "problem"):
            for B in (500, 1000, 2000):
                line = f"{t:>5}{s:>9}{B:>7} |"
                for b in (A, B_):
                    c = cell(b, t, s, B)
                    mis = None if c["prec"] is None else c["loose"] - c["prec"]
                    line += f" {pc(c['prec']):>9}{pc(mis):>11}{pc(c['pt']):>10}{pc(c['other']):>9}{pc(c['blank']):>8} |"
                print(line)

    print("\nS2. can the node test see the event at all -- share of an event's interior cores that")
    print("    pass, with NO connectivity required (per arm, on that arm's own cores)")
    print(f"{'type':>5}{'region':>9} | " + " | ".join(f"{b['label']:>44}" for b in (A, B_)))
    print(f"{'':>14} | " + " | ".join(
        f"{'events':>8}{'median svlen':>14}{'median pass':>13}{'>=50%':>8}" for _ in (A, B_)))
    for t in ("DEL", "DUP"):
        for s in ("unique", "problem"):
            ids = [r for r, (tt, sl) in A["ev"].items() if tt == t and klass(A["shp"][r]) == s]
            line = f"{t:>5}{s:>9} |"
            for b in (A, B_):
                mine = [i for i in ids if i in b["ev"]]
                sz = int(np.median([b["ev"][i][1] for i in mine]))
                v = np.array([b["passfrac"][i] for i in mine])
                line += f" {len(mine):>8}{sz:>14,}{np.median(v):>13.2f}{np.mean(v >= 0.5):>8.0%} |"
            print(line)

    print("\nS3. why a false run gets long: stickiness of the node test along the graph")
    print(f"{'':>7}{'DEL P(nbr passes)':>20}{'DUP P(nbr passes)':>20}{'nodes below floor':>20}")
    for b in (A, B_):
        bf = "{:.1%}".format(b["below_floor"] / b["K"])
        print(f"{b['label']:>7}{b['contig']['DEL']:>20.3f}{b['contig']['DUP']:>20.3f}{bf:>20}")

    print("\nS4. run length of the calls that meet no target event, in bp they claim")
    print(f"{'type':>5}{'region':>9} | " + " | ".join(f"{b['label']:>34}" for b in (A, B_)))
    print(f"{'':>14} | " + " | ".join(f"{'n':>9}{'p50':>8}{'p90':>8}{'max':>9}" for _ in (A, B_)))
    for t in ("DEL", "DUP"):
        for s in ("unique", "problem"):
            line = f"{t:>5}{s:>9} |"
            for b in (A, B_):
                L = np.array([(r["n"] - 1) * b["stride"] + b["corelen"]
                              for r in b["rows"][t] if klass(r["sh"]) == s and not r["hit"]])
                line += (f" {len(L):>9,}{np.percentile(L,50):>8.0f}{np.percentile(L,90):>8.0f}"
                         f"{L.max():>9.0f} |") if len(L) else f" {'0':>9}{'-':>8}{'-':>8}{'-':>9} |"
            print(line)

    print("\nS5. fragments per found event -- mirror twin = 2 of these, more is fragmentation.")
    print("    Distrust the MEAN: it is dominated by a couple of ~1 Mb deletions.")
    print(f"{'type':>5}{'region':>9}{'bp':>7} | " +
          " | ".join(f"{b['label']:>21}" for b in (A, B_)))
    print(f"{'':>21} | " + " | ".join(f"{'median':>10}{'mean':>11}" for _ in (A, B_)))
    for t in ("DEL", "DUP"):
        for s in ("unique", "problem"):
            for B in (1000, 2000):
                line = f"{t:>5}{s:>9}{B:>7} |"
                for b in (A, B_):
                    c = cell(b, t, s, B)
                    md = "-" if c["permed"] is None else f"{c['permed']:.1f}"
                    mn = "-" if c["perloc"] is None else f"{c['perloc']:.1f}"
                    line += f" {md:>10}{mn:>11} |"
                print(line)

    print(f"\nS6. SENSITIVITY -- the same four tables with the class rule swapped to "
          f"'more than {ALT:.0%} of a locus's cores are multi-copy'. The unique tables barely move;")
    print("    the problem tables are where the definition does the work.")
    for t in ("DEL", "DUP"):
        for s in ("unique", "problem"):
            table(t, s, A, B_, alt=True)


if __name__ == "__main__":
    main()
