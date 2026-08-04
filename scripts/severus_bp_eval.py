#!/usr/bin/env python3
"""Score a Severus VCF at the JUNCTION level against the simulator's own novel adjacencies.

Truvari scores variants; this scores breakpoints, which is what the graph actually produces. The
truth is insilicoSV's `adjacencies.bed` -- the pairs of reference loci the tumour puts next to each
other and the normal does not -- so a deletion contributes ONE junction, an inversion TWO, and a
complex event as many as its grammar creates. Insertions have no adjacency row (their sequence has no
reference interval), so they are matched separately on position and length.

A call matches a truth junction when BOTH of its breakpoints land within TOL bp of the truth's two
breakpoints, on the right contigs, in either pairing order. Nothing about orientation is required:
getting the strand wrong but both coordinates right is still a located breakpoint, and conflating the
two would hide which one failed.

usage: severus_bp_eval.py <severus.vcf> <truth.vcf> <adjacencies.bed> [tolerance=100]
"""
import collections
import re
import sys

MINSV = 50          # Severus's own floor; anything under it is not in its output by construction
BANDS = ((50, 1000), (1000, 10000), (10000, 100000), (100000, 10 ** 12))


def info(s, key, default=""):
    m = re.search(r"(?:^|;)" + key + r"=([^;]*)", s)
    return m.group(1) if m else default


def load_truth(vcf):
    """svid -> (svtype, svlen, chrom, pos, end) for every SV at or above Severus's size floor."""
    sv, ins = {}, []
    for line in open(vcf):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        t, sl = info(f[7], "SVTYPE"), info(f[7], "SVLEN")
        if not t:
            continue
        L = int(sl or 0)
        if L < MINSV:
            continue
        rec = (t, L, f[0], int(f[1]), int(info(f[7], "END") or f[1]))
        sv[info(f[7], "SVID") or f[2]] = rec
        if t == "INS":
            ins.append(rec)
    return sv, ins


def load_adj(bed, sv):
    """One entry per truth junction of an SV that is big enough for Severus to report."""
    out = []
    for line in open(bed):
        c = line.rstrip("\n").split("\t")
        if len(c) < 10 or c[9] not in sv:
            continue
        out.append((c[0], int(c[1]) + 1, c[3], int(c[4]) + 1, c[9]))
    return out


BND = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def load_calls(vcf):
    """Every Severus record as a junction (chrom1, pos1, chrom2, pos2), plus insertions apart."""
    jn, ins, seen_mate = [], [], set()
    for line in open(vcf):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        if f[6] not in ("PASS", "."):
            continue
        t = info(f[7], "SVTYPE")
        chrom, pos = f[0], int(f[1])
        if t == "INS" and not BND.search(f[4]):
            ins.append((chrom, pos, int(info(f[7], "SVLEN") or 0), f[2]))
            continue
        # A bracket ALT names the other side outright, so take it whenever it is there rather than
        # only when SVTYPE says BND -- the COLO829 truth set writes DEL/DUP/INV/INS that way and
        # gives no END at all, and reading those as untyped junctions is what makes them scoreable.
        m = BND.search(f[4])
        if m:
            mate = info(f[7], "MATE_ID") or info(f[7], "MATEID")   # Severus spells it MATE_ID
            if f[2] in seen_mate:                      # a breakend pair is ONE junction
                continue
            if mate:
                seen_mate.add(mate)
            jn.append((chrom, pos, m.group(1), int(m.group(2)), t, f[2]))
            continue
        end = int(info(f[7], "END") or 0)
        if end:
            jn.append((chrom, pos, chrom, end, t, f[2]))
    return jn, ins


def near(a, b, tol):
    return abs(a - b) <= tol


def match(call, truth, tol):
    c1, p1, c2, p2 = call[:4]
    t1, q1, t2, q2, _ = truth
    return ((c1 == t1 and c2 == t2 and near(p1, q1, tol) and near(p2, q2, tol)) or
            (c1 == t2 and c2 == t1 and near(p1, q2, tol) and near(p2, q1, tol)))


def band(n):
    for lo, hi in BANDS:
        if lo <= n < hi:
            return f"{lo}-{hi if hi < 10**12 else 'inf'}"
    return "?"


def main():
    svcf, tvcf, bed = sys.argv[1:4]
    tol = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    sv, truth_ins = load_truth(tvcf)
    adj = load_adj(bed, sv)
    calls, call_ins = load_calls(svcf)

    hit_truth, hit_call = set(), set()
    for i, a in enumerate(adj):
        for j, c in enumerate(calls):
            if match(c, a, tol):
                hit_truth.add(i); hit_call.add(j)
    # insertions: position and length, not a pair of breakpoints
    ins_hit, ins_call_hit = set(), set()
    for i, (t, L, ch, p, _e) in enumerate(truth_ins):
        for j, (cc, cp, cl, _id) in enumerate(call_ins):
            if cc == ch and near(cp, p, tol) and (cl == 0 or abs(cl - L) <= max(100, 0.5 * L)):
                ins_hit.add(i); ins_call_hit.add(j)

    print(f"tolerance +/-{tol} bp   truth junctions {len(adj)} (from {len(sv)} SVs >= {MINSV} bp)"
          f"   truth insertions {len(truth_ins)}")
    print(f"severus calls: {len(calls)} junctions + {len(call_ins)} insertions\n")

    per = collections.defaultdict(lambda: [0, 0])
    for i, a in enumerate(adj):
        t, L = sv[a[4]][0], sv[a[4]][1]
        for k in ((t, "all"), (t, band(L))):
            per[k][0] += 1
            per[k][1] += int(i in hit_truth)
    for i, r in enumerate(truth_ins):
        for k in (("INS", "all"), ("INS", band(r[1]))):
            per[k][0] += 1
            per[k][1] += int(i in ins_hit)

    print(f"  {'type':<6s} {'size band':<14s} {'truth jn':>9s} {'recovered':>10s} {'recall':>8s}")
    for (t, b) in sorted(per, key=lambda k: (k[0], k[1] != "all", k[1])):
        n, h = per[(t, b)]
        print(f"  {t:<6s} {b:<14s} {n:>9d} {h:>10d} {100*h/n:>7.1f}%")

    fp_j, fp_i = len(calls) - len(hit_call), len(call_ins) - len(ins_call_hit)
    tot_c, tot_h = len(calls) + len(call_ins), len(hit_call) + len(ins_call_hit)
    print(f"\n  calls matching a truth junction : {tot_h}/{tot_c} = {100*tot_h/tot_c if tot_c else 0:.1f}% precision")
    print(f"  unmatched: {fp_j} junction calls + {fp_i} insertion calls")
    un = collections.Counter(c[4] for j, c in enumerate(calls) if j not in hit_call)
    print(f"  unmatched junction calls by SVTYPE: {dict(un)}")


if __name__ == "__main__":
    main()
