#!/usr/bin/env python3
"""Verify outputs/svp/gt/*.jsonl against the simulator's own tumour genome.

For each sampled truth record we rebuild the ALT allele exactly the way tools/gt_cores does
(reference flank + variant + reference flank) and look for the 120 bp centred on the junction in
insilicoSV's out.fa. If the splice point or the reference span were off by even one base, the probe
would not be there. This is what caught the INS insertion point being POS instead of POS-1.

usage: verify_gt_cores.py [hard|complex] [n_per_type=40]
"""
import bisect
import collections
import json
import random
import sys

W = 300           # flank each side of the junction
PROBE = 60        # half-width of the probe taken from the junction
RC = str.maketrans("ACGTN", "TGCAN")


def faidx(path):
    ix = {}
    for ln in open(path + ".fai"):
        f = ln.split()
        ix[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return ix


def fetch(path, ix, chrom, s0, n):
    length, off, lb, lw = ix[chrom]
    s0 = max(0, s0)
    n = min(n, length - s0)
    fs = off + (s0 // lb) * lw + (s0 % lb)
    nb = n + (n // lb + 2) * (lw - lb)
    with open(path, "rb") as f:
        f.seek(fs)
        raw = f.read(nb)
    return raw.replace(b"\n", b"").replace(b"\r", b"")[:n].upper().decode()


def insseqs(vcf):
    out = {}
    for ln in open(vcf):
        if ln[0] == "#":
            continue
        f = ln.split("\t")
        if "INSSEQ=" in f[7]:
            out[f[2]] = f[7].split("INSSEQ=")[1].split(";")[0]
    return out


def shift_table(recs):
    """Cumulative ref->tumour offset per contig: an INS/DUP pushes later bases right, a DEL left."""
    by = collections.defaultdict(list)
    for r in recs:
        by[r["chrom"]].append(r)
    cum = {}
    for c, v in by.items():
        v.sort(key=lambda r: r["pos"])
        tot, arr = 0, []
        for r in v:
            arr.append((r["pos"], tot))
            t = r["svtype"]
            if t in ("INS", "DUP"):
                tot += abs(r["svlen"])
            elif t == "DEL":
                tot -= abs(r["svlen"])
        cum[c] = arr
    return cum


def alt_allele(r, t, ref, rix, ins, s0, e0):
    left = fetch(ref, rix, r["chrom"], s0 - W, W)
    right = fetch(ref, rix, r["chrom"], e0, W)
    seg = fetch(ref, rix, r["chrom"], s0, e0 - s0)
    if t == "DEL":
        return left + right
    if t == "INS":
        return left + ins.get(r["id"], "") + right
    if t == "INV":
        return left + seg.translate(RC)[::-1] + right
    if t == "DUP":
        return left + seg + seg + right
    if t in ("SNV", "INDEL"):
        return left + r.get("_alt", "") + right
    return None


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "hard"
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    ref, tum = "simulation/ref.fasta", f"simulation/truth/{name}/out.fa"
    vcf = f"simulation/truth/{name}/out.vcf"
    rix, tix = faidx(ref), faidx(tum)
    recs = [json.loads(l) for l in open(f"outputs/svp/gt/{name}.jsonl")]

    alts = {}
    for ln in open(vcf):
        if ln[0] == "#":
            continue
        f = ln.split("\t")
        if not f[4].startswith("<"):
            alts[f[2]] = f[4].upper()
    for r in recs:
        r["_alt"] = alts.get(r["id"], "")

    ins = insseqs(vcf)
    cum = shift_table(recs)

    def off_at(c, pos):
        arr = cum[c]
        i = bisect.bisect_left(arr, (pos, -10**9))
        return arr[i - 1][1] if i > 0 else 0

    random.seed(11)
    byt = collections.defaultdict(list)
    for r in recs:
        byt[r["svtype"]].append(r)

    # Two things the naive probe cannot see, and both are properties of the TRUTH, not of gt_cores:
    #  - the offset model ignores translocations/dispersed dups, so the predicted tumour window can be
    #    wrong; fall back to the whole contig.
    #  - a locus inside a larger inversion appears reverse-complemented; try the RC probe too.
    print(f"=== {name}: ALT junction probes looked up in {tum}")
    print(f"{'type':<8}{'tested':>8}{'window':>8}{'contig':>8}{'rc':>6}{'ABSENT':>8}")
    cache = {}
    tot = collections.Counter()
    for t in ("SNV", "INDEL", "INS", "DEL", "DUP", "INV"):
        if not byt[t]:
            continue
        c = collections.Counter()
        for r in random.sample(byt[t], min(per, len(byt[t]))):
            s0, e0 = r["span"]
            tc = r["chrom"] + "_hapA"
            if tc not in tix:
                continue
            alt = alt_allele(r, t, ref, rix, ins, s0, e0)
            if alt is None:
                continue
            probe = alt[W - PROBE:W + PROBE]
            c["tested"] += 1
            tpos = s0 + off_at(r["chrom"], r["pos"])
            if probe in fetch(tum, tix, tc, max(0, tpos - 30000), 60000):
                c["window"] += 1
                continue
            if tc not in cache:
                cache.clear()
                cache[tc] = fetch(tum, tix, tc, 0, tix[tc][0])
            if probe in cache[tc]:
                c["contig"] += 1
            elif probe.translate(RC)[::-1] in cache[tc]:
                c["rc"] += 1
            else:
                c["absent"] += 1
                print(f"         ABSENT {r['id']} {t} {r['chrom']}:{r['pos']} svlen={r['svlen']}")
        tot.update(c)
        print(f"{t:<8}{c['tested']:>8}{c['window']:>8}{c['contig']:>8}{c['rc']:>6}{c['absent']:>8}")
    found = tot["window"] + tot["contig"] + tot["rc"]
    print(f"{'TOTAL':<8}{tot['tested']:>8}{tot['window']:>8}{tot['contig']:>8}{tot['rc']:>6}{tot['absent']:>8}"
          f"   -> {found}/{tot['tested']} junctions confirmed in the tumour genome")


if __name__ == "__main__":
    main()
