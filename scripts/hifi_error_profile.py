#!/usr/bin/env python3
"""Is homopolymer compression justified on THIS data, or is it folklore?

svp compresses homopolymers before hashing a core, on the stated grounds that long-read indel errors
are largely run-length slips. That claim is checkable on the reads themselves.

It cannot be checked against the reference: a germline SNP happens about once per 1000 bp and the
measured read error rate here is 0.018%, so most differences from the reference are real variants,
not errors. So the reads are compared to EACH OTHER -- at a locus with enough depth, a difference one
read carries alone is an error and a difference half or all of them carry is a variant.

Two numbers decide it: what share of errors are indels rather than substitutions, and how much more
often those indel errors land inside a homopolymer run than chance would put them. Substitutions are
the control -- they have no reason to prefer a run, so if they show the same enrichment the measure
is picking up something else.

usage: hifi_error_profile.py <normal.bam> [chrom=chr1] [ref.fa]
"""
import collections
import sys

import pysam

BAM = sys.argv[1]
REG = sys.argv[2] if len(sys.argv) > 2 else "chr1"
REF = sys.argv[3] if len(sys.argv) > 3 else "/home/datasets/references/genomes/hg38/hg38.fa"
START, END = 2_000_000, 7_000_000
MIN_DEPTH = 10
SINGLETON = 0.20            # a difference this rare among the reads covering it is an error

ref = pysam.FastaFile(REF)
bam = pysam.AlignmentFile(BAM)
RS = ref.fetch(REG, START, END).upper()


def hp_run(pos):
    """length of the homopolymer run of the reference base at `pos` (absolute coord)"""
    i = pos - START
    if i < 1 or i >= len(RS) - 1:
        return 0
    b = RS[i]
    if b not in "ACGT":
        return 0
    L, j = 1, i - 1
    while j >= 0 and RS[j] == b:
        L += 1; j -= 1
    j = i + 1
    while j < len(RS) and RS[j] == b:
        L += 1; j += 1
    return L


# pass 1: every difference, and how many reads at that locus carry it
ins = collections.Counter(); dele = collections.Counter(); sub = collections.Counter()
depth = collections.Counter()
nreads = 0
for a in bam.fetch(REG, START, END):
    if a.is_secondary or a.is_unmapped or a.is_supplementary or a.mapping_quality < 20:
        continue
    nreads += 1
    rp = a.reference_start
    qs = a.query_sequence
    if qs is None:
        continue
    qp = 0
    for op, ln in a.cigartuples:
        if op in (0, 7, 8):
            for k in range(ln):
                p = rp + k
                if START <= p < END:
                    depth[p] += 1
                    if qs[qp + k] != RS[p - START]:
                        sub[p] += 1
            rp += ln; qp += ln
        elif op == 1:
            if START <= rp < END:
                ins[(rp, ln)] += 1
            qp += ln
        elif op == 2:
            if START <= rp < END:
                dele[(rp, ln)] += 1
            rp += ln
        elif op == 4:
            qp += ln
        elif op == 3:
            rp += ln

print(f"  {nreads:,} reads, {REG}:{START:,}-{END:,}")


def errors(counter, keyed_by_pos):
    """a difference carried by a small minority of the reads that cover it"""
    out = []
    for key, n in counter.items():
        p = key[0] if keyed_by_pos else key
        d = depth.get(p, 0)
        if d >= MIN_DEPTH and n / d <= SINGLETON:
            out.append((key, n, d))
    return out


e_sub = errors(sub, False)
e_ins = errors(ins, True)
e_del = errors(dele, True)
n_sub, n_ins, n_del = len(e_sub), len(e_ins), len(e_del)
tot = n_sub + n_ins + n_del
print(f"\n  ERRORS (a difference seen in <={100*SINGLETON:.0f}% of the reads over that locus):")
print(f"    substitution {n_sub:7,d}  ({100*n_sub/tot:.1f}%)")
print(f"    insertion    {n_ins:7,d}  ({100*n_ins/tot:.1f}%)")
print(f"    deletion     {n_del:7,d}  ({100*n_del/tot:.1f}%)")
print(f"    indel total  {n_ins+n_del:7,d}  ({100*(n_ins+n_del)/tot:.1f}%)")

# background: how much of this window is inside a homopolymer of each length
bg = collections.Counter()
step = 37                                            # sample the window, not every base
for p in range(START, END, step):
    bg[min(hp_run(p), 8)] += 1
bgtot = sum(bg.values())

for lab, errs in (("indel", e_ins + e_del), ("substitution", [((k, 0), n, d) for k, n, d in e_sub])):
    h = collections.Counter()
    for key, _n, _d in errs:
        h[min(hp_run(key[0]), 8)] += 1
    t = sum(h.values())
    if not t:
        continue
    print(f"\n  homopolymer context of the {lab} errors (reference run length):")
    print(f"    {'run':>5s} {'of errors':>10s} {'of genome':>11s} {'enrichment':>12s}")
    for L in range(1, 9):
        he, be = 100 * h[L] / t, 100 * bg[L] / bgtot
        print(f"    {L if L < 8 else '8+':>5} {he:9.1f}% {be:10.1f}% {he/be if be else 0:11.1f}x")
    ge3 = sum(h[L] for L in range(3, 9)) / t
    bge3 = sum(bg[L] for L in range(3, 9)) / bgtot
    print(f"    run>=3: {100*ge3:.1f}% of errors, {100*bge3:.1f}% of the genome  -> {ge3/bge3:.1f}x")
