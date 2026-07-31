#!/usr/bin/env python3
"""Are the simulated reads' indel errors homopolymer-length slips, or are they uniform?

The assemblies come out at ~130 indels/100 kb with indels outnumbering mismatches ~30:1, which LOOKS
like the classic long-read homopolymer signature. This measures it instead of inferring it: align a
sample of reads to the reference, take every indel from the CIGAR, and ask whether it sits inside a
homopolymer run of the reference.

The number only means something against a null: a random indel would land in a run of length >= L at
whatever rate the reference itself is made of such runs. So the reference's own run-length composition
is measured too, and the two are reported side by side. Enrichment = observed / expected.

usage: homopolymer_errors.py <reads.fastq> <ref.fasta> <n_reads> <threads> <outdir>
"""
import collections
import os
import subprocess
import sys


def sample_reads(src, dst, n):
    """First n records of a FASTQ. They are already in shard order, i.e. not position-sorted."""
    with open(src) as fi, open(dst, "w") as fo:
        for i, line in enumerate(fi):
            if i >= n * 4:
                break
            fo.write(line)


def load_ref(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line[0] == ">":
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def run_len_at(seq, pos):
    """Length of the homopolymer run of seq that covers pos (1 if the base differs from both sides)."""
    if pos < 0 or pos >= len(seq):
        return 0
    b = seq[pos]
    if b not in "ACGT":
        return 0
    i = pos
    while i > 0 and seq[i - 1] == b:
        i -= 1
    j = pos
    while j + 1 < len(seq) and seq[j + 1] == b:
        j += 1
    return j - i + 1


def ref_run_composition(seqs, sample_every=97):
    """What fraction of reference BASES sit in a run of length >= L. This is the null: a uniformly
    placed indel lands in such a run at exactly this rate."""
    hist = collections.Counter()
    total = 0
    for s in seqs.values():
        i, n = 0, len(s)
        while i < n:
            b = s[i]
            j = i
            while j + 1 < n and s[j + 1] == b:
                j += 1
            L = j - i + 1
            if b in "ACGT":
                hist[L] += L          # every base of the run is in a run of that length
                total += L
            i = j + 1
    return hist, total


def main():
    reads, ref, n, threads, outdir = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
    os.makedirs(outdir, exist_ok=True)
    sub = os.path.join(outdir, "sample.fastq")
    paf = os.path.join(outdir, "sample.paf")

    if not os.path.exists(sub):
        print(f"sampling {n:,} reads ...", flush=True)
        sample_reads(reads, sub, n)
    if not os.path.exists(paf):
        print("aligning with minimap2 ...", flush=True)
        with open(paf, "w") as fo:
            subprocess.run(["minimap2", "-x", "map-hifi", "-t", threads, "-c", "--cs=long",
                            "--secondary=no", ref, sub], stdout=fo, check=True)

    print("loading reference ...", flush=True)
    seqs = load_ref(ref)

    # --- walk the cs tag: =MATCH  *xy(sub)  +ins  -del ---
    ins = collections.Counter()      # run length at the insertion point -> count
    dele = collections.Counter()
    sub_hist = collections.Counter()
    ins_len, del_len = collections.Counter(), collections.Counter()
    ins_expand = collections.Counter()   # insertions that provably lengthen a homopolymer run
    aligned = 0
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        tname, tstart = f[5], int(f[7])
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None or tname not in seqs:
            continue
        aligned += 1
        s = seqs[tname]
        pos = tstart
        i = 0
        while i < len(cs):
            op = cs[i]
            j = i + 1
            while j < len(cs) and cs[j] not in "=*+-":
                j += 1
            arg = cs[i + 1:j]
            if op == "=":
                pos += len(arg)
            elif op == "*":
                sub_hist[run_len_at(s, pos)] += 1
                pos += 1
            elif op == "+":                      # inserted bases, reference does not advance
                ins[run_len_at(s, pos - 1)] += 1
                ins_len[len(arg)] += 1
                # An insertion's "run length at the preceding base" is a weak proxy: the alignment is
                # free to place it on either side of a run. The direct test is whether the inserted
                # bases are all one base that CONTINUES the run on either side -- i.e. the read
                # lengthened a homopolymer. That is unambiguous however the aligner broke the tie.
                a = arg.upper()
                if a and len(set(a)) == 1 and (
                        (pos > 0 and s[pos - 1] == a[0]) or (pos < len(s) and s[pos] == a[0])):
                    ins_expand[max(run_len_at(s, pos - 1), run_len_at(s, pos))] += 1
            elif op == "-":                      # deleted reference bases
                dele[run_len_at(s, pos)] += 1
                del_len[len(arg)] += 1
                pos += len(arg)
            i = j

    print("measuring reference run composition ...", flush=True)
    ref_hist, ref_total = ref_run_composition(seqs)

    def frac_ge(hist, L):
        tot = sum(hist.values())
        return sum(v for k, v in hist.items() if k >= L) / tot if tot else 0.0

    n_ins, n_del, n_sub = sum(ins.values()), sum(dele.values()), sum(sub_hist.values())
    print(f"\naligned reads: {aligned:,}")
    print(f"substitutions {n_sub:,} | insertions {n_ins:,} | deletions {n_del:,}"
          f" | indel:sub = {(n_ins + n_del) / max(n_sub, 1):.1f}:1")

    print(f"\n{'run len':>9}{'reference':>12}{'insertions':>12}{'deletions':>12}"
          f"{'subs':>10}{'ins enrich':>12}{'del enrich':>12}")
    for L in (1, 2, 3, 4, 5, 6, 8, 10):
        r = frac_ge(ref_hist, L)
        fi_, fd, fs = frac_ge(ins, L), frac_ge(dele, L), frac_ge(sub_hist, L)
        print(f"{'>=' + str(L):>9}{r:>12.2%}{fi_:>12.2%}{fd:>12.2%}{fs:>10.2%}"
              f"{fi_ / r if r else 0:>12.1f}x{fd / r if r else 0:>11.1f}x")

    n_exp = sum(ins_expand.values())
    print(f"\nhomopolymer-EXPANDING insertions (inserted base continues a run): {n_exp:,} / {n_ins:,}"
          f" = {n_exp / max(n_ins, 1):.1%} of all insertions")

    print(f"\nindel length distribution (bp: insertions / deletions)")
    for L in sorted(set(ins_len) | set(del_len))[:10]:
        print(f"  {L:>3}: {ins_len[L]:>10,} / {del_len[L]:>10,}")
    big_i = sum(v for k, v in ins_len.items() if k >= 50)
    big_d = sum(v for k, v in del_len.items() if k >= 50)
    print(f"  >=50bp: {big_i:,} / {big_d:,}   <- these are the ones an SV caller would see")


if __name__ == "__main__":
    main()
