#!/usr/bin/env python3
"""Subtract a normal's calls from a tumour's, so a germline caller can be scored on somatic truth.

Sniffles2, cuteSV and pbsv have no notion of a control sample: run on a tumour they report every
germline SV the genome carries alongside the somatic ones. Scoring that against a somatic truth set
would count thousands of correct germline calls as false positives and measure nothing about the
tool. So each is run twice, and a tumour call survives here only if the normal has no call of its own
within TOL bp of the same breakpoint. That is the same tumour-minus-normal step a somatic caller does
internally, done outside so the comparison is like for like.

Matching is on POSITION, not type: a caller that gets the locus right but the class wrong has still
seen the germline event, and letting a type mismatch push it through would manufacture a somatic call
out of a germline one.

usage: somatic_diff.py <tumour.vcf> <normal.vcf> <out.vcf> [tolerance=500]
"""
import re
import sys

BND = re.compile(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]')


def info(s, key, default=""):
    m = re.search(rf'(?:^|;){key}=([^;]*)', s)
    return m.group(1) if m else default


def breakends(line):
    """every reference locus this record touches, as (chrom, pos)"""
    f = line.rstrip('\n').split('\t')
    out = [(f[0], int(f[1]))]
    m = BND.search(f[4])
    if m:
        out.append((m.group(1), int(m.group(2))))
    else:
        end = info(f[7], 'END')
        if end.isdigit():
            out.append((f[0], int(end)))
    return out


def main():
    tum, norm, out = sys.argv[1], sys.argv[2], sys.argv[3]
    tol = int(sys.argv[4]) if len(sys.argv) > 4 else 500

    # normal breakends bucketed by kb, so a lookup only touches the three buckets a tolerance can span
    seen = {}
    for line in open(norm):
        if line[0] == '#':
            continue
        for c, p in breakends(line):
            seen.setdefault((c, p // 1000), []).append(p)

    def in_normal(c, p):
        for d in (-1, 0, 1):
            for q in seen.get((c, p // 1000 + d), ()):
                if abs(p - q) <= tol:
                    return True
        return False

    kept = total = 0
    with open(out, 'w') as o:
        for line in open(tum):
            if line[0] == '#':
                o.write(line)
                continue
            total += 1
            # one end matching the normal is enough to call the whole record germline
            if any(in_normal(c, p) for c, p in breakends(line)):
                continue
            o.write(line)
            kept += 1
    print(f"  somatik fark: {total} tumor cagrisi -> {kept} (normalde +/-{tol} bp karsiligi olmayan)")


if __name__ == '__main__':
    main()
