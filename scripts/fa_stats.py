#!/usr/bin/env python3
"""Contig count / total length / N50 / NG50 of a FASTA, one line per file.

usage: fa_stats.py <genome_size> <a.fa> [b.fa ...]
"""
import sys


def lengths(path):
    out, n = [], 0
    for line in open(path):
        if line[0] == ">":
            if n:
                out.append(n)
            n = 0
        else:
            n += len(line.strip())
    if n:
        out.append(n)
    return sorted(out, reverse=True)


def nx(ls, target):
    h = 0
    for v in ls:
        h += v
        if h >= target:
            return v
    return 0


def main():
    g = int(sys.argv[1])
    print(f"{'file':<34}{'contigs':>9}{'total(Mb)':>11}{'N50(Mb)':>9}{'NG50(Mb)':>10}{'max(Mb)':>9}")
    for p in sys.argv[2:]:
        ls = lengths(p)
        if not ls:
            print(f"{p.split('/')[-1]:<34}{'EMPTY':>9}")
            continue
        tot = sum(ls)
        print(f"{p.split('/')[-1]:<34}{len(ls):>9}{tot/1e6:>11.2f}"
              f"{nx(ls, tot/2)/1e6:>9.2f}{nx(ls, g/2)/1e6:>10.2f}{ls[0]/1e6:>9.2f}")


if __name__ == "__main__":
    main()
