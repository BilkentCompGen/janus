#!/usr/bin/env python3
import sys, collections, argparse

ap = argparse.ArgumentParser()
ap.add_argument('tumor'); ap.add_argument('normal')
ap.add_argument('--somatic-tag', action='store_true')
ap.add_argument('--min-len', type=int, default=50)
ap.add_argument('--dist', type=int, default=1000)
a = ap.parse_args()

def info(f):
    d = {}
    for kv in f[7].split(';'):
        k, _, v = kv.partition('=')
        d[k] = v
    return d

def svtype(f, d):
    if 'SVTYPE' in d: return d['SVTYPE']
    return 'INS' if len(f[4]) > len(f[3]) else 'DEL'

def is_sv(f, d):
    if f[4].startswith('<') or '[' in f[4] or ']' in f[4]: return True
    if 'SVLEN' in d:
        try: return abs(int(d['SVLEN'].split(',')[0])) >= a.min_len
        except ValueError: pass
    return abs(len(f[4]) - len(f[3])) >= a.min_len

normal = collections.defaultdict(list)
for l in open(a.normal):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t'); d = info(f)
    if is_sv(f, d): normal[(f[0], svtype(f, d))].append(int(f[1]))

for l in open(a.tumor):
    if l[0] == '#': sys.stdout.write(l); continue
    f = l.rstrip('\n').split('\t'); d = info(f)
    if f[6] not in ('PASS', '.') or not is_sv(f, d): continue
    if a.somatic_tag and 'SOMATIC' not in d: continue
    if any(abs(p - int(f[1])) <= a.dist for p in normal[(f[0], svtype(f, d))]): continue
    sys.stdout.write(l)
