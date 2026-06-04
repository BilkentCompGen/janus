#!/usr/bin/env python3
# Coverage histogram by provenance class — confirms (or refutes) that TRUE_PHANTOM sigs
# are error-derived (concentrated at low coverage, near min_mult=3) vs genuine sigs
# (centered around the per-sample read depth, ~30x).

import struct
import sys
from collections import Counter

COUNTS = sys.argv[1] if len(sys.argv) > 1 else 'outputs/svp/v2_test/counts.bin'
PROV   = sys.argv[2] if len(sys.argv) > 2 else '/tmp/sig_provenance.bin'

with open(PROV, 'rb') as f:
    K1 = struct.unpack('<Q', f.read(8))[0]
    prov_raw = f.read(K1 * 16)

with open(COUNTS, 'rb') as f:
    K2 = struct.unpack('<Q', f.read(8))[0]
    sig_raw = f.read(K2 * 16)

assert K1 == K2, f'mismatch K {K1} != {K2}'
K = K1
print(f'K = {K:,}', file=sys.stderr)

# Cov bins.
def bin_cov(c):
    if c == 0: return '0'
    if c <= 4: return f'{c}'
    if c <= 9: return '5-9'
    if c <= 19: return '10-19'
    if c <= 39: return '20-39'
    return '40+'

ncov_dist = {cls: Counter() for cls in ['BOTH', 'N_only', 'T_only', 'PHANTOM']}
tcov_dist = {cls: Counter() for cls in ['BOTH', 'N_only', 'T_only', 'PHANTOM']}
maxcov_dist = {cls: Counter() for cls in ['BOTH', 'N_only', 'T_only', 'PHANTOM']}
count_by = Counter()
sum_n = {c:0 for c in ncov_dist}; sum_t = {c:0 for c in ncov_dist}

for i in range(K):
    sig, n, t = struct.unpack_from('<QII', sig_raw, i*16)
    nf, nr, tf, tr = struct.unpack_from('<IIII', prov_raw, i*16)
    in_n = (nf+nr) > 0; in_t = (tf+tr) > 0
    if   in_n and in_t: cls = 'BOTH'
    elif in_n:          cls = 'N_only'
    elif in_t:          cls = 'T_only'
    else:               cls = 'PHANTOM'
    ncov_dist[cls][bin_cov(n)] += 1
    tcov_dist[cls][bin_cov(t)] += 1
    maxcov_dist[cls][bin_cov(max(n, t))] += 1
    count_by[cls] += 1
    sum_n[cls] += n; sum_t[cls] += t

bins = ['0', '1', '2', '3', '4', '5-9', '10-19', '20-39', '40+']

print('\n== Coverage histograms by provenance class ==')
for cls in ['BOTH', 'N_only', 'T_only', 'PHANTOM']:
    tot = max(1, count_by[cls])
    mean_n = sum_n[cls] / tot
    mean_t = sum_t[cls] / tot
    print(f'\n[{cls}]  n_sigs = {count_by[cls]:,}   mean ncov = {mean_n:.1f}   mean tcov = {mean_t:.1f}')
    print(f'  bin       ncov            tcov            max(n,t)')
    for b in bins:
        n = ncov_dist[cls].get(b, 0); t = tcov_dist[cls].get(b, 0); m = maxcov_dist[cls].get(b, 0)
        print(f'  {b:>6s}    {n:>9,} ({100*n/tot:5.2f}%)   {t:>9,} ({100*t/tot:5.2f}%)   {m:>9,} ({100*m/tot:5.2f}%)')
