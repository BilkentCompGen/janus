#!/usr/bin/env python3
# For each truth-bearing color=2 super-node, classify its bubble shape AND record the
# super-node's chain_len + block_len (a proxy for SV "size" of the T-allele). Are the
# detour cases concentrated at small block_len (small DELs/SNPs/sub-LCP-window) or
# spread across all sizes?

import struct
import sys
from collections import Counter, defaultdict

MODEA  = sys.argv[1] if len(sys.argv) > 1 else 'outputs/svp/v2_test/modea.bin'
COUNTS = sys.argv[2] if len(sys.argv) > 2 else 'outputs/svp/v2_test/counts.bin'
PROV   = sys.argv[3] if len(sys.argv) > 3 else '/tmp/sig_provenance.bin'

with open(PROV, 'rb') as f:
    K = struct.unpack('<Q', f.read(8))[0]
    prov_raw = f.read(K * 16)
sig_class = bytearray(K)
for i in range(K):
    nf, nr, tf, tr = struct.unpack_from('<IIII', prov_raw, i*16)
    n = (nf+nr) > 0; t = (tf+tr) > 0
    if n and t: sig_class[i] = 0
    elif n: sig_class[i] = 1
    elif t: sig_class[i] = 2
    else:   sig_class[i] = 3

with open(MODEA, 'rb') as f:
    data = f.read()
off = 0
S = struct.unpack_from('<Q', data, off)[0]; off += 8
sn_color = bytearray(S); sn_chain_len = [0]*S; sn_block_len = [0]*S
for i in range(S):
    cl, bl, rl, c = struct.unpack_from('<IIIB', data, off)
    sn_color[i] = c; sn_chain_len[i] = cl; sn_block_len[i] = bl
    off += 16
tc = struct.unpack_from('<Q', data, off)[0]; off += 8
chain_ids_raw = data[off:off+tc*4]; off += tc*4
to = struct.unpack_from('<Q', data, off)[0]; off += 8 + to*2
tb = struct.unpack_from('<Q', data, off)[0]; off += 8 + tb
tr = struct.unpack_from('<Q', data, off)[0]; off += 8 + tr*12
E = struct.unpack_from('<Q', data, off)[0]; off += 8

node2super = [0xFFFFFFFF] * K
super_first = [0]*S; cur = 0
for s in range(S):
    super_first[s] = cur; cur += sn_chain_len[s]

# Find truth-bearing color=2 super-nodes (≥1 sig with sig_class==2).
def super_first_truth(sid):
    p = super_first[sid]
    for j in range(sn_chain_len[sid]):
        nid = struct.unpack_from('<I', chain_ids_raw, (p+j)*4)[0]
        if sig_class[nid] == 2: return True
    return False

truth_super = [s for s in range(S) if sn_color[s] == 2 and super_first_truth(s)]
print(f'truth-bearing color=2 super-nodes = {len(truth_super):,}', file=sys.stderr)

# Build adjacency
out_adj = [[] for _ in range(S)]; in_adj = [[] for _ in range(S)]
for k in range(E):
    u, v = struct.unpack_from('<II', data, off + k*8)
    out_adj[u].append(v); in_adj[v].append(u)

# For each truth-bearing color=2 super-node, find brackets and classify.
NODE_BUDGET = 200
def find_brackets(sn):
    lefts = set(); rights = set()
    stack = [sn]; visited = {sn}; budget = NODE_BUDGET
    while stack and budget > 0:
        u = stack.pop(); budget -= 1
        for v in in_adj[u]:
            if sn_color[v] == 3: lefts.add(v)
            elif v not in visited: visited.add(v); stack.append(v)
    stack = [sn]; visited = {sn}; budget = NODE_BUDGET
    while stack and budget > 0:
        u = stack.pop(); budget -= 1
        for v in out_adj[u]:
            if sn_color[v] == 3: rights.add(v)
            elif v not in visited: visited.add(v); stack.append(v)
    return lefts, rights

def has_color1_path(s1, s2, max_hops=15):
    if s1 == s2: return False
    stack = [(s1, 0)]; visited = {s1}; budget = NODE_BUDGET
    while stack and budget > 0:
        u, d = stack.pop(); budget -= 1
        if d >= max_hops: continue
        for v in out_adj[u]:
            if v == s2: return True
            if sn_color[v] == 3 or v in visited: continue
            if sn_color[v] != 1: continue
            visited.add(v); stack.append((v, d+1))
    return False

def has_color3_path(s1, s2, max_hops=8):
    if s1 == s2: return False
    stack = [(s1, 0)]; visited = {s1}; budget = NODE_BUDGET
    while stack and budget > 0:
        u, d = stack.pop(); budget -= 1
        if d >= max_hops: continue
        for v in out_adj[u]:
            if v == s2 and sn_color[v] == 3: return True
            if v in visited or sn_color[v] != 3: continue
            visited.add(v); stack.append((v, d+1))
    return False

import random
random.seed(42)
sample = truth_super if len(truth_super) < 8000 else random.sample(truth_super, 8000)

def bl_bin(bp):
    if bp < 50: return '<50'
    if bp < 100: return '50-99'
    if bp < 200: return '100-199'
    if bp < 500: return '200-499'
    if bp < 5000: return '500-4999'
    if bp < 50000: return '5K-50K'
    return '50K+'

# Tabulate by bubble shape × block_len bin
table = defaultdict(Counter)
class_count = Counter()
total_by_blbin = Counter()
for sn in sample:
    lefts, rights = find_brackets(sn)
    if not lefts or not rights:
        cls = 'no_brackets'
    else:
        found_std = False; found_detour = False; found_orphan = False
        for L in lefts:
            for R in rights:
                if has_color1_path(L, R): found_std = True
                elif has_color3_path(L, R): found_detour = True
                else: found_orphan = True
        if found_std: cls = 'std_both'
        elif found_detour: cls = 'detour'
        else: cls = 'orphan'
    bl = sn_block_len[sn]; cl_chain = sn_chain_len[sn]
    bin_key = bl_bin(bl)
    table[cls][bin_key] += 1
    class_count[cls] += 1
    total_by_blbin[bin_key] += 1

print('\n== Truth-bearing color=2 super-node bubble shape × block_len ==')
print(f'  sample size = {len(sample):,}')
print(f'  class totals: {dict(class_count)}')
order = ['<50', '50-99', '100-199', '200-499', '500-4999', '5K-50K', '50K+']
print(f'\n{"block_len bin":<12s}  {"total":>7s}  {"std_both":>9s}  {"detour":>7s}  {"orphan":>7s}  {"no_br":>5s}')
for b in order:
    tot = total_by_blbin[b]
    if tot == 0: continue
    std = table['std_both'][b]; det = table['detour'][b]; orph = table['orphan'][b]; nb = table['no_brackets'][b]
    print(f'{b:<12s}  {tot:>7,}  {std:>5,} ({100*std/tot:4.1f}%)  {det:>3,} ({100*det/tot:4.1f}%)  {orph:>3,} ({100*orph/tot:4.1f}%)  {nb:>5,}')

# Also report mean block_len per class
print('\n== mean block_len (T-allele super-node) per class ==')
total_bp = defaultdict(int); count = defaultdict(int)
for sn in sample:
    lefts, rights = find_brackets(sn)
    if not lefts or not rights:
        cls = 'no_brackets'
    else:
        found_std = False; found_detour = False; found_orphan = False
        for L in lefts:
            for R in rights:
                if has_color1_path(L, R): found_std = True
                elif has_color3_path(L, R): found_detour = True
                else: found_orphan = True
        if found_std: cls = 'std_both'
        elif found_detour: cls = 'detour'
        else: cls = 'orphan'
    total_bp[cls] += sn_block_len[sn]
    count[cls] += 1
for cls in ['std_both', 'detour', 'orphan', 'no_brackets']:
    if count[cls] == 0: continue
    print(f'  {cls:15s} mean block_len = {total_bp[cls]/count[cls]:>7.1f} bp   (n = {count[cls]})')
