#!/usr/bin/env python3
# For the 17K N_only_no_T and 18K T_only_no_N pairs, check whether a parallel
# color=3-only path exists between S1 and S2 (the "tumor/normal backbone" detour
# that bypasses the allele). If yes, the pair is actually a DEL/INS bubble with
# the opposite allele implicit in the backbone — not an "orphan".

import struct
import sys
from collections import defaultdict, deque, Counter

MODEA = sys.argv[1] if len(sys.argv) > 1 else 'outputs/svp/v2_test/modea.bin'

with open(MODEA, 'rb') as f:
    data = f.read()

off = 0
S = struct.unpack_from('<Q', data, off)[0]; off += 8
color = bytearray(S)
for i in range(S):
    cl, bl, rl, c = struct.unpack_from('<IIIB', data, off)
    color[i] = c
    off += 16
tc = struct.unpack_from('<Q', data, off)[0]; off += 8 + tc*4
to = struct.unpack_from('<Q', data, off)[0]; off += 8 + to*2
tb = struct.unpack_from('<Q', data, off)[0]; off += 8 + tb
tr = struct.unpack_from('<Q', data, off)[0]; off += 8 + tr*12
E = struct.unpack_from('<Q', data, off)[0]; off += 8

out_adj = [[] for _ in range(S)]
for k in range(E):
    u, v = struct.unpack_from('<II', data, off + k*8)
    out_adj[u].append(v)
print(f'loaded S={S:,} E={E:,}', file=sys.stderr)

# Re-derive the N_only_no_T and T_only_no_N pair list (cheap to recompute).
NODE_BUDGET = 500
n_only_pairs = []
t_only_pairs = []
std_both_pairs = []
for s1 in range(S):
    if color[s1] != 3: continue
    closures = defaultdict(lambda: {'c1':0,'c2':0,'direct':False})
    for x in out_adj[s1]:
        if color[x] == 3:
            closures[x]['direct'] = True
            continue
        branch_c = color[x]
        visited = {x}
        stack = [x]
        closed = set()
        budget = NODE_BUDGET
        while stack and budget > 0:
            n = stack.pop(); budget -= 1
            for nx in out_adj[n]:
                if color[nx] == 3:
                    closed.add(nx)
                else:
                    if nx in visited: continue
                    visited.add(nx); stack.append(nx)
        for s2 in closed:
            if branch_c == 1: closures[s2]['c1'] += 1
            else:             closures[s2]['c2'] += 1
    for s2, v in closures.items():
        has1 = v['c1'] > 0; has2 = v['c2'] > 0
        if has1 and has2: std_both_pairs.append((s1, s2))
        elif has1 and not has2 and not v['direct']: n_only_pairs.append((s1, s2))
        elif has2 and not has1 and not v['direct']: t_only_pairs.append((s1, s2))

print(f'N_only_no_T pairs = {len(n_only_pairs):,}', file=sys.stderr)
print(f'T_only_no_N pairs = {len(t_only_pairs):,}', file=sys.stderr)
print(f'std_both pairs    = {len(std_both_pairs):,}', file=sys.stderr)

# For each pair, BFS from S1 through color==3 only and see if S2 is reachable.
def color3_path_exists(s1, s2, max_hops=8):
    if s1 == s2: return False
    bfs = deque([(s1, 0)])
    visited = {s1}
    while bfs:
        n, d = bfs.popleft()
        if d >= max_hops: continue
        for nx in out_adj[n]:
            if nx == s2 and color[nx] == 3:
                return True
            if nx in visited: continue
            if color[nx] != 3: continue
            visited.add(nx)
            bfs.append((nx, d+1))
    return False

def color_path_exists(s1, s2, allowed, max_hops=8):
    if s1 == s2: return False
    bfs = deque([(s1, 0)])
    visited = {s1}
    while bfs:
        n, d = bfs.popleft()
        if d >= max_hops: continue
        for nx in out_adj[n]:
            if nx == s2:
                return True
            if nx in visited: continue
            if color[nx] not in allowed: continue
            visited.add(nx)
            bfs.append((nx, d+1))
    return False

# Sample analysis: check the first 5000 of each
import random
random.seed(42)
def analyze(pairs, label):
    sample = pairs if len(pairs) < 5000 else random.sample(pairs, 5000)
    has_c3 = 0
    no_c3 = 0
    for s1, s2 in sample:
        if color3_path_exists(s1, s2):
            has_c3 += 1
        else:
            no_c3 += 1
    print(f'\n[{label}] sample {len(sample):,}: has color=3 path = {has_c3:,} ({100*has_c3/len(sample):.1f}%), no color=3 path = {no_c3:,} ({100*no_c3/len(sample):.1f}%)')

analyze(n_only_pairs, 'N_only_no_T → check color=3 bypass (DEL-like?)')
analyze(t_only_pairs, 'T_only_no_N → check color=3 bypass (INS-like?)')
analyze(std_both_pairs, 'std_both → check if also has color=3 bypass (variant + germline polymorphism?)')
