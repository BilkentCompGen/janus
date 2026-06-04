#!/usr/bin/env python3
# Follow-up: dump structure of the cycles, multi-path super-bubbles, and orphan-allele
# (N-only / T-only without direct edge) cases.

import struct
import sys
from collections import defaultdict, Counter, deque

MODEA = sys.argv[1] if len(sys.argv) > 1 else 'outputs/svp/v2_test/modea.bin'

with open(MODEA, 'rb') as f:
    data = f.read()

off = 0
S = struct.unpack_from('<Q', data, off)[0]; off += 8
color = bytearray(S)
chain_len = [0]*S
block_len = [0]*S
for i in range(S):
    cl, bl, rl, c = struct.unpack_from('<IIIB', data, off)
    color[i] = c
    chain_len[i] = cl
    block_len[i] = bl
    off += 16

tc = struct.unpack_from('<Q', data, off)[0]; off += 8 + tc*4
to = struct.unpack_from('<Q', data, off)[0]; off += 8 + to*2
tb = struct.unpack_from('<Q', data, off)[0]; off += 8 + tb
tr = struct.unpack_from('<Q', data, off)[0]; off += 8 + tr*12
E = struct.unpack_from('<Q', data, off)[0]; off += 8

out_adj = [[] for _ in range(S)]
in_adj = [[] for _ in range(S)]
for k in range(E):
    u, v = struct.unpack_from('<II', data, off + k*8)
    out_adj[u].append(v)
    in_adj[v].append(u)

def fmt(n):
    return f'sn{n}(c={color[n]},chain={chain_len[n]},bp={block_len[n]})'

def neighborhood(node, hops=2):
    out = [f'  {fmt(node)}']
    out.append(f'    out:')
    for v in out_adj[node][:10]:
        out.append(f'      → {fmt(v)}')
        if hops > 1:
            for vv in out_adj[v][:5]:
                out.append(f'         → {fmt(vv)}')
    out.append(f'    in:')
    for u in in_adj[node][:10]:
        out.append(f'      ← {fmt(u)}')
    return '\n'.join(out)

print('=== CYCLE-IN-ALLELE cases (5 examples) ===')
for s1, s2 in [(240219, 463286), (240219, 581638), (259241, 463286), (302506, 578747)]:
    print(f'\n--- (S1={s1}, S2={s2}) ---')
    print(f'S1: {fmt(s1)}, out_deg={len(out_adj[s1])}, in_deg={len(in_adj[s1])}')
    print(f'S2: {fmt(s2)}, out_deg={len(out_adj[s2])}, in_deg={len(in_adj[s2])}')
    # BFS from s1 through color != 3 to s2, dump path nodes
    parent = {s1: None}
    bfs = deque([s1])
    found = False
    while bfs and not found:
        u = bfs.popleft()
        for v in out_adj[u]:
            if v in parent: continue
            if v == s2 or color[v] != 3 or u == s1:
                parent[v] = u
                if v == s2:
                    found = True; break
                bfs.append(v)
    if found:
        path = []
        cur = s2
        while cur is not None:
            path.append(cur)
            cur = parent.get(cur)
        path.reverse()
        print(f'  path ({len(path)} nodes): ' + ' → '.join(fmt(n) for n in path))
    # Also dump 1-hop nbhd of S1 and S2
    print(f'  S1 nbhd: out_colors={Counter(color[v] for v in out_adj[s1])}, in_colors={Counter(color[u] for u in in_adj[s1])}')
    print(f'  S2 nbhd: out_colors={Counter(color[v] for v in out_adj[s2])}, in_colors={Counter(color[u] for u in in_adj[s2])}')

print('\n\n=== N_only_no_T cases (5 examples) ===')
for s1, s2 in [(24, 321119), (319, 503312), (319, 374452), (344, 351344), (344, 338189)]:
    print(f'\n--- (S1={s1}, S2={s2}) ---')
    print(f'S1: {fmt(s1)}, out_deg={len(out_adj[s1])}, in_deg={len(in_adj[s1])}')
    print(f'S2: {fmt(s2)}, out_deg={len(out_adj[s2])}, in_deg={len(in_adj[s2])}')
    print(f'  S1 out: {[fmt(v) for v in out_adj[s1]]}')
    print(f'  S2 in:  {[fmt(u) for u in in_adj[s2]]}')
    # Does S1 have a direct color=3 edge to anywhere? Where?
    direct_targets = [v for v in out_adj[s1] if color[v] == 3]
    print(f'  S1 direct color-3 targets: {[fmt(v) for v in direct_targets[:5]]}')

print('\n\n=== T_only_no_N cases (5 examples) ===')
for s1, s2 in [(10, 95523), (40, 96577), (40, 182122), (147, 688444), (159, 159789)]:
    print(f'\n--- (S1={s1}, S2={s2}) ---')
    print(f'S1: {fmt(s1)}, out_deg={len(out_adj[s1])}, in_deg={len(in_adj[s1])}')
    print(f'S2: {fmt(s2)}, out_deg={len(out_adj[s2])}, in_deg={len(in_adj[s2])}')
    print(f'  S1 out: {[fmt(v) for v in out_adj[s1]]}')
    print(f'  S2 in:  {[fmt(u) for u in in_adj[s2]]}')
    direct_targets = [v for v in out_adj[s1] if color[v] == 3]
    print(f'  S1 direct color-3 targets: {[fmt(v) for v in direct_targets[:5]]}')

print('\n\n=== MULTI_PATH cases (5 examples) ===')
for s1, s2 in [(2353, 217479), (4446, 386954), (16064, 878812), (16923, 678894), (18655, 615391)]:
    print(f'\n--- (S1={s1}, S2={s2}) ---')
    print(f'S1: {fmt(s1)}, out_deg={len(out_adj[s1])}, in_deg={len(in_adj[s1])}')
    print(f'S2: {fmt(s2)}, out_deg={len(out_adj[s2])}, in_deg={len(in_adj[s2])}')
    print(f'  S1 out: {[fmt(v) for v in out_adj[s1]]}')
    print(f'  S2 in:  {[fmt(u) for u in in_adj[s2]]}')
