#!/usr/bin/env python3
# Inspect direct_only pairs — why does Mode A leave two adjacent color=3 super-nodes
# unmerged? Check the branching structure at their boundary.

import struct
import sys
from collections import defaultdict, Counter

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
    color[i] = c; chain_len[i] = cl; block_len[i] = bl
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

# Distribution of out_deg / in_deg over color=3 super-nodes
out_deg_dist = Counter()
in_deg_dist = Counter()
both1 = 0
for u in range(S):
    if color[u] != 3: continue
    od = len(out_adj[u]); idg = len(in_adj[u])
    out_deg_dist[min(od, 10)] += 1
    in_deg_dist[min(idg, 10)] += 1
    if od == 1 and idg == 1: both1 += 1

print('color=3 super-node out_deg distribution (cap 10):')
for k in sorted(out_deg_dist): print(f'  od={k}: {out_deg_dist[k]:>10,}')
print('color=3 super-node in_deg distribution (cap 10):')
for k in sorted(in_deg_dist): print(f'  id={k}: {in_deg_dist[k]:>10,}')
print(f'\ncolor=3 super-nodes with out_deg==1 AND in_deg==1: {both1:,}')

# For each color=3 → color=3 direct edge u→v, classify why it didn't merge:
#   (a) out_deg[u] > 1 (u branches outwards)
#   (b) in_deg[v] > 1 (v has multiple predecessors)
#   (c) both
#   (d) neither — should be 0; if not, Mode A bug
no_branch = 0; only_u_branches = 0; only_v_branches = 0; both_branch = 0
for u in range(S):
    if color[u] != 3: continue
    od = len(out_adj[u])
    for v in out_adj[u]:
        if color[v] != 3: continue
        idv = len(in_adj[v])
        if od > 1 and idv > 1: both_branch += 1
        elif od > 1: only_u_branches += 1
        elif idv > 1: only_v_branches += 1
        else: no_branch += 1
print(f'\nReason for color3→color3 direct edges remaining unmerged:')
print(f'  u branches AND v joins (both >1):  {both_branch:,}')
print(f'  u branches only (od>1, idv=1):     {only_u_branches:,}')
print(f'  v joins only   (od=1, idv>1):      {only_v_branches:,}')
print(f'  neither branches (Mode A bug?):    {no_branch:,}')

# Dump 5 direct_only examples with neighborhood for diagram
print('\n=== Example direct_only triples (S1 → S2 with neighborhoods) ===')
seen = 0
for u in range(S):
    if seen >= 5: break
    if color[u] != 3: continue
    if len(out_adj[u]) < 2: continue
    for v in out_adj[u]:
        if color[v] != 3: continue
        if len(in_adj[v]) < 2: continue
        print(f'\n--- direct: sn{u}(chain={chain_len[u]},bp={block_len[u]}) → sn{v}(chain={chain_len[v]},bp={block_len[v]}) ---')
        print(f'  sn{u} out_deg={len(out_adj[u])} → {[f"sn{x}(c{color[x]},bp{block_len[x]})" for x in out_adj[u][:8]]}')
        print(f'  sn{v} in_deg={len(in_adj[v])}  ← {[f"sn{x}(c{color[x]},bp{block_len[x]})" for x in in_adj[v][:8]]}')
        seen += 1
        break
