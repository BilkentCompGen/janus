#!/usr/bin/env python3
# modea_audit.py — classify the bubble shapes that actually exist in modea.bin so we
# know which of the "theoretical" Mode-B exceptional cases (D-only, I-only, multi-path,
# super-bubble, cycle-in-allele) are real and how often.

import struct
import sys
from collections import defaultdict, Counter

MODEA = sys.argv[1] if len(sys.argv) > 1 else 'outputs/svp/v2_test/modea.bin'

with open(MODEA, 'rb') as f:
    data = f.read()
print(f'mmap modea.bin {len(data)/1e9:.2f} GB', file=sys.stderr)

off = 0
S = struct.unpack_from('<Q', data, off)[0]; off += 8
color = bytearray(S)
for i in range(S):
    cl, bl, rl, c = struct.unpack_from('<IIIB', data, off)
    color[i] = c
    off += 16
print(f'S super-nodes = {S:,}  color counts: {dict(Counter(color))}', file=sys.stderr)

tc = struct.unpack_from('<Q', data, off)[0]; off += 8 + tc*4
to = struct.unpack_from('<Q', data, off)[0]; off += 8 + to*2
tb = struct.unpack_from('<Q', data, off)[0]; off += 8 + tb
tr = struct.unpack_from('<Q', data, off)[0]; off += 8 + tr*12
E = struct.unpack_from('<Q', data, off)[0]; off += 8
print(f'super-edges E = {E:,}', file=sys.stderr)

out_adj = [[] for _ in range(S)]
for k in range(E):
    u, v = struct.unpack_from('<II', data, off + k*8)
    out_adj[u].append(v)
print('built adjacency', file=sys.stderr)

# For each color=3 super-node S1, forward-BFS through color != 3 successors.
# Stop at any color=3 node S2 (record as bubble closure). Per allele-branch
# track which "branch color" it started with (color of S1's immediate non-anchor
# out-neighbor). Also record any direct S1->S2 color3->color3 edges.

# Aggregations:
NODE_BUDGET = 500   # per-S1 walk node cap (truncate exotic super-bubbles, count separately)

pair_info = defaultdict(lambda: {'c1':0,'c2':0,'mix':0,'direct':False,'inter_anchor':False,'cycle':False})
truncated_starts = 0
n_S1_with_outedge = 0
n_direct_anchor_edges = 0
direct_pairs = set()
n_intermediate_anchor_paths = 0

# track examples for each exceptional class (cap collected to 10 each)
examples = defaultdict(list)
def add_example(tag, payload):
    if len(examples[tag]) < 10:
        examples[tag].append(payload)

for s1 in range(S):
    if color[s1] != 3:
        continue
    starts_nontrivial = False
    for x in out_adj[s1]:
        if color[x] == 3:
            n_direct_anchor_edges += 1
            direct_pairs.add((s1, x))
            pair_info[(s1, x)]['direct'] = True
            continue
        starts_nontrivial = True
    if not starts_nontrivial and not out_adj[s1]:
        continue
    n_S1_with_outedge += 1 if starts_nontrivial else 0

    # BFS per allele-branch (start color)
    for x in out_adj[s1]:
        if color[x] == 3:
            continue
        branch_c = color[x]                 # 1 or 2
        visited = {x}
        # BFS — but emit closure once per s2 hit; also dedupe by branch_c
        stack = [x]
        nodes = 0
        closed = set()
        saw_intermediate_anchor = False
        cycle_hit = False
        while stack and nodes < NODE_BUDGET:
            n = stack.pop()
            nodes += 1
            for nx in out_adj[n]:
                if nx == s1:
                    cycle_hit = True
                    continue
                if color[nx] == 3:
                    closed.add(nx)
                else:
                    if nx in visited:
                        cycle_hit = True
                        continue
                    visited.add(nx)
                    stack.append(nx)
        if nodes >= NODE_BUDGET:
            truncated_starts += 1
            add_example('truncated', (s1, x, branch_c, nodes))
        for s2 in closed:
            key = (s1, s2)
            if branch_c == 1: pair_info[key]['c1'] += 1
            else:             pair_info[key]['c2'] += 1
            if cycle_hit:     pair_info[key]['cycle'] = True

# Now classify each (s1, s2) closure pair
classes = Counter()
for (s1, s2), v in pair_info.items():
    has1 = v['c1'] > 0
    has2 = v['c2'] > 0
    direct = v['direct']
    cycle = v['cycle']
    if has1 and has2:
        cls = 'std_both'
        if v['c1'] > 1 or v['c2'] > 1:
            cls = 'multi_path_both'
        if direct:
            cls = 'std_both_plus_direct'
    elif has1 and not has2:
        cls = 'N_only_plus_direct' if direct else 'N_only_no_T'
    elif has2 and not has1:
        cls = 'T_only_plus_direct' if direct else 'T_only_no_N'
    else:
        # no allele path; just a direct anchor edge
        cls = 'direct_only'
    if cycle and cls != 'direct_only':
        cls = cls + '_cycle'
    classes[cls] += 1
    add_example(cls, (s1, s2, v))

print('\n== Bucket counts (per unique (S1, S2) anchor pair) ==')
total = sum(classes.values())
for cls, n in sorted(classes.items(), key=lambda kv: -kv[1]):
    print(f'  {cls:30s} = {n:>10,}  ({100.0*n/total:6.2f}%)')

print(f'\nTotal closures = {total:,}')
print(f'Truncated walks (>{NODE_BUDGET} nodes) = {truncated_starts:,}')
print(f'Direct color3->color3 edges = {n_direct_anchor_edges:,}  (distinct pairs: {len(direct_pairs):,})')

print('\n== examples (first 5 of each class) ==')
for cls, exs in examples.items():
    print(f'  [{cls}]')
    for e in exs[:5]:
        print(f'     {e}')
