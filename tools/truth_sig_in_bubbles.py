#!/usr/bin/env python3
# For each truth somatic sig that landed in a color=2 super-node, identify the (S1, S2)
# anchor pair it belongs to and the bubble shape. Reports whether real DELs/INSs are in
# std_both bubbles (canonical) or T_only_no_N orphans (no parallel N path).
#
# Pipeline:
#   1. Load truth somatic sigs S from a precomputed file (or from add_sigs via the C++ tool).
#   2. Load modea.bin, counts.bin → node_to_super, super-edges.
#   3. For each color=2 super-node sn containing a truth sig:
#        - Look at sn's in-neighbors: which anchors (color=3) can reach sn?
#        - Look at sn's out-neighbors: which anchors does sn reach?
#        - For each (S1, S2) bracketing sn: check if a parallel color=1 path exists.
#   4. Tabulate.

import struct
import sys
import subprocess
import os
import tempfile
from collections import defaultdict, Counter

MODEA = sys.argv[1] if len(sys.argv) > 1 else 'outputs/svp/v2_test/modea.bin'
COUNTS = sys.argv[2] if len(sys.argv) > 2 else 'outputs/svp/v2_test/counts.bin'
NFA = sys.argv[3] if len(sys.argv) > 3 else 'simulation/ref.fasta'
TFA = sys.argv[4] if len(sys.argv) > 4 else 'simulation/truth/wgs/out.fa'

# Build a one-shot extractor C++ binary or just rerun truth_sig_localizer with stderr
# capture. To avoid recomputing |S|, dump it via a side-channel binary. Cheaper:
# write a small companion tool that prints the truth set.
# For this session: rerun truth extraction in-process via a C++ helper that dumps
# (sig, node_id) pairs. But we already ran truth_sig_localizer; instead, let's spawn
# a small companion that just dumps truth somatic sig hex list to stdout.

# Simplest: regenerate truth sigs by calling a tool that dumps them
DUMP = '/tmp/truth_sigs.bin'
if not os.path.exists(DUMP) or os.path.getsize(DUMP) == 0:
    print('Need truth-sig dump — running C++ helper truth_sig_dump...', file=sys.stderr)
    sys.exit('Run tools/truth_sig_dump first to dump |S| to ' + DUMP)

# Load truth sigs (one u64 per record)
with open(DUMP, 'rb') as f:
    raw = f.read()
truth = set()
for i in range(0, len(raw), 8):
    truth.add(struct.unpack_from('<Q', raw, i)[0])
print(f'loaded |S| = {len(truth):,} from {DUMP}', file=sys.stderr)

# Load counts.bin (header + SigCount[K] only; skip bytes section)
with open(COUNTS, 'rb') as f:
    K = struct.unpack('<Q', f.read(8))[0]
    sigs_raw = f.read(K * 16)
sig2node = {}
sig_arr = []
for i in range(K):
    s, n, t = struct.unpack_from('<QII', sigs_raw, i*16)
    sig2node[s] = i
    sig_arr.append((n, t))
print(f'loaded counts.bin K = {K:,}', file=sys.stderr)

# Load modea.bin
with open(MODEA, 'rb') as f:
    data = f.read()
off = 0
S = struct.unpack_from('<Q', data, off)[0]; off += 8
sn_color = bytearray(S)
sn_chain_len = [0]*S
for i in range(S):
    cl, bl, rl, c = struct.unpack_from('<IIIB', data, off)
    sn_color[i] = c
    sn_chain_len[i] = cl
    off += 16
tc = struct.unpack_from('<Q', data, off)[0]; off += 8
chain_ids_raw = data[off:off+tc*4]; off += tc*4
to = struct.unpack_from('<Q', data, off)[0]; off += 8 + to*2
tb = struct.unpack_from('<Q', data, off)[0]; off += 8 + tb
tr = struct.unpack_from('<Q', data, off)[0]; off += 8 + tr*12
E = struct.unpack_from('<Q', data, off)[0]; off += 8

# Build node→super mapping
node2super = [0xFFFFFFFF] * K
p = 0
for s in range(S):
    m = sn_chain_len[s]
    for j in range(m):
        nid = struct.unpack_from('<I', chain_ids_raw, p*4)[0]
        node2super[nid] = s
        p += 1

# Build super-edges
out_adj = [[] for _ in range(S)]
in_adj  = [[] for _ in range(S)]
for k in range(E):
    u, v = struct.unpack_from('<II', data, off + k*8)
    out_adj[u].append(v); in_adj[v].append(u)
print(f'loaded super-graph: S={S:,} E={E:,}', file=sys.stderr)

# Find truth-sig-containing super-nodes
truth_super = set()
for s in truth:
    nid = sig2node.get(s)
    if nid is None: continue
    sid = node2super[nid]
    if sn_color[sid] == 2: truth_super.add(sid)
print(f'color=2 super-nodes containing >=1 truth sig: {len(truth_super):,}', file=sys.stderr)

# For each truth-bearing color=2 super-node, find its (S1, S2) bracket and classify
# the bubble.
NODE_BUDGET = 200

def find_brackets(sn):
    # Walk backward from sn through color != 3 to find color=3 left anchors S1.
    lefts = set()
    stack = [sn]; visited = {sn}; budget = NODE_BUDGET
    while stack and budget > 0:
        u = stack.pop(); budget -= 1
        for v in in_adj[u]:
            if sn_color[v] == 3: lefts.add(v)
            elif v not in visited: visited.add(v); stack.append(v)
    # Walk forward from sn through color != 3 to find color=3 right anchors S2.
    rights = set()
    stack = [sn]; visited = {sn}; budget = NODE_BUDGET
    while stack and budget > 0:
        u = stack.pop(); budget -= 1
        for v in out_adj[u]:
            if sn_color[v] == 3: rights.add(v)
            elif v not in visited: visited.add(v); stack.append(v)
    return lefts, rights

def has_color1_path(s1, s2, max_hops=15):
    # From s1, walk through color=1 (and other non-3) to see if any path closes at s2.
    if s1 == s2: return False
    stack = [(s1, 0)]; visited = {s1}; budget = NODE_BUDGET
    while stack and budget > 0:
        u, d = stack.pop(); budget -= 1
        if d >= max_hops: continue
        for v in out_adj[u]:
            if v == s2: return True
            if sn_color[v] == 3 or v in visited: continue
            if sn_color[v] != 1: continue  # only color=1 path
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
            if v in visited: continue
            if sn_color[v] != 3: continue
            visited.add(v); stack.append((v, d+1))
    return False

verdict = Counter()
sample = list(truth_super)
# audit a sample for speed
import random
random.seed(42)
if len(sample) > 5000: sample = random.sample(sample, 5000)
for sn in sample:
    lefts, rights = find_brackets(sn)
    # For each (L, R) pair: classify
    if not lefts or not rights:
        verdict['no_brackets'] += 1; continue
    found_std = False; found_orphan = False; found_detour = False
    for L in lefts:
        for R in rights:
            n_present = has_color1_path(L, R)
            if n_present: found_std = True
            elif has_color3_path(L, R): found_detour = True
            else: found_orphan = True
    if found_std: verdict['std_both'] += 1
    elif found_detour: verdict['T_only_with_c3_detour'] += 1
    else: verdict['T_only_no_N_no_detour'] += 1

print(f'\n== Truth-bearing color=2 super-nodes — bubble shape (sample of {len(sample):,}) ==')
total = sum(verdict.values())
for k, v in sorted(verdict.items(), key=lambda kv: -kv[1]):
    print(f'  {k:30s} = {v:>7,}  ({100*v/total:5.2f}%)')
