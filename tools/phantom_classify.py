#!/usr/bin/env python3
# Classify every Mode-A bubble shape (std_both / N_only_no_T / T_only_no_N) by what
# its constituent non-anchor (color != 3) super-nodes are made of, using per-sig
# genomic provenance.
#
# Provenance buckets per sig (from sig_provenance.bin):
#   GENUINE_BOTH       — n>0 and t>0 in assemblies (would be color=3 in a clean world)
#   GENUINE_N_ONLY     — only normal genome
#   GENUINE_T_ONLY     — only tumor genome    (this is exactly |S|)
#   TRUE_PHANTOM       — 0 in all 4 assembly buckets
#   REPETITIVE         — any single (sample,strand) count > 1
#
# A super-node may contain multiple sigs of different buckets; we record a vector
# count per super-node and assign a primary class.

import struct
import sys
from collections import defaultdict, Counter

MODEA  = sys.argv[1] if len(sys.argv) > 1 else 'outputs/svp/v2_test/modea.bin'
COUNTS = sys.argv[2] if len(sys.argv) > 2 else 'outputs/svp/v2_test/counts.bin'
PROV   = sys.argv[3] if len(sys.argv) > 3 else '/tmp/sig_provenance.bin'

# Load provenance (parallel to counts.bin order).
with open(PROV, 'rb') as f:
    K = struct.unpack('<Q', f.read(8))[0]
    raw = f.read(K * 16)
print(f'[ph] loaded provenance: K = {K:,}', file=sys.stderr)
# Each record: 4 × u32. Stored in counts.bin sig order, so index == node_id.

def prov_at(i):
    return struct.unpack_from('<IIII', raw, i*16)  # (n_fwd, n_rc, t_fwd, t_rc)

# Classify each sig's provenance bucket.
sig_class = bytearray(K)   # 0=GENUINE_BOTH, 1=N_only, 2=T_only, 3=phantom, 4=repetitive_extra
CLASS = {0:'GENUINE_BOTH', 1:'GENUINE_N_ONLY', 2:'GENUINE_T_ONLY', 3:'TRUE_PHANTOM'}
is_repetitive = bytearray(K)
for i in range(K):
    nf, nr, tf, tr = prov_at(i)
    n = (nf + nr) > 0
    t = (tf + tr) > 0
    if   n and t: sig_class[i] = 0
    elif n:       sig_class[i] = 1
    elif t:       sig_class[i] = 2
    else:         sig_class[i] = 3
    if max(nf, nr, tf, tr) > 1: is_repetitive[i] = 1
print('[ph] per-sig provenance class:', dict(Counter(sig_class)), file=sys.stderr)
print(f'[ph] repetitive sigs: {sum(is_repetitive):,}', file=sys.stderr)

# Load modea.bin to get node_id -> super_id mapping and super-edges.
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

# Build node->super_id and super->[node_ids]
node2super = [0xFFFFFFFF] * K
super_nodes = [[] for _ in range(S)]
p = 0
for s in range(S):
    m = sn_chain_len[s]
    for j in range(m):
        nid = struct.unpack_from('<I', chain_ids_raw, p*4)[0]
        node2super[nid] = s
        super_nodes[s].append(nid)
        p += 1
print(f'[ph] super-graph loaded: S = {S:,}, E = {E:,}', file=sys.stderr)

out_adj = [[] for _ in range(S)]
in_adj  = [[] for _ in range(S)]
for k in range(E):
    u, v = struct.unpack_from('<II', data, off + k*8)
    out_adj[u].append(v); in_adj[v].append(u)

# Per-super-node sig-class summary.
def super_class(sid):
    cnt = Counter()
    rep = 0
    for nid in super_nodes[sid]:
        cnt[sig_class[nid]] += 1
        if is_repetitive[nid]: rep += 1
    return cnt, rep

# Re-derive the bubble pair classification (same logic as modea_audit.py).
NODE_BUDGET = 500
pair_info = defaultdict(lambda: {'c1':0,'c2':0,'direct':False})
for s1 in range(S):
    if sn_color[s1] != 3: continue
    for x in out_adj[s1]:
        if sn_color[x] == 3:
            pair_info[(s1, x)]['direct'] = True
            continue
        branch_c = sn_color[x]
        visited = {x}; stack = [x]; budget = NODE_BUDGET
        closed = set()
        while stack and budget > 0:
            n = stack.pop(); budget -= 1
            for nx in out_adj[n]:
                if sn_color[nx] == 3: closed.add(nx)
                elif nx not in visited: visited.add(nx); stack.append(nx)
        for s2 in closed:
            key = (s1, s2)
            if branch_c == 1: pair_info[key]['c1'] += 1
            else:             pair_info[key]['c2'] += 1

# Bucketize pairs.
pair_class = {}
for (s1, s2), v in pair_info.items():
    has1 = v['c1'] > 0; has2 = v['c2'] > 0; direct = v['direct']
    if has1 and has2:
        cls = 'std_both_plus_direct' if direct else 'std_both'
    elif has1 and not has2:
        cls = 'N_only_plus_direct' if direct else 'N_only_no_T'
    elif has2 and not has1:
        cls = 'T_only_plus_direct' if direct else 'T_only_no_N'
    else:
        cls = 'direct_only'
    pair_class[(s1, s2)] = cls

print(f'[ph] pairs classified: {len(pair_class):,}', file=sys.stderr)

# For each non-direct_only pair, walk between (S1, S2) again and collect the color≠3
# super-nodes along all paths. Record their provenance class.

def interior_supernodes(s1, s2):
    """Return set of color != 3 super-nodes on any path from s1's color != 3 successors
    that ultimately closes at s2."""
    # Forward reachability set from s1 (through color != 3) AND backward reachability
    # set from s2 (through color != 3). Intersection = "interior".
    fwd = set()
    stack = [x for x in out_adj[s1] if sn_color[x] != 3]
    fwd.update(stack)
    budget = NODE_BUDGET
    while stack and budget > 0:
        u = stack.pop(); budget -= 1
        for v in out_adj[u]:
            if sn_color[v] != 3 and v not in fwd:
                fwd.add(v); stack.append(v)
    bwd = set()
    stack = [x for x in in_adj[s2] if sn_color[x] != 3]
    bwd.update(stack)
    budget = NODE_BUDGET
    while stack and budget > 0:
        u = stack.pop(); budget -= 1
        for v in in_adj[u]:
            if sn_color[v] != 3 and v not in bwd:
                bwd.add(v); stack.append(v)
    return fwd & bwd

# For tractability we sample.
import random
random.seed(42)
SAMPLE_PER_CLASS = 3000

samples_by_class = defaultdict(list)
for k, v in pair_class.items():
    samples_by_class[v].append(k)
for c in samples_by_class:
    if len(samples_by_class[c]) > SAMPLE_PER_CLASS:
        samples_by_class[c] = random.sample(samples_by_class[c], SAMPLE_PER_CLASS)

results = {}
for cls, pairs in samples_by_class.items():
    # Per interior super-node, record what its sigs are made of.
    sn_class_counts = Counter()           # primary class of each super-node
    sn_with_phantom = 0
    sn_with_repetitive = 0
    sn_pure_phantom = 0
    sn_pure_repetitive = 0
    total_super = 0
    sigs_total = 0
    sigs_phantom = 0
    sigs_repetitive = 0
    pure_phantom_pairs = 0
    pure_repetitive_pairs = 0
    sample_with_truth = 0
    pure_genuine_pairs = 0
    for (s1, s2) in pairs:
        interior = interior_supernodes(s1, s2)
        if not interior: continue
        pair_has_phantom_all = True
        pair_has_repetitive_all = True
        pair_has_truth_T = False
        pair_has_only_genuine_TorN = True
        for sid in interior:
            cnt, rep = super_class(sid)
            total_super += 1
            sigs_total += sum(cnt.values())
            sigs_phantom += cnt.get(3, 0)
            sigs_repetitive += rep
            # primary class of super-node = majority
            primary = max(cnt, key=cnt.get)
            sn_class_counts[primary] += 1
            if cnt.get(3, 0) > 0: sn_with_phantom += 1
            if rep > 0:           sn_with_repetitive += 1
            if cnt.get(3, 0) == sum(cnt.values()): sn_pure_phantom += 1
            else:                 pair_has_phantom_all = False
            if rep == sum(cnt.values()): sn_pure_repetitive += 1
            else:                 pair_has_repetitive_all = False
            if cnt.get(2, 0) > 0: pair_has_truth_T = True
            if cnt.get(3, 0) > 0: pair_has_only_genuine_TorN = False
        if pair_has_truth_T: sample_with_truth += 1
        if pair_has_phantom_all:    pure_phantom_pairs += 1
        if pair_has_repetitive_all: pure_repetitive_pairs += 1
        if pair_has_only_genuine_TorN: pure_genuine_pairs += 1
    results[cls] = {
        'n_pairs':            len(pairs),
        'sample_with_truth_T_sig': sample_with_truth,
        'total_interior_super_nodes': total_super,
        'super_node_class_dist': dict(sn_class_counts),
        'super_nodes_with_phantom_sig': sn_with_phantom,
        'super_nodes_with_repetitive_sig': sn_with_repetitive,
        'super_nodes_PURE_phantom':    sn_pure_phantom,
        'super_nodes_PURE_repetitive': sn_pure_repetitive,
        'pairs_PURE_phantom':    pure_phantom_pairs,
        'pairs_PURE_repetitive': pure_repetitive_pairs,
        'pairs_NO_phantom_sigs': pure_genuine_pairs,
        'total_sigs':         sigs_total,
        'phantom_sigs':       sigs_phantom,
        'repetitive_sigs':    sigs_repetitive,
    }

print('\n== Phantom-bubble attribution ==')
for cls in ['std_both', 'std_both_plus_direct', 'N_only_no_T', 'T_only_no_N',
            'N_only_plus_direct', 'T_only_plus_direct', 'direct_only']:
    if cls not in results: continue
    r = results[cls]
    print(f'\n[{cls}]   sample n_pairs = {r["n_pairs"]:,}')
    npairs = max(1, r['n_pairs'])
    print(f'  pairs with >=1 truth-T-sig super-node = {r["sample_with_truth_T_sig"]:>6,}  ({100*r["sample_with_truth_T_sig"]/npairs:5.1f}%)')
    print(f'  pairs with NO TRUE_PHANTOM sigs       = {r["pairs_NO_phantom_sigs"]:>6,}  ({100*r["pairs_NO_phantom_sigs"]/npairs:5.1f}%)')
    print(f'  pairs where ALL interior sigs are TRUE_PHANTOMS    = {r["pairs_PURE_phantom"]:>6,}  ({100*r["pairs_PURE_phantom"]/npairs:5.1f}%)')
    print(f'  pairs where ALL interior sigs are REPETITIVE       = {r["pairs_PURE_repetitive"]:>6,}  ({100*r["pairs_PURE_repetitive"]/npairs:5.1f}%)')
    print(f'  total interior super-nodes (across all pairs)      = {r["total_interior_super_nodes"]:,}')
    cdist = r['super_node_class_dist']
    s = max(1, sum(cdist.values()))
    print(f'  super-node primary class dist:')
    for k, v in sorted(cdist.items()):
        print(f'      {CLASS.get(k,k):20s} = {v:>8,}  ({100*v/s:5.1f}%)')
    sn_total = max(1, r['total_interior_super_nodes'])
    print(f'  super-nodes touching >=1 TRUE_PHANTOM sig  = {r["super_nodes_with_phantom_sig"]:,}  ({100*r["super_nodes_with_phantom_sig"]/sn_total:5.1f}%)')
    print(f'  super-nodes touching >=1 REPETITIVE sig    = {r["super_nodes_with_repetitive_sig"]:,}  ({100*r["super_nodes_with_repetitive_sig"]/sn_total:5.1f}%)')
    st = max(1, r['total_sigs'])
    print(f'  total interior sigs = {r["total_sigs"]:,}; TRUE_PHANTOM = {r["phantom_sigs"]:,} ({100*r["phantom_sigs"]/st:5.2f}%); REPETITIVE = {r["repetitive_sigs"]:,} ({100*r["repetitive_sigs"]/st:5.2f}%)')
