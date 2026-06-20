#!/usr/bin/env python3
"""Node-class composition + read-error depth distribution for a built graph.bin + its annot.bin.
Reads only the small header arrays (color/ncov/tcov/truthsrc) via memmap, not the payload.
Usage: depth_analysis.py <graph.bin> <annot.bin> [label]"""
import sys, numpy as np
G, A = sys.argv[1], sys.argv[2]
label = sys.argv[3] if len(sys.argv) > 3 else G
K = int(np.fromfile(G, dtype=np.uint64, count=1)[0])
color = np.memmap(G, dtype=np.uint8,  mode='r', offset=8,      shape=(K,))
ncov  = np.memmap(G, dtype=np.uint32, mode='r', offset=8+K,    shape=(K,))
tcov  = np.memmap(G, dtype=np.uint32, mode='r', offset=8+5*K,  shape=(K,))
Ka = int(np.fromfile(A, dtype=np.uint64, count=1, offset=4)[0])
assert Ka == K, (K, Ka)
ts = np.asarray(np.memmap(A, dtype=np.uint8, mode='r', offset=12, shape=(K,)))
REF, OUT = 1 | 2, 4 | 8
err     = ts == 0
both    = ((ts & REF) > 0) & ((ts & OUT) > 0)
refonly = ((ts & REF) > 0) & ((ts & OUT) == 0)
somatic = ((ts & OUT) > 0) & ((ts & REF) == 0)
print(f"==== {label} ====")
print(f"K = {K:,} nodes")
for nm, m in [('both(backbone)', both), ('ref_only', refonly), ('somatic', somatic), ('read_error', err)]:
    print(f"  {nm:16s} {int(m.sum()):>12,}  {100*m.mean():6.3f}%")
ec = np.asarray(color[err]); en = np.asarray(ncov[err]).astype(np.int64); et = np.asarray(tcov[err]).astype(np.int64)
depth = np.where(ec == 1, en, np.where(ec == 2, et, np.maximum(en, et)))
ne = int(err.sum())
print(f"\nread-error depth (surviving-colour count), n={ne:,}:")
print(f"  median {int(np.median(depth))}  mean {depth.mean():.1f}  p90 {int(np.percentile(depth,90))}  p99 {int(np.percentile(depth,99))}  max {int(depth.max())}")
print(f"  ==10(floor) {100*(depth==10).mean():.1f}%   <=12 {100*(depth<=12).mean():.1f}%   >=50 {100*(depth>=50).mean():.1f}%")
for c, nm in [(1, 'normal-only'), (2, 'tumor-only'), (3, 'shared')]:
    mc = ec == c
    if mc.any(): print(f"    err c{c} {nm:11s} {int(mc.sum()):>10,} ({100*mc.mean():4.1f}%) median depth {int(np.median(depth[mc]))}")
# somatic vs colour-2 error depth (the blind-precision lever)
sm = somatic
if sm.any():
    sd = np.asarray(np.maximum(np.asarray(ncov[sm]), np.asarray(tcov[sm]))).astype(np.int64)
    st = np.asarray(tcov[sm]).astype(np.int64)
    print(f"\nsomatic node tumor-depth: median {int(np.median(st))} mean {st.mean():.1f}  (vs colour-2 error median {int(np.median(et[ec==2])) if (ec==2).any() else 0})")
