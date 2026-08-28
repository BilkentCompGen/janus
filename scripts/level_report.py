#!/usr/bin/env python3
"""Full per-level analysis (identical metrics for L5 and L4 -> apples-to-apples).
Usage: level_report.py <graph.bin> <trace_dir> <label>
Reads <trace_dir>/{variants.tsv,diag.tsv,bubbles.tsv,error_clusters.tsv,annot.bin}."""
import sys, os, json, collections
import numpy as np

GRAPH, TDIR, LABEL = sys.argv[1], sys.argv[2], sys.argv[3]
def rows(p):
    if not os.path.exists(p): return []
    f=open(p); head=f.readline().rstrip("\n").split("\t")
    return [dict(zip(head,l.rstrip("\n").split("\t"))) for l in f]

print(f"\n################  {LABEL}  ################")

# ---- node classes + error depth (from graph.bin + annot.bin) ----
K=int(np.fromfile(GRAPH,dtype=np.uint64,count=1)[0])
color=np.memmap(GRAPH,dtype=np.uint8,mode='r',offset=8,shape=(K,))
ncov =np.memmap(GRAPH,dtype=np.uint32,mode='r',offset=8+K,shape=(K,))
tcov =np.memmap(GRAPH,dtype=np.uint32,mode='r',offset=8+5*K,shape=(K,))
ts=np.asarray(np.memmap(TDIR+"/annot.bin",dtype=np.uint8,mode='r',offset=12,shape=(K,)))
REF,OUT=3,12
err=ts==0; both=((ts&REF)>0)&((ts&OUT)>0); refonly=((ts&REF)>0)&((ts&OUT)==0); somatic=((ts&OUT)>0)&((ts&REF)==0)
E=int(np.fromfile(GRAPH,dtype=np.uint64,count=1,offset=8+9*K)[0])
print(f"\n[graph]  K={K:,} nodes  E={E:,} edges")
print(f"[node classes]")
for nm,m in [('both(backbone)',both),('ref_only',refonly),('somatic',somatic),('read_error',err)]:
    print(f"   {nm:16s} {int(m.sum()):>12,}  {100*m.mean():6.3f}%")

# ---- marking-1: raw variant recall + by type + colour agreement ----
V=rows(TDIR+"/variants.tsv")
byt=collections.defaultdict(lambda:[0,0]); agree=[0,0]; tot=[0,0]
for v in V:
    t=v['type']; byt[t][1]+=1; tot[1]+=1
    if v['found']=='1':
        byt[t][0]+=1; tot[0]+=1; agree[1]+=1; agree[0]+=(v['best_color']=='2')
print(f"\n[marking-1 RAW recall (located by ref-neighbour attribution)]  {tot[0]:,}/{tot[1]:,} = {100*tot[0]/tot[1]:.2f}%")
for t in ['SNP','DEL','INS','DUP','INV']:
    if t in byt: a,b=byt[t]; print(f"   {t:4s} {a:>6,}/{b:<6,}  {100*a/b:5.1f}%")
print(f"   somatic<->colour-2 agreement: {agree[0]:,}/{agree[1]:,} = {100*agree[0]/max(1,agree[1]):.3f}%")

# ---- diag: TRUE recall decomposition ----
D={d['svid']:d['bucket'] for d in rows(TDIR+"/diag.tsv") if d['found']=='0'}
located=tot[0]
gap=sum(1 for v in V if v['found']=='0' and D.get(v['svid'])=='sig_present_attribution_gap')
absent_buckets=collections.Counter(D.get(v['svid']) for v in V if v['found']=='0' and D.get(v['svid'])!='sig_present_attribution_gap')
absent=sum(absent_buckets.values())
TT=tot[1]
print(f"\n[TRUE recall decomposition (diag)]")
print(f"   located (found=1)            {located:>7,}  {100*located/TT:5.2f}%")
print(f"   sig_present (attribution gap){gap:>7,}  {100*gap/TT:5.2f}%   <- signature IS in graph")
print(f"   truly_absent                 {absent:>7,}  {100*absent/TT:5.2f}%")
print(f"   ===> TRUE RECALL (located+sig_present)/total = {100*(located+gap)/TT:.2f}%   (raw was {100*located/TT:.2f}%)")
print(f"   truly_absent buckets: "+" ".join(f"{k}={v}" for k,v in absent_buckets.most_common() if k))
# true recall by type
print(f"   by type (located | gap | absent | true_recall):")
for t in ['SNP','DEL','INS','DUP','INV']:
    vt=[v for v in V if v['type']==t]
    loc=sum(1 for v in vt if v['found']=='1'); g=sum(1 for v in vt if v['found']=='0' and D.get(v['svid'])=='sig_present_attribution_gap')
    ab=len(vt)-loc-g
    if vt: print(f"     {t:4s} {loc:>6,} | {g:>5,} | {ab:>4,} | {100*(loc+g)/len(vt):5.1f}%")

# ---- marking-2: error clusters + depth ----
EC=rows(TDIR+"/error_clusters.tsv")
ec=np.asarray(color[err]); en=np.asarray(ncov[err]).astype(np.int64); et=np.asarray(tcov[err]).astype(np.int64)
depth=np.where(ec==1,en,np.where(ec==2,et,np.maximum(en,et))); ne=int(err.sum())
print(f"\n[marking-2 read errors]  {ne:,} nodes  in {len(EC):,} clusters")
print(f"   depth: median {int(np.median(depth))}  mean {depth.mean():.1f}  p99 {int(np.percentile(depth,99))}  max {int(depth.max())}")
print(f"   ==10 {100*(depth==10).mean():.1f}%   <=12 {100*(depth<=12).mean():.1f}%   >=50 {100*(depth>=50).mean():.1f}%")
for c,nm in [(1,'normal-only'),(2,'tumor-only'),(3,'shared')]:
    mc=ec==c
    if mc.any(): print(f"     c{c} {nm:11s} {int(mc.sum()):>9,} ({100*mc.mean():4.1f}%) median {int(np.median(depth[mc]))}")

# ---- marking-3: bubbles ----
B=rows(TDIR+"/bubbles.tsv")
if B:
    clean=sum(1 for b in B if b['clean']=='1'); cyc=sum(1 for b in B if b['kind']=='cycle')
    matched=sum(1 for b in B if b['matched_svid']!='NA'); bk=collections.Counter(b['bucket'] for b in B)
    print(f"\n[marking-3 bubbles/cycles]  {len(B):,} total  clean {clean:,} ({100*clean/len(B):.1f}%)  contaminated {len(B)-clean:,}  cycles {cyc}  matched-a-variant {matched:,} ({100*matched/len(B):.1f}%)")
    print(f"   buckets: "+" ".join(f"{k}={v:,}" for k,v in bk.items()))

# ---- somatic vs colour-2 error depth (precision lever) ----
if somatic.any():
    st=np.asarray(tcov[somatic]).astype(np.int64)
    print(f"\n[precision lever] somatic tumor-depth median {int(np.median(st))} mean {st.mean():.1f}  vs colour-2 error median {int(np.median(et[ec==2])) if (ec==2).any() else 0}")
