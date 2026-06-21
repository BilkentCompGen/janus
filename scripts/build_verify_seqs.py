#!/usr/bin/env python3
"""Independent verification: hand-construct alt-allele junction/body/interior sequences for a
sample of truth variants via a SEPARATE code path (python + samtools faidx, not the C++
build_alleles), so `trace probeseq` re-derives each region's node colours independently.

Emits  <out>/verify_seqs.tsv     (label <TAB> SEQUENCE  -> feed to `trace probeseq`)
       <out>/verify_targets.json  (label -> {type,svid,svlen,role,expect})

role/expectation:
  seam (junc/alt/left/right)  -> expect a colour-2 / out-only somatic core (the novel junction)
  body  (dup_body)            -> expect colour-3 SHARED (duplicated unit already in normal)
  delbody (large del segment) -> expect colour-1 NORMAL-only (deleted, absent from tumor)
  interior (inv interior)     -> expect colour-3 SHARED (inverted body, same up to strand)

Usage: build_verify_seqs.py <ref.fasta> <out.vcf> <outdir> [CTX=600]
"""
import sys, os, subprocess, collections, json

REF, VCF, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
CTX = int(sys.argv[4]) if len(sys.argv) > 4 else 600
SAMTOOLS = os.path.expanduser("~/.local/bin/samtools")
if not os.path.exists(SAMTOOLS): SAMTOOLS = "samtools"

# chrom lengths from .fai
chrlen = {}
with open(REF + ".fai") as f:
    for ln in f:
        c = ln.split("\t"); chrlen[c[0]] = int(c[1])

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]

def info_get(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="): return kv[len(key)+1:]
    return ""

# pick a deterministic, size-stratified sample of variants
buckets = collections.defaultdict(list)   # bucket-key -> list of dicts
with open(VCF) as f:
    for ln in f:
        if ln[0] == "#": continue
        fl = ln.rstrip("\n").split("\t")
        if len(fl) < 8: continue
        chrom = fl[0]
        if chrom not in chrlen: continue
        pos0 = int(fl[1]) - 1
        info = fl[7]; svt = info_get(info, "SVTYPE")
        svlen = int(info_get(info, "SVLEN") or 0)
        endv = info_get(info, "END"); end0 = int(endv) if endv else pos0
        if svt not in ("DEL","INS","DUP","INV"): svt, svlen, end0 = "SNP", 0, pos0
        v = dict(chrom=chrom, pos0=pos0, end0=end0, svlen=svlen, svid=fl[2],
                 alt=fl[4], inseq=info_get(info, "INSSEQ"))
        # need room for flanks
        if pos0 < CTX + 20 or end0 + CTX + 20 > chrlen[chrom]: continue
        if svt == "SNP":
            if len(fl[3]) == 1 and len(fl[4]) == 1 and fl[4] in "ACGT": buckets["SNP"].append(v)
        elif svt == "DEL":
            buckets["DEL_s" if svlen < 25 else ("DEL_l" if svlen >= 1000 else "DEL_m")].append(v)
        elif svt == "INS":
            buckets["INS_s" if svlen < 25 else ("INS_l" if svlen >= 1000 else "INS_m")].append(v)
        elif svt == "DUP":
            if svlen >= 1000: buckets["DUP_l"].append(v)
        elif svt == "INV":
            buckets["INV"].append(v)

TAKE = {"SNP":5,"DEL_s":5,"DEL_l":5,"INS_s":5,"INS_l":5,"DUP_l":6,"INV":6}
picked = {k: buckets[k][:n] for k, n in TAKE.items()}

# batch-fetch all needed 0-based [a,b) regions via one samtools call; map by request order
regions = []                 # (chrom, a, b)
def need(chrom, a, b):
    a = max(0, a); regions.append((chrom, a, b)); return len(regions) - 1
plan = []                    # (label, meta, builder) builder = fn(fetched-list)->seq
def add_targets():
    for k, vs in picked.items():
        for v in vs:
            c, p, e = v["chrom"], v["pos0"], v["end0"]
            L = need(c, p - CTX, p)
            if k == "SNP":
                R = need(c, p + 1, p + CTX)
                plan.append((f"SNP_{v['svid']}_alt", dict(type="SNP",svid=v['svid'],svlen=0,role="seam",expect="c2"),
                             (lambda F,L=L,R=R,b=v['alt']: F[L] + b + F[R])))
            elif k.startswith("DEL"):
                R = need(c, e, e + CTX)
                plan.append((f"DEL_{v['svid']}_junc", dict(type="DEL",svid=v['svid'],svlen=v['svlen'],role="seam",expect="c2"),
                             (lambda F,L=L,R=R: F[L] + F[R])))
                if k == "DEL_l":
                    B = need(c, p, min(e, p + 1500))
                    plan.append((f"DEL_{v['svid']}_body", dict(type="DEL",svid=v['svid'],svlen=v['svlen'],role="delbody",expect="c1"),
                                 (lambda F,B=B: F[B])))
            elif k.startswith("INS"):
                R = need(c, p, p + CTX)
                plan.append((f"INS_{v['svid']}_alt", dict(type="INS",svid=v['svid'],svlen=v['svlen'],role="seam",expect="c2"),
                             (lambda F,L=L,R=R,ins=v['inseq']: F[L] + ins + F[R])))
            elif k == "DUP_l":
                JL = need(c, e - CTX, e); JR = need(c, p, p + CTX)
                plan.append((f"DUP_{v['svid']}_junc", dict(type="DUP",svid=v['svid'],svlen=v['svlen'],role="seam",expect="c2"),
                             (lambda F,JL=JL,JR=JR: F[JL] + F[JR])))
                B = need(c, p, min(e, p + 1500))
                plan.append((f"DUP_{v['svid']}_body", dict(type="DUP",svid=v['svid'],svlen=v['svlen'],role="body",expect="c3"),
                             (lambda F,B=B: F[B])))
            elif k == "INV":
                IE = need(c, e - CTX, e); IP = need(c, p, p + CTX); RR = need(c, e, e + CTX)
                plan.append((f"INV_{v['svid']}_left", dict(type="INV",svid=v['svid'],svlen=v['svlen'],role="seam",expect="c2"),
                             (lambda F,L=L,IE=IE: F[L] + rc(F[IE]))))
                plan.append((f"INV_{v['svid']}_right", dict(type="INV",svid=v['svid'],svlen=v['svlen'],role="seam",expect="c2"),
                             (lambda F,IP=IP,RR=RR: rc(F[IP]) + F[RR])))
                Bi = need(c, p, min(e, p + 1500))
                plan.append((f"INV_{v['svid']}_interior", dict(type="INV",svid=v['svid'],svlen=v['svlen'],role="interior",expect="c3"),
                             (lambda F,Bi=Bi: F[Bi])))
add_targets()

# one samtools faidx call for all regions (1-based inclusive: [a+1, b])
reg_args = [f"{c}:{a+1}-{b}" for (c, a, b) in regions]
res = subprocess.run([SAMTOOLS, "faidx", REF, *reg_args], capture_output=True, text=True)
if res.returncode != 0:
    sys.stderr.write(res.stderr); sys.exit(1)
# parse FASTA records IN ORDER -> list aligned with `regions`
F = []; cur = []
for ln in res.stdout.splitlines():
    if ln.startswith(">"):
        if cur: F.append("".join(cur)); cur = []
    else:
        cur.append(ln.strip())
if cur: F.append("".join(cur))
assert len(F) == len(regions), (len(F), len(regions))
F = [s.upper() for s in F]

os.makedirs(OUT, exist_ok=True)
targets = {}
with open(os.path.join(OUT, "verify_seqs.tsv"), "w") as out:
    for label, meta, builder in plan:
        seq = builder(F)
        if len(seq) < 60: continue
        out.write(f"{label}\t{seq}\n")
        targets[label] = meta
with open(os.path.join(OUT, "verify_targets.json"), "w") as out:
    json.dump(targets, out, indent=0)
print(f"wrote {len(targets)} verification sequences for {sum(len(v) for v in picked.values())} variants -> {OUT}/verify_seqs.tsv")
for k in TAKE: print(f"   {k:7s} sampled {len(picked[k])}")
