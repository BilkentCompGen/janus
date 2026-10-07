import argparse, bisect, random, re, subprocess, sys, zlib

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--fai", default="/home/yuio/snape/ref/hg38_parmask.fa.fai")
ap.add_argument("--vntr", default="/home/yuio/snape/ref/vntr.trf.bed")
ap.add_argument("--events", type=int, default=2000); ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--dominant", default="none"); ap.add_argument("--plan-only", action="store_true")
ap.add_argument("--region", default="")
A = ap.parse_args()
rng = random.Random(A.seed)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CG = re.compile(r"(\d+)([MIDNSHP=X])")
TYPES = ["DEL", "INS", "DUP", "INV", "BND"]
MARGIN = 300; REACH = 45000

def rc(s): return s[::-1].translate(COMP)
def half(name): return zlib.crc32(name.encode()) & 1
def draw(name, k): return zlib.crc32((name + "#" + str(k)).encode()) / 4294967296.0

class Read:
    __slots__ = ("name", "s", "e", "seq", "cg", "bl", "mq")
    def __init__(self, f):
        self.name = f[0]; self.s = int(f[3]) - 1; self.mq = int(f[4]); self.seq = f[9]; self.cg = f[5]; self.bl = None
        self.e = self.s + sum(int(n) for n, op in CG.findall(f[5]) if op in "M=XDN")
    @property
    def blocks(self):
        if self.bl is None:
            r, q, b = self.s, 0, []
            for n, op in CG.findall(self.cg):
                n = int(n)
                if op in "M=X": b.append((r, q, n, 1)); r += n; q += n
                elif op in "IS": q += n
                elif op in "DN": b.append((r, q, n, 0)); r += n
            self.bl = b
        return self.bl
    def off(self, x):
        for r, q, n, m in self.blocks:
            if x < r: return q
            if x < r + n: return q + (x - r if m else 0)
        return self.blocks[-1][1] + (self.blocks[-1][2] if self.blocks[-1][3] else 0)
    def qs(self): return self.blocks[0][1]
    def qe(self): b = self.blocks[-1]; return b[1] + (b[2] if b[3] else 0)

CACHE = {}
def region(c, a, b):
    k = (c, a // 20000, b // 20000)
    if k not in CACHE:
        if len(CACHE) > 8: CACHE.clear()
        o = subprocess.run(["samtools", "view", "-F", "0x904", A.bam, f"{c}:{max(1, a // 20000 * 20000)}-{(b // 20000 + 1) * 20000}"], capture_output=True, text=True).stdout
        CACHE[k] = [Read(f) for f in (l.split("\t", 10) for l in o.splitlines()) if zlib.crc32(f[0].encode()) & 1]
        CACHE[k] = [r for r in CACHE[k] if r.e > r.s and half(r.name) == 1]
    return CACHE[k]

def walk(c, p, d, n):
    out = []; got = 0
    for _ in range(5):
        if got >= n: break
        pool = [r for r in region(c, p - REACH, p + REACH) if r.mq >= 20 and (r.s + 50 <= p and r.e > p + 200 if d == "+" else r.s < p - 200 and r.e >= p + 50 - 50)]
        pool = [r for r in pool if (r.e - p if d == "+" else p - r.s) >= 500]
        if not pool: break
        r = max(rng.sample(pool, min(4, len(pool))), key=lambda r: r.e - p if d == "+" else p - r.s)
        if d == "+":
            q = r.off(p); take = r.seq[q:min(r.qe(), q + n - got)]; out.append(take); got += len(take); p = r.e
        else:
            q = r.off(p); lo = max(r.qs(), q - (n - got)); take = rc(r.seq[lo:q]); out.append(take); got += len(take); p = r.s
    return "".join(out)

def junction_read(r, own, other, ins, n):
    c, p, d = own
    if d == "-": return r.seq[:r.off(p)] + ins + walk(other[0], other[1], other[2], n)
    return rc(walk(other[0], other[1], other[2], n)) + rc(ins) + r.seq[r.off(p):]

def ends(ev):
    t, c, a, c2, b, ins = ev["type"], ev["c"], ev["a"], ev["c2"], ev["b"], ev.get("ins", "")
    if t == "DEL": return [((c, a, "-"), (c, b, "+"), "")]
    if t == "INS": return [((c, a, "-"), (c, a, "+"), ins)]
    if t == "DUP": return [((c, b, "-"), (c, a, "+"), "")]
    if t == "INV": return [((c, a, "-"), (c, b, "-"), ""), ((c, a, "+"), (c, b, "+"), "")]
    return [((c, a, ev["d1"]), (c2, b, ev["d2"]), "")]

def spans(r, p): return r.s + MARGIN <= p <= r.e - MARGIN

def apply(ev, k, edits):
    t, c, a, b, f = ev["type"], ev["c"], ev["a"], ev["b"], ev["f"]; made = 0
    same = t != "BND"
    if same:
        for r in region(c, a - REACH, b + REACH):
            if r.name in edits: continue
            u = draw(r.name, k); sa, sb = spans(r, a), spans(r, b)
            if t == "INS":
                if sa and u < f: edits[r.name] = r.seq[:r.off(a)] + ev["ins"] + r.seq[r.off(a):]; made += 1
                continue
            if sa and sb:
                if u < f:
                    qa, qb = r.off(a), r.off(b)
                    edits[r.name] = r.seq[:qa] + r.seq[qb:] if t == "DEL" else r.seq[:qb] + r.seq[qa:] if t == "DUP" else r.seq[:qa] + rc(r.seq[qa:qb]) + r.seq[qb:]
                    made += 1
                continue
            if t == "DEL" and r.s >= a - MARGIN and r.e <= b + MARGIN and not sa and not sb:
                if u < f: edits[r.name] = None
    if same and t == "INS": return made
    opts = {}
    for e1, e2, ins in ends(ev):
        opts.setdefault((e1[0], e1[1]), []).append((e1, e2, ins)); opts.setdefault((e2[0], e2[1]), []).append((e2, e1, rc(ins)))
    for (oc, op), ol in opts.items():
        for r in region(oc, op - REACH, op + REACH):
            if r.name in edits or not spans(r, op): continue
            if same and spans(r, a) and spans(r, b): continue
            u = draw(r.name, k)
            if u >= f: continue
            if len(ol) == 1 and u >= f / 2: edits[r.name] = None; continue
            own, other, sq = ol[zlib.crc32(r.name.encode()) >> 1 & 1] if len(ol) > 1 else ol[0]
            n = len(r.seq) - r.off(op) if own[2] == "-" else r.off(op)
            s = junction_read(r, own, other, sq, n)
            if len(s) >= 1000: edits[r.name] = s; made += 1
            else: edits[r.name] = None
    return made

def depth(c, p): return sum(1 for r in region(c, p - REACH, p + REACH) if spans(r, p) and r.mq >= 20)

def plan():
    chrom = [(l.split()[0], int(l.split()[1])) for l in open(A.fai)][:23]
    if A.region: chrom = [x for x in chrom if x[0] in A.region.split(",")]
    tot = sum(n for _, n in chrom); vn = []
    for l in open(A.vntr):
        c, x, y = l.split()[:3]
        if int(y) - int(x) >= 150 and any(c == q for q, _ in chrom): vn.append((c, int(x), int(y)))
    w = {t: 1.0 for t in TYPES}
    if A.dominant in w: w[A.dominant] = 6.0
    used = {}; evs = []
    def pos():
        x = rng.randrange(tot)
        for c, n in chrom:
            if x < n: return c, max(100000, min(n - 100000, x))
            x -= n
    def free(c, p):
        l = used.setdefault(c, []); i = bisect.bisect(l, p)
        return (i == 0 or p - l[i - 1] > 60000) and (i == len(l) or l[i] - p > 60000)
    def take(c, p): bisect.insort(used.setdefault(c, []), p)
    tries = 0
    while len(evs) < A.events and tries < A.events * 30:
        tries += 1
        t = rng.choices(TYPES, [w[x] for x in TYPES])[0]
        f = rng.choices([0.5, 0.35, 0.2, 0.1, 1.0], [5, 2, 2, 1, 1])[0]
        if t in ("DEL", "DUP", "INS") and rng.random() < 0.3:
            c, x, y = rng.choice(vn); n = min(y - x, int(50 * 10 ** (rng.random() * 1.5)))
            a = rng.randrange(x, y - n + 1); b = a + n
        else:
            c, a = pos(); n = int(50 * 10 ** (rng.random() * (3.0 if t != "INV" else 2.7))) if t != "INV" else int(300 * 10 ** (rng.random() * 2.2)); b = a + n
        ev = dict(type=t, c=c, a=a, c2=c, b=b, f=f)
        if t == "INS":
            ev["b"] = a; n = min(n, 6000); k = rng.random()
            if k < 0.4: ev["ins"] = "".join(rng.choice("ACGT") for _ in range(n))
            else:
                c3, p3 = pos(); s = walk(c3, p3, "+", n)
                if len(s) < n: continue
                ev["ins"] = s
            ev["len"] = n
        elif t == "BND":
            c2, b = pos()
            if c2 == c: continue
            ev.update(c2=c2, b=b, d1=rng.choice("+-"), d2=rng.choice("+-"), len=0)
        else: ev["len"] = n
        if not free(c, a) or (t == "BND" and not free(ev["c2"], b)): continue
        if depth(c, a) < 8 or depth(ev["c2"], ev["b"]) < 8: continue
        take(c, a); take(ev["c2"], ev["b"]); evs.append(ev)
    return evs

evs = plan(); edits = {}
with open(A.out + ".truth.vcf", "w") as o:
    o.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for k, ev in enumerate(evs):
        ev["made"] = apply(ev, k, edits)
        if ev["type"] == "BND":
            br = "[" if ev["d2"] == "+" else "]"; alt = f"N{br}{ev['c2']}:{ev['b'] + 1}{br}" if ev["d1"] == "-" else f"{br}{ev['c2']}:{ev['b'] + 1}{br}N"
            info = "SVTYPE=BND"
        else:
            alt = f"<{ev['type']}>"; info = f"SVTYPE={ev['type']};END={ev['b'] + 1};SVLEN={ev['len']}"
        o.write(f"{ev['c']}\t{ev['a'] + 1}\tsim{k}\tN\t{alt}\t.\t{'PASS' if ev['made'] >= 1 else 'NOREADS'}\t{info};AF={ev['f']};READS={ev['made']}\n")
sys.stderr.write(f"[implant] {len(evs)} events, {sum(1 for e in evs if e['made'] >= 3)} with >= 3 junction reads, {len(edits)} reads edited or removed\n")
if A.plan_only: sys.exit(0)
HD = "@HD\tVN:1.6\tSO:unsorted\n"
P = [subprocess.Popen(["samtools", "view", "-b", "-@", "3", "--output-fmt-option", "level=1", "-o", f"{A.out}.{x}.bam", "-"], stdin=subprocess.PIPE) for x in ("n", "t")]
for p in P: p.stdin.write(HD.encode())
src = subprocess.Popen(["samtools", "view", "-@", "4", "-F", "0x900", A.bam] + (A.region.split(",") if A.region else []), stdout=subprocess.PIPE)
n = [0, 0, 0, 0]
for l in src.stdout:
    f = l.split(b"\t", 10); name = f[0].decode(); h = half(name); s = f[9]
    if h and name in edits:
        e = edits[name]
        if e is None: n[2] += 1; continue
        s = e.encode(); n[3] += 1
    P[h].stdin.write(f[0] + b"\t4\t*\t0\t0\t*\t*\t0\t0\t" + s + b"\t*\n"); n[h] += 1
for p in P: p.stdin.close(); p.wait()
sys.stderr.write(f"[implant] normal half {n[0]} reads, tumor half {n[1]} reads ({n[3]} edited, {n[2]} removed)\n")
