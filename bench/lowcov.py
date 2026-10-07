import collections, mmap, re, struct, subprocess, sys
S = "/home/yuio/snape"
p, D, vcf = sys.argv[1], sys.argv[2], sys.argv[3]
exec(open(f"{S}/bench/score.py").read().split("T, C = recs")[0])
T = recs(f"{S}/hifi/truth/castle_without_snape/{p}/None_minda_ensemble.vcf")
ENC = {"A": 0, "C": 1, "G": 2}
COMP = str.maketrans("ACGTN", "TGCAN")

def pack(s):
    b = bytearray(8)
    for i, c in enumerate(s[:32]): b[i >> 2] |= ENC.get(c, 3) << (6 - 2 * (i & 3))
    return bytes(b)

def keys(path):
    f = open(path, "rb"); m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    n = struct.unpack_from("<Q", m, 8)[0]; at = 32; out = {}
    for i in range(n):
        tag = struct.unpack_from("<I", m, at)[0]; ln = tag & 0x7fffffff
        if ln >= 32:
            k = pack(m[at + 4:at + 36].decode().upper()) if tag >> 31 else bytes(m[at + 4:at + 12])
            out[(ln, k)] = i
        at += 4 + (ln if tag >> 31 else (ln + 3) // 4)
    return out

KT, KN = keys(f"{D}/t.svb"), keys(f"{D}/n.svb")
CG = re.compile(r"(\d+)([MIDNSHP=X])")

def reads(bam, c, pos, K, w=300):
    sig, span = set(), set()
    o = subprocess.run(["samtools", "view", bam, f"{c}:{max(1, pos - w)}-{pos + w}"], capture_output=True, text=True).stdout
    for l in o.splitlines():
        f = l.split("\t", 11); flag = int(f[1]); s = f[9]
        if flag & 0x100 or len(s) < 32: continue
        k = (len(s), pack(s[-32:][::-1].translate(COMP)) if flag & 16 else pack(s[:32]))
        i = K.get(k)
        if i is None: continue
        r = int(f[3]); ops = CG.findall(f[5]); hit = False
        if ops and ops[0][1] in "SH" and int(ops[0][0]) >= 200 and abs(r - pos) <= w: hit = True
        for n, op in ops:
            n = int(n)
            if op in "ID" and n >= 40 and r - w <= pos <= r + (n if op == "D" else 0) + w: hit = True
            if op in "MD=XN": r += n
        if ops and ops[-1][1] in "SH" and int(ops[-1][0]) >= 200 and abs(r - pos) <= w: hit = True
        span.add(i)
        if hit: sig.add(i)
    return sig, span

R = []
for l in open(vcf):
    if l[0] == "#": continue
    f = l.rstrip("\n").split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
    R.append(((info["SVTYPE"], f[0], int(f[1]), m.group(1), int(m.group(2))) if m else (info["SVTYPE"], f[0], int(f[1]), f[0], int(info.get("END", f[1]))), f[6]))
ix = collections.defaultdict(list)
for r, flt in R: ix[r[1]].append((r, flt)); ix[r[3]].append((r, flt))
aln = collections.defaultdict(list)
for l in open(f"{D}/ev.paf"):
    f = l.split("\t", 9); aln[f[5]].append((int(f[7]), int(f[8]), f[0].split("|")[0]))
pieces = collections.defaultdict(set)
for l in open(f"{D}/ev.fa"):
    if l[0] == ">":
        s, i, k = l[1:].split()[0].split("_"); pieces[int(i)].add(l[1:].split()[0])
placed = set(l[1:].split()[0] for l in open(f"{D}/ev.on.fa") if l[0] == ">")
sites = collections.Counter(l.split("\t")[0] for l in open(f"{D}/ev.tsv"))

def near(c, pos, w=500):
    return set(q for a, b, q in aln.get(c, []) if a - w <= pos <= b + w)

out = open(f"{D}/lowcov.tsv", "w")
out.write("type\tchr1\tpos1\tchr2\tpos2\toutcome\ttumor_sig\ttumor_span\tnormal_sig\tstage\n")
tab = collections.Counter()
for t in T:
    h = [flt for r, flt in ix.get(t[1], []) + ix.get(t[3], []) if hit(t, r)]
    res = "PASS" if "PASS" in h else ("LOW" if h else "none")
    s1, sp1 = reads(f"{S}/hifi/{p}/t.bam", t[1], t[2], KT); s2, sp2 = reads(f"{S}/hifi/{p}/t.bam", t[3], t[4], KT)
    n1, _ = reads(f"{S}/hifi/{p}/n.bam", t[1], t[2], KN); n2, _ = reads(f"{S}/hifi/{p}/n.bam", t[3], t[4], KN)
    sig = s1 | s2; nsig = n1 | n2; span = sp1 | sp2
    stage = "-"
    if res == "none":
        cov = near(t[1], t[2]) | near(t[3], t[4])
        own = set(q for i in sig for q in pieces.get(i, ()))
        if len(sig) < 3: stage = "support<3"
        elif own & cov: stage = "call: piece of a supporting read aligned here, no record"
        elif cov: stage = "call/select: other piece aligned here, none from supporting reads" if not own else "place: supporting read's piece not aligned here; other piece is"
        elif own: stage = "place/align: supporting read has a piece" + (" (placed elsewhere)" if own & placed else " (unplaced)")
        else: stage = "no site on any supporting read (graph/prune)"
    b = "0-2" if len(sig) < 3 else "3-4" if len(sig) < 5 else "5-9" if len(sig) < 10 else "10+"
    tab[(b, res if res != "none" else "none: " + stage)] += 1
    out.write(f"{t[0]}\t{t[1]}\t{t[2]}\t{t[3]}\t{t[4]}\t{res}\t{len(sig)}\t{len(span)}\t{len(nsig)}\t{stage}\n")
out.close()
print(p, "truth", len(T))
for b in ("0-2", "3-4", "5-9", "10+"):
    tot = sum(v for (x, _), v in tab.items() if x == b)
    print(f"tumor signature reads {b}: {tot} events")
    for (x, k), v in sorted(tab.items(), key=lambda z: -z[1]):
        if x == b: print(f"    {v:5d}  {k}")
