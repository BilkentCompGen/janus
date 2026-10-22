import sys, io, contextlib
SCRIPTS = {}
SCRIPTS["measure"] = r'''#!/usr/bin/env python3
import sys, os, re, json, glob, datetime, subprocess, platform

D = os.path.abspath(sys.argv[1]); W = os.path.join(D, "work"); M = os.path.join(D, "monitor")
for f in ("start.txt", "total.txt", "note.txt"):
    if os.path.exists(os.path.join(D, f)): os.replace(os.path.join(D, f), os.path.join(M, f))
KEEP = {"calls.vcf"}

def iostat():
    out, t, first = [], None, True
    for l in open(os.path.join(M, "iostat.txt")):
        f = l.split()
        if len(f) == 3 and f[2] in ("AM", "PM") and "/" in f[0]:
            t = datetime.datetime.strptime(" ".join(f), "%m/%d/%Y %I:%M:%S %p").timestamp()
        elif f and f[0].startswith("nvme") and t is not None:
            if first: first = False; continue
            out.append((t, float(f[2]), float(f[8]), float(f[-1])))
    return out

def temps():
    out = []
    for l in open(os.path.join(M, "temp.txt")):
        f = l.split()
        if len(f) < 3: continue
        c = float(f[3].strip("+°C")) if len(f) > 3 else None
        out.append((int(f[0]), int(f[1]) / 1000, int(f[2]) / 1000, c))
    return out

IO, TP = iostat(), temps()

def window(a, b):
    io = [r for r in IO if a < r[0] <= b + 5]
    tp = [r for r in TP if a <= r[0] <= b + 5]
    if not io: rd = wr = ut = None
    else: rd = round(sum(r[1] for r in io) * 5 / 1000, 2); wr = round(sum(r[2] for r in io) * 5 / 1000, 2); ut = round(sum(r[3] for r in io) / len(io), 1)
    mx = lambda i: max((r[i] for r in tp if r[i] is not None), default=None)
    return rd, wr, ut, mx(1), mx(2), mx(3)

def report(n):
    for f in (os.path.join(W, "log." + n), os.path.join(D, "logs", "log." + n)):
        if os.path.exists(f):
            return [l.strip() for l in open(f, errors="replace") if l.startswith("[") and not l.startswith("[M::")]
    return []

steps = []
for l in open(os.path.join(W, "steps.txt")):
    f = l.split()
    if len(f) < 8 or f[2] != "s": continue
    n = f[0]; e = os.path.join(W, "done." + n)
    end = os.path.getmtime(e if os.path.exists(e) else os.path.join(W, "log." + n)); a = end - float(f[1])
    s = dict(step=n, start=round(a, 2), end=round(end, 2), wall_s=float(f[1]), peak_rss_kib=int(f[3]),
             user_s=float(f[5]), sys_s=float(f[7]))
    s["cpu_s"] = round(s["user_s"] + s["sys_s"], 2)
    s["read_gb"], s["written_gb"], s["util_pct"], s["drive_c"], s["controller_c"], s["cpu_c"] = window(a, end)
    s["report"] = report(n)
    steps.append(s)
prev = steps[-1]["end"]

if os.path.exists(os.path.join(M, "total.txt")):
    start = float(open(os.path.join(M, "start.txt")).read())
    tot = open(os.path.join(M, "total.txt")).read().split()
    total = dict(wall_s=float(tot[1]), peak_rss_kib=int(tot[4]), user_s=float(tot[7]), sys_s=float(tot[10]))
else:
    start = steps[0]["start"]
    total = dict(wall_s=round(sum(s["wall_s"] for s in steps), 2), peak_rss_kib=max(s["peak_rss_kib"] for s in steps),
                 user_s=round(sum(s["user_s"] for s in steps), 2), sys_s=round(sum(s["sys_s"] for s in steps), 2))
total["cpu_s"] = round(total["user_s"] + total["sys_s"], 2)
total["step_wall_s"] = round(sum(s["wall_s"] for s in steps), 2)
total["peak_step"] = max(steps, key=lambda s: s["peak_rss_kib"])["step"]
total["read_gb"], total["written_gb"], total["util_pct"], total["drive_c"], total["controller_c"], total["cpu_c"] = window(start, prev)

files = sorted(((os.path.getsize(p), os.path.relpath(p, W)) for p in glob.glob(os.path.join(W, "*")) if os.path.isfile(p)), reverse=True)
files = [dict(file=f, bytes=b, kept=f in KEEP or f.startswith("log.")) for b, f in files if not f.startswith("done.")]
GONE = "deleted" if os.path.exists(os.path.join(M, "params.json")) and json.load(open(os.path.join(M, "params.json"))).get("intermediates") == "deleted" else "trash"
total["intermediate_gb"] = round(sum(f["bytes"] for f in files) / 1e9, 2)

cols = ["step", "start", "end", "wall_s", "peak_rss_kib", "user_s", "sys_s", "cpu_s", "read_gb", "written_gb", "util_pct",
        "drive_c", "controller_c", "cpu_c"]
with open(os.path.join(D, "steps.tsv"), "w") as o:
    o.write("\t".join(cols) + "\n")
    for s in steps + [dict(total, step="total", start=round(start, 2), end=round(prev, 2))]:
        o.write("\t".join("" if s.get(c) is None else str(s[c]) for c in cols) + "\n")
with open(os.path.join(D, "files.tsv"), "w") as o:
    o.write("file\tbytes\tkept\n")
    for f in files: o.write(f"{f['file']}\t{f['bytes']}\t{'yes' if f['kept'] else GONE}\n")

S = "/home/yuio/janus"
P = json.load(open(os.path.join(M, "params.json"))) if os.path.exists(os.path.join(M, "params.json")) else {}
pair = P.get("pair", os.path.basename(D))
src = lambda s: next(f for f in (f"{S}/hifi/{pair}/{s}.bam", f"{S}/hifi/{pair}/{s}.fastq.gz") if os.path.exists(f))
reads = {s: dict(file=os.path.realpath(src(s)), bytes=os.path.getsize(src(s))) for s in "nt"}
for st in steps:
    if st["step"] in ("pack_n", "pack_t"):
        for l in st["report"]:
            m = re.match(r"\[pack\] (\d+) reads, (\d+) bases", l)
            if m: reads[st["step"][-1]].update(reads=int(m.group(1)), bases=int(m.group(2)))
cpu = next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")), "")
mem = int(next(l.split()[1] for l in open("/proc/meminfo") if l.startswith("MemTotal")))
commit = subprocess.run(["git", "-C", S, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
note = open(os.path.join(M, "note.txt")).read().strip() if os.path.exists(os.path.join(M, "note.txt")) else ""
json.dump(dict(pair=pair, run=os.path.basename(D), params=P, commit=commit, note=note, started=datetime.datetime.fromtimestamp(start).isoformat(timespec="seconds"),
               threads=P.get("threads", 16), cores="taskset 0-%d" % (P.get("threads", 16) - 1), machine=dict(cpu=cpu, mem_kib=mem, kernel=platform.release()),
               reads=reads, total=total, steps=steps, files=files),
          open(os.path.join(D, "run.json"), "w"), indent=1)
print(f"{pair}: {total['wall_s']} s wall, {total['cpu_s']} CPU s, peak {total['peak_rss_kib'] / 2**20:.1f} GiB ({total['peak_step']}), "
      f"intermediates {total['intermediate_gb']} GB")
'''
SCRIPTS["tool"] = r'''#!/usr/bin/env python3
import sys, os, json
D = os.path.abspath(sys.argv[1])
P = json.load(open(f"{D}/params.json"))
wall, rss, us, sy, ex = open(f"{D}/time.txt").read().split()[-5:]
t = dict(wall_s=float(wall), host_peak_kib=int(rss), host_user_s=float(us), host_sys_s=float(sy), exit=int(ex))
cg = {}
if os.path.exists(f"{D}/cgroup.txt"):
    for l in open(f"{D}/cgroup.txt"):
        k, v = l.split()
        cg[k] = int(v)
t["docker_user_s"] = cg.get("user_usec", 0) / 1e6
t["docker_sys_s"] = cg.get("system_usec", 0) / 1e6
t["docker_peak_kib"] = cg.get("memory_peak", 0) // 1024
t["user_s"] = round(t["host_user_s"] + t["docker_user_s"], 2)
t["sys_s"] = round(t["host_sys_s"] + t["docker_sys_s"], 2)
t["cpu_s"] = round(t["user_s"] + t["sys_s"], 2)
t["docker_anon_peak_kib"] = cg.get("anon_peak", 0) // 1024
t["peak_rss_kib"] = max(t["host_peak_kib"], t["docker_anon_peak_kib"] if "anon_peak" in cg else t["docker_peak_kib"])
files = []
O = f"{D}/out"
for root, _, fs in os.walk(O):
    for f in fs:
        p = os.path.join(root, f)
        if os.path.isfile(p) and not os.path.islink(p):
            files.append(dict(file=os.path.relpath(p, O), bytes=os.path.getsize(p), kept=os.path.relpath(p, O) in ("calls.vcf", "log.txt")))
for k in ("bam", "bam_n", "bam_t"):
    if k in P and os.path.exists(P[k]): files.append(dict(file=os.path.basename(P[k]), bytes=os.path.getsize(P[k]), kept=False))
files.sort(key=lambda f: -f["bytes"])
calls = None
if os.path.exists(f"{O}/calls.vcf"):
    calls = sum(1 for l in open(f"{O}/calls.vcf") if not l.startswith("#"))
json.dump(dict(params=P, total=t, calls=calls, files=files), open(f"{D}/run.json", "w"), indent=1)
with open(f"{D}/files.tsv", "w") as o:
    o.write("file\tbytes\tkept\n")
    for f in files: o.write(f"{f['file']}\t{f['bytes']}\t{'yes' if f['kept'] else 'deleted'}\n")
print(f"{P.get('tool')} {P.get('pair')} {P.get('coverage_x') or 'full'}: exit {t['exit']}, {t['wall_s']:.0f} s wall, "
      f"{t['cpu_s']:.0f} CPU s, peak {t['peak_rss_kib'] / 2**20:.1f} GiB, calls {calls}")
'''
SCRIPTS["mix"] = r'''import mmap, struct, sys
out = sys.argv[1]; parts = [(a.split(":")[0], int(a.split(":")[1]), int(a.split(":")[2])) for a in sys.argv[2:]]
n = total = longest = 0; off = []
with open(out, "wb") as o:
    o.write(b"\0" * 32); pos = 32
    for src, skip, take in parts:
        with open(src, "rb") as f:
            m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            p, seen, got, end = 32, 0, 0, len(m)
            while p < end and got < take:
                tag = struct.unpack_from("<I", m, p)[0]; ln = tag & 0x7fffffff
                sz = 4 + (ln if tag >> 31 else (ln + 3) // 4)
                if seen >= skip:
                    o.write(m[p:p + sz]); off.append(pos); pos += sz; got += ln; total += ln; n += 1; longest = max(longest, ln)
                seen += ln; p += sz
            if got < take: sys.exit(f"{src}: only {got} of {take} bases after skipping {skip}")
    off.append(pos)
    o.seek(0); o.write(b"SVB1" + struct.pack("<IQQ", 1, n, total) + b"\0" * 8)
with open(out + ".idx", "wb") as o:
    o.write(b"SVBI" + struct.pack("<IQQ", 1, n, longest) + b"\0" * 8)
    o.write(struct.pack("<%dQ" % len(off), *off))
print(out, n, "reads", total, "bases")
'''
SCRIPTS["score"] = r'''import gzip, re, sys
TOL = 500
ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")

def recs(path):
    op = gzip.open if path.endswith(".gz") else open
    out, mates = [], set()
    for ln in op(path, "rt"):
        if ln[0] == "#": continue
        f = ln.rstrip("\n").split("\t")
        if f[6] not in ("PASS", "."): continue
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        t, c, p = info.get("SVTYPE"), f[0], int(float(f[1]))
        m = ALT.search(f[4])
        if m:
            if info.get("MATEID", info.get("MATE_ID")) in mates: continue
            mates.add(f[2])
            out.append((t, c, p, m.group(1), int(m.group(2))))
        else:
            out.append((t, c, p, c, int(float(info.get("END", p)))))
    return out

def hit(a, b):
    for x in ((b[1], b[2], b[3], b[4]), (b[3], b[4], b[1], b[2])):
        if a[1] == x[0] and a[3] == x[2] and abs(a[2] - x[1]) <= TOL and abs(a[4] - x[3]) <= TOL: return True
    if a[0] == "INS" or b[0] == "INS":
        return a[1] == b[1] and abs(a[2] - b[2]) <= TOL
    return False

T, C = recs(sys.argv[1]), recs(sys.argv[2])
tp_t = sum(any(hit(t, c) for c in C) for t in T)
tp_c = sum(any(hit(c, t) for t in T) for c in C)
r, p = tp_t / len(T), (tp_c / len(C) if C else 0)
print("truth %d calls %d recall %.3f precision %.3f f1 %.3f" % (len(T), len(C), r, p, 2 * r * p / (r + p) if r + p else 0))
'''
SCRIPTS["verdict"] = r'''import sys, time, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
exec(SCRIPTS["score"].split("T, C = recs")[0])
sys.path.insert(0, "/home/yuio/janus/src")
from janus import features as JF, classify as JC
H = "/home/yuio/janus/bench"; TR = "/home/yuio/janus/hifi/truth"
TRUTH = {p: f"{H}/eval_sets/{p}/all_without/None_minda_ensemble.vcf" for p in ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]}
TRUTH.update(hg008=f"{TR}/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz",
             colo829=f"{TR}/colo829/truth.vcf")
RUNS = __import__("os").environ.get("RUNS", f"{H}/runs/full"); OUTV = __import__("os").environ.get("OUTV", f"{H}/verdict")
TRAIN = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]
C = __import__("os").environ.get("COV", f"{H}/runs/coverage"); U = __import__("os").environ.get("PUR", f"{H}/runs/purity"); O = f"{OUTV}/coverage_matched"
DEPTHS = [5, 10, 15, 20, 30]; PUR = [20, 40, 60, 80, 100]
HGB = dict(JC.HGB); CUTS = JC.CUTS
TT = {p: recs(TRUTH[p]) for p in TRAIN}
BY = {}
for p in TRAIN:
    BY[p] = {}
    for t in TT[p]: BY[p].setdefault(t[1], []).append(t); BY[p].setdefault(t[3], []).append(t)

def crecords(path, p):
    L, X, K, R = [], [], [], []
    for l in open(path):
        if l[0] == "#": continue
        f = l.split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
        K.append((f[0], m.group(1) if m else f[0], int(f[1]), int(m.group(2)) if m else int(f[1]) + int(info["SVLEN"]), info["SVTYPE"], info["O"]))
        R.append((info["SVTYPE"], f[0], int(f[1]), m.group(1), int(m.group(2))) if m else (info["SVTYPE"], f[0], int(f[1]), f[0], int(info["END"])))
        X.append([float(v) for v in info["F"].split(",")]); L.append(l)
    y = np.array([any(hit(r, t) for t in BY[p].get(r[1], [])) for r in R], float)
    keep = np.array(["\tDUPL\t" not in l for l in L])
    return L, np.array(X), K, y, keep

def score(v, p):
    Cc = recs(v); bc = {}
    for c in Cc: bc.setdefault(c[1], []).append(c); bc.setdefault(c[3], []).append(c)
    tpc = sum(any(hit(c, t) for t in BY[p].get(c[1], [])) for c in Cc)
    tpt = sum(any(hit(t, c) for c in bc.get(t[1], []) + bc.get(t[3], [])) for t in TT[p])
    P = tpc / max(len(Cc), 1); R = tpt / max(len(TT[p]), 1)
    return len(Cc), tpc, tpt, P, R, 2 * P * R / max(P + R, 1e-9)

mark, f1, pct = JF.duplicates, JC.f1, JF.pct

t0 = time.time()
import os, json, pickle
os.makedirs(O, exist_ok=True); os.makedirs(f"{OUTV}/purity_matched", exist_ok=True)
FULL = [p for p in TRUTH if os.path.exists(f"{RUNS}/{p}/calls.vcf")]
for p in FULL:
    if p in TT: continue
    TT[p] = recs(TRUTH[p]); BY[p] = {}
    for t in TT[p]: BY[p].setdefault(t[1], []).append(t); BY[p].setdefault(t[3], []).append(t)

def sets(p):
    s = {"full": f"{RUNS}/{p}/calls.vcf"}
    s.update({f"{n}x": f"{C}/{p}_{n}x/calls.vcf" for n in DEPTHS}); s.update({f"{q}p": f"{U}/{p}_{q}p/calls.vcf" for q in PUR})
    return {k: f for k, f in s.items() if os.path.exists(f)}

def feats(path, p):
    L, X, K, y, keep = crecords(path, p)
    return L, JF.features(X, f"{os.path.dirname(path)}/t.cnt.hist", JF.minimum(X)), K, y, keep

def fitm(X, y):  # the bench varies HGB and fitm (supp); the default is JC.fit_models
    return [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(5)]

prob = JC.predict

def fit(train):
    ks = [k for k in D if k[0] in train and k[1] != "100p"]
    m = JC.train([(k[0], D[k][1][D[k][4]], D[k][3][D[k][4]]) for k in ks], fitm=lambda X, y: fitm(X, y))
    return m["models"], m["cuts"], m["cut"], m["train_s"], m["n"]

D = {(p, k): feats(f, p) for p in FULL for k, f in sets(p).items()}
log = open(f"{OUTV}/log.txt", "w"); folds = {}; fits = {}
cov = open(f"{O}/scores.tsv", "w"); cov.write("pair\tdepth\tcut\ttrain_pairs\tpass\ttruth\ttp_calls\ttp_truth\tprecision\trecall\tf1\n")
pur = open(f"{OUTV}/purity_matched/scores.tsv", "w"); pur.write("pair\tpurity\tcut\ttrain_pairs\tpass\ttruth\ttp_calls\ttp_truth\tprecision\trecall\tf1\n")
for p in FULL:
    tr = tuple(q for q in TRAIN if q != p)
    if tr not in fits: fits[tr] = fit(tr)
    ms, cuts, c0, ts, nt = fits[tr]
    if len(tr) == len(TRAIN): pickle.dump(dict(models=ms, cuts=cuts, cut=c0, train=tr), open(f"{OUTV}/model.pkl", "wb"))
    for q, k in [k for k in D if k[0] == p]:
        L, X, K, y, keep = D[(q, k)]; m = int(round(X[0, 56])); cut = cuts.get(m, c0)
        t = time.time(); pr = prob(ms, X); ti = time.time() - t; d = mark(K, pr)
        v = f"{OUTV}/{p}.vcf" if k == "full" else f"{O}/{p}_{k}.vcf" if k[-1] == "x" else f"{OUTV}/purity_matched/{p}_{k}.vcf"
        with open(v, "w") as o:
            o.write("".join(l for l in open(sets(p)[k]) if l[0] == "#"))
            for l, r, u in zip(L, pr, d):
                g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if r >= cut else "LOW"; g[7] += ";VP=%.4f" % r; o.write("\t".join(g))
        if k == "full":
            folds[p] = dict(train=list(tr), cut=round(float(cut), 2), train_s=round(ts, 3), train_records=nt, infer_s=round(ti, 4), records=len(L),
                            prob=[round(float(x), 4) for x in pr], label=[int(x) for x in y], dupl=[int(x) for x in d], features=len(X[0]))
            s = f"{p} {'held-out' if p in TRAIN else 'not in training'} PASS {int(((pr >= cut) & ~d).sum())} cut {cut:.2f} minimum {m} DUPL {int(d.sum())} records {len(L)}"
            print(s); log.write(s + "\n")
        else:
            n, tpc, tpt, pp, rr, ff = score(v, p)
            row = f"{cut:.2f}\t{','.join(tr)}\t{n}\t{len(TT[p])}\t{tpc}\t{tpt}\t{pp:.4f}\t{rr:.4f}\t{ff:.4f}\n"
            if k[-1] == "x": cov.write(f"{p}\t{k[:-1]}\t{row}")
            else: pur.write(f"{p}\t{k[:-1]}\t{row}")
s = f"verdict (one model for all depths and purities on depth-normalised counts and four count-stage context values, threshold per chosen minimum from held-out predictions inside the training pairs; each CASTLE pair held out; {len(FULL)} pairs, 5 seeds) {time.time() - t0:.1f} s"
print(s); log.write(s + "\n")
json.dump(folds, open(f"{OUTV}/folds.json", "w"))
'''
SCRIPTS["sweep"] = r'''import os, sys, glob, itertools, hashlib, numpy as np
H = "/home/yuio/janus/bench"; R = f"{H}/runs"; W = f"{R}/sweep"
CASTLE = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]
N0 = [5, 10, 15, 20, 25]; NT = [20, 30, 40, 50, 60]; GRID = list(itertools.product(N0, NT)) + [(70, 40), (10, 5), (70, 5)]

def exps(p):
    e = [("full", p)] + [("coverage", f"{p}_{n}x") for n in (5, 10, 15, 20, 30)] + [("purity", f"{p}_{q}p") for q in (20, 40, 60, 80, 100)]
    return [(k, n) for k, n in e if os.path.exists(f"{R}/{k}/{n}/t.cnt.hist")]

def hist(d):
    return {l.split()[0]: np.array(l.split()[1:], float) for l in open(f"{d}/t.cnt.hist")}

def median(h):
    n = h["shared_normal"] + h["normal_only"]; tot = n[3:].sum(); run = 0; m = 3
    while m + 1 < len(n) and 2 * (run + n[m]) < tot: run += n[m]; m += 1
    return m

def p95(h):
    t = h["tumor_only"]; tot = t[3:].sum(); run = 0; p = 3
    while p + 1 < len(t) and 20 * (run + t[p]) < 19 * tot: run += t[p]; p += 1
    return p

def pair(d):
    return open(f"{d}/n.cnt", "rb").read(27)[26], open(f"{d}/t.cnt", "rb").read(27)[26]

def runs(k, n):
    out = {}
    for d in [f"{R}/{k}/{n}"] + sorted(glob.glob(f"{W}/runs/{k}/{n}/*/")):
        d = d.rstrip("/")
        if all(os.path.exists(f"{d}/{f}") for f in ("calls.vcf", "n.cnt", "t.cnt", "t.cnt.hist")): out.setdefault(pair(d), d)
    return out

def choice(k, n, a, b):
    rs = runs(k, n); hs = {}
    for (m0, mt), d in rs.items(): hs.setdefault(m0, hist(d))
    m0 = 2 if median(next(iter(hs.values()))) < a else 3
    if m0 not in hs: return m0, 0, None
    mt = 3 if p95(hs[m0]) >= b else 2
    return m0, mt, rs.get((m0, mt))

mode = sys.argv[1]
if mode == "plan":
    for k, n in exps(sys.argv[2]):
        want = {}
        for a, b in GRID:
            m0, mt, d = choice(k, n, a, b)
            if not d: want.setdefault((m0, mt), (a, b))
        for (m0, mt), (a, b) in sorted(want.items()): print(k, n, a, b, m0, mt or "?")
elif mode == "check":
    k, n, d = sys.argv[2:5]; m0, mt = pair(d)
    ref = [x for (a, b), x in runs(k, n).items() if a == m0 and x != d]
    same = [all((hist(d)[g] == hist(x)[g]).all() for g in hist(d)) for x in ref]
    print(f"{k}/{n} {os.path.basename(d)}: minimum read support {m0}/{mt}" + (f", count spectra {'identical to' if all(same) else 'DIFFERENT FROM'} the earlier run with m0 {m0}" if ref else ""))
elif mode == "configs":
    os.makedirs(W, exist_ok=True); keys = {}; P = open(f"{W}/pairs.tsv", "w"); P.write("theta0\tthetat\texperiment\tm0\tmt\trun\n")
    for a, b in GRID:
        sel = {}
        for p in CASTLE:
            for k, n in exps(p):
                m0, mt, d = choice(k, n, a, b); sel[f"{k}/{n}"] = d; P.write(f"{a}\t{b}\t{k}/{n}\t{m0}\t{mt}\t{d}\n")
        if not all(sel.values()): print(f"theta {a}/{b}: missing runs for {[e for e, d in sel.items() if not d]}", file=sys.stderr); continue
        key = hashlib.sha1("".join(f"{e}={d};" for e, d in sorted(sel.items())).encode()).hexdigest()[:10]; keys[(a, b)] = key
        for e, d in sel.items():
            l = f"{W}/conf/{key}/{e}"; os.makedirs(os.path.dirname(l), exist_ok=True)
            if not os.path.islink(l): os.symlink(d, l)
    P.close()
    with open(f"{W}/combos.tsv", "w") as o:
        o.write("theta0\ttheta_t\tconfig\n"); o.write("".join(f"{a}\t{b}\t{k}\n" for (a, b), k in keys.items()))
    print("\n".join(sorted(set(keys.values()))))
elif mode == "rescore":
    C = sys.argv[2]; S = {}
    WITH = ["Severus", "SAVANA", "nanomonsv", "Sniffles2", "dysgu", "SVision-pro", "cuteSV", "DeBreak", "kled", "SVDSS", "janus", "Severus ONT", "SAVANA ONT", "nanomonsv ONT", "Sniffles2 ONT", "SvABA", "GRIPSS", "Manta"]
    LOO = {"snape": "janus", "severus_phased": "Severus", "savana_phased": "SAVANA", "nanomonsv": "nanomonsv", "dysgu": "dysgu", "svision-pro": "SVision-pro", "sniffles2": "Sniffles2", "cutesv": "cuteSV", "debreak": "DeBreak", "kled": "kled", "svdss": "SVDSS"}
    TOOL = {"snape": "janus", "severus": "Severus", "savana": "SAVANA", "nanomonsv": "nanomonsv"}
    SETS = ["with", "without", "loo", "30x", "20x", "15x", "10x", "5x", "100p", "80p", "60p", "40p", "20p"]
    rows = lambda d: [l.split("\t") for l in open(f"{d}/results/overall.tsv").read().split("\n")[1:] if l]
    for p in CASTLE:
        for s, names in (("with", WITH), ("without", [n for n in WITH if n != "janus"])):
            r = rows(f"{C}/eval_sets/{p}/all_{s}"); assert len(r) == len(names), (p, s)
            for n, x in zip(names, r): S[(s, n, p)] = [float(v) for v in x[4:7]]
        S[("without", "janus", p)] = [float(v) for v in rows(f"{C}/eval_sets/{p}/all_without_snape")[0][4:7]]
        for k, n in LOO.items(): S[("loo", n, p)] = [float(v) for v in rows(f"{C}/eval_loo/{p}/{k}/score")[0][4:7]]
    for l in open(f"{C}/coverage_eval/scores.tsv"):
        e, t, n, pr, rc, f = l.split(); p, k = e.split("_"); S[(k, TOOL[t], p)] = [float(pr), float(rc), float(f)]
    print("set\tcaller\t" + "\t".join(CASTLE) + "\tP\tR\tF1")
    for s in SETS:
        for n in dict.fromkeys(k[1] for k in S if k[0] == s):
            v = [S.get((s, n, p)) for p in CASTLE]
            print(f"{s}\t{n}\t" + "\t".join(f"{x[2]:.4f}" if x else "NA" for x in v) + "\t" + "\t".join(f"{np.mean([x[i] for x in v if x]):.4f}" for i in range(3)))
elif mode == "rescores":
    SETS = ["with", "without", "loo", "30x", "20x", "15x", "10x", "5x", "100p", "80p", "60p", "40p", "20p"]
    print("theta0\ttheta_t\tconfig\t" + "\t".join(f"{s}\t{s}_rank" if s in SETS[:3] else s for s in SETS) + "\twith_runner_up")
    for l in open(f"{W}/combos.tsv").read().split("\n")[1:]:
        if not l or not os.path.exists(f"{W}/conf/{l.split()[2]}/rescore.tsv"): continue
        a, b, key = l.split("\t"); E = {}
        for r in open(f"{W}/conf/{key}/rescore.tsv").read().split("\n")[1:]:
            if r: x = r.split("\t"); E.setdefault(x[0], {})[x[1]] = float(x[-1])
        rank = lambda s: 1 + sum(v > E[s]["janus"] for n, v in E[s].items() if n != "janus")
        nxt = max((v, n) for n, v in E["with"].items() if n != "janus")
        print(f"{a}\t{b}\t{key}\t" + "\t".join(f"{E[s]['janus']:.4f}\t{rank(s)}" if s in SETS[:3] else f"{E[s]['janus']:.4f}" for s in SETS) + f"\t{nxt[1]} {nxt[0]:.4f}")
elif mode == "table":
    rows = []
    for l in open(f"{W}/combos.tsv").read().split("\n")[1:]:
        if not l: continue
        a, b, key = l.split("\t")
        for f in sorted(glob.glob(f"{W}/conf/{key}/score/*/*/minda/results/overall.tsv")):
            c, e = f.split("/")[-5:-3]; v = open(f).read().split("\n")[1].split("\t")
            rows.append((int(a), int(b), key, e, c, int(v[7]), float(v[4]), float(v[5]), float(v[6])))
    with open(f"{W}/scores.tsv", "w") as o:
        o.write("theta0\ttheta_t\tconfig\texperiment\tconsensus\tcalls\tprecision\trecall\tf1\n")
        for r in sorted(rows): o.write("%d\t%d\t%s\t%s\t%s\t%d\t%.4f\t%.4f\t%.4f\n" % r)
    with open(f"{W}/summary.tsv", "w") as o:
        o.write("theta0\ttheta_t\tconfig\tfull_f1_without_janus\tfull_f1_with_janus\tcoverage_f1_without_janus\tpurity_f1_without_janus\n")
        for (a, b, key) in sorted({r[:3] for r in rows}):
            kind = lambda e: "full" if "_" not in e else "coverage" if e.endswith("x") else "purity"
            m = lambda c, k: np.mean([r[8] for r in rows if r[:3] == (a, b, key) and r[4] == c and kind(r[3]) == k])
            o.write(f"{a}\t{b}\t{key}\t{m('without', 'full'):.4f}\t{m('with', 'full'):.4f}\t{m('without', 'coverage'):.4f}\t{m('without', 'purity'):.4f}\n")
    print(open(f"{W}/summary.tsv").read())
'''
SCRIPTS["subtract"] = r'''#!/usr/bin/env python3
import sys, collections, argparse

ap = argparse.ArgumentParser()
ap.add_argument('tumor'); ap.add_argument('normal')
ap.add_argument('--somatic-tag', action='store_true')
ap.add_argument('--min-len', type=int, default=50)
ap.add_argument('--dist', type=int, default=1000)
a = ap.parse_args()

def info(f):
    d = {}
    for kv in f[7].split(';'):
        k, _, v = kv.partition('=')
        d[k] = v
    return d

def svtype(f, d):
    if 'SVTYPE' in d: return d['SVTYPE']
    return 'INS' if len(f[4]) > len(f[3]) else 'DEL'

def is_sv(f, d):
    if f[4].startswith('<') or '[' in f[4] or ']' in f[4]: return True
    if 'SVLEN' in d:
        try: return abs(int(d['SVLEN'].split(',')[0])) >= a.min_len
        except ValueError: pass
    return abs(len(f[4]) - len(f[3])) >= a.min_len

normal = collections.defaultdict(list)
for l in open(a.normal):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t'); d = info(f)
    if is_sv(f, d): normal[(f[0], svtype(f, d))].append(int(f[1]))

for l in open(a.tumor):
    if l[0] == '#': sys.stdout.write(l); continue
    f = l.rstrip('\n').split('\t'); d = info(f)
    if f[6] not in ('PASS', '.') or not is_sv(f, d): continue
    if a.somatic_tag and 'SOMATIC' not in d: continue
    if any(abs(p - int(f[1])) <= a.dist for p in normal[(f[0], svtype(f, d))]): continue
    sys.stdout.write(l)
'''
SCRIPTS["savana"] = r'''import sys
import multiprocessing.process as mp
from savana.savana import main

join = mp.BaseProcess.join


def reap(self, timeout=None):
    join(self, timeout)
    if timeout is not None and self.is_alive():
        self.kill()
        join(self)


mp.BaseProcess.join = reap
sys.exit(main())
'''
SCRIPTS["simtruth"] = r'''import sys, json, re, bisect, os
G = "/home/yuio/janus/simulation/phylo/genomes"
OUT = sys.argv[1] if len(sys.argv) > 1 else "/home/yuio/janus/bench/runs/sim/truth"
plan = json.load(open(f"{G}/tumor/plan.json"))
PAR = plan["parent"]


def base(c):
    return "_".join(c.split("_")[:2])


def load(path, segs):
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        m = re.search(r"cg:Z:(\d+)M$", f[12]) if len(f) > 12 else None
        if not m: continue
        q, qs, qe, st, t, ts, te = base(f[0]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8])
        segs.setdefault(q, []).append((qs, qe, st, t, ts, te))


PAF = {}


def paf(node):
    if node not in PAF:
        segs = {}
        if node == "germline":
            for c in os.listdir(f"{G}/germline"):
                for h in ("A", "B"):
                    p = f"{G}/germline/{c}/sim.hap{h}.paf"
                    if os.path.exists(p): load(p, segs)
        else:
            load(f"{G}/tumor/{node}/sim.paf", segs)
        for k in segs: segs[k].sort()
        PAF[node] = {k: (v, [x[0] for x in v]) for k, v in segs.items()}
    return PAF[node]


def step(node, c, p):
    d = paf(node).get(base(c))
    if not d: return None
    v, starts = d
    i = bisect.bisect_right(starts, p) - 1
    if i < 0: i, p, snap = 0, v[0][0], True
    qs, qe, st, t, ts, te = v[i]
    snap = False
    if p >= qe:
        if i + 1 < len(v) and v[i + 1][0] - p < p - (qe - 1):
            qs, qe, st, t, ts, te = v[i + 1]; p = qs
        else:
            p = qe - 1
        snap = True
    return (t, ts + (p - qs) if st == "+" else te - 1 - (p - qs), snap)


def lift(node, c, p):
    snapped = False
    while True:
        r = step(node, c, p)
        if r is None: return None
        c, p, s = r
        snapped |= s
        if node == "germline": return c, p, snapped
        node = PAR[node] or "germline"


def chain(n):
    out = []
    while n:
        out.append(n); n = PAR[n]
    return out


rows = [l.rstrip("\n").split("\t") for l in open(f"{G}/tumor/mutations.tsv")]
hdr, rows = rows[0], [r for r in rows[1:] if r[5] == "SV"]
lifted, failed, snaps = {}, 0, 0
for r in rows:
    mid, birth, c, pos, end, svt, svlen = r[0], r[1], r[2], int(r[3]), int(r[4]), r[6], int(r[7])
    info = dict(x.split("=", 1) for x in r[10].split(";") if "=" in x)
    start = PAR[birth] or "germline"
    a, b = lift(start, c, pos), lift(start, c, end)
    if not a or not b: failed += 1; continue
    J = []
    if svt in ("dDUP", "rTRA"):
        tg = lift(start, info["TARGET_CHROM"], int(info["TARGET"]))
        if not tg: failed += 1; continue
        J = [(tg[0], tg[1], a[0], a[1]), (b[0], b[1], tg[0], tg[1])]
        snaps += tg[2]
    elif svt == "INS":
        J = [(a[0], a[1], a[0], a[1])]
    else:
        J = [(a[0], a[1], b[0], b[1])]
    snaps += a[2] or b[2]
    lifted[mid] = dict(birth=birth, type=svt, len=svlen, j=J)
for node in plan["nodes"]:
    p = f"{G}/tumor/{node}/sim.vcf"
    if not os.path.exists(p): continue
    start = PAR[node] or "germline"; parts = {}
    for l in open(p):
        if l[0] == "#" or "SVTYPE=rTRA" not in l: continue
        x = l.split("\t"); info = dict(kv.split("=", 1) for kv in x[7].split(";") if "=" in kv)
        parts.setdefault(info["SVID"].replace("Imported_0_", ""), []).append((x[0], int(x[1]), int(info["END"]), int(info["SVLEN"])))
    for mid, ps in parts.items():
        if mid not in lifted: continue
        L = [(lift(start, c, a), lift(start, c, e), n) for c, a, e, n in ps]
        if any(x[0] is None or x[1] is None for x in L): continue
        full = [x for x in L if x[2] > 0]; empty = [x for x in L if x[2] == 0]
        J = []
        if len(full) == 1 and empty:
            (a1, a2, _), tg = full[0], empty[0][0]
            J = [(tg[0], tg[1], a1[0], a1[1]), (a2[0], a2[1], tg[0], tg[1]), (tg[0], tg[1], tg[0], tg[1]), (a1[0], a1[1], a2[0], a2[1])]
        elif len(full) == 2:
            (a1, a2, _), (b1, b2, _) = full
            J = [(x[0], x[1], y[0], y[1]) for x in (a1, a2) for y in (b1, b2)] + [(a1[0], a1[1], a2[0], a2[1]), (b1[0], b1[1], b2[0], b2[1])]
        if J: lifted[mid]["j"] = J
for mid, d in lifted.items():
    j = d["j"]
    if d["type"] == "dDUP": d["j"] = j + [(j[0][0], j[0][1], j[0][0], j[0][1])]
    if d["type"] == "DUP": d["j"] = j + [(j[0][0], j[0][1], j[0][0], j[0][1]), (j[0][2], j[0][3], j[0][2], j[0][3])]
leaves = sorted(n for n in plan["nodes"] if n.startswith("T"))
os.makedirs(OUT, exist_ok=True)
for t in leaves:
    anc = set(chain(t))
    with open(f"{OUT}/{t}.tsv", "w") as o:
        for mid, d in lifted.items():
            if d["birth"] in anc:
                for k, j in enumerate(d["j"]):
                    o.write(f"{mid}\t{d['type']}\t{d['len']}\t{k}\t{j[0]}\t{j[1]}\t{j[2]}\t{j[3]}\n")
json.dump(dict(parent=PAR, leaves=leaves), open(f"{OUT}/tree.json", "w"))
n = len(rows)
print(f"SVs {n}, lifted {len(lifted)} ({len(lifted) / n:.3%}), failed {failed}, snapped {snaps}")
for t in leaves:
    k = len(set(l.split("\t")[0] for l in open(f"{OUT}/{t}.tsv")))
    print(t, k, end="  ")
print()
'''
SCRIPTS["simscore"] = r'''import sys, re, bisect
TOL, MIN = 500, 50
ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def truth(path):
    sv, J = {}, []
    for l in open(path):
        mid, typ, ln, k, c1, p1, c2, p2 = l.rstrip("\n").split("\t")
        if typ not in ("dDUP", "rTRA") and int(ln) < MIN: continue
        sv[mid] = typ
        J.append((c1, int(p1), c2, int(p2), mid))
    return sv, J


def calls(path, passonly=True):
    out = []
    for l in open(path):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t")
        if passonly and f[6] not in ("PASS", "."): continue
        info = dict(x.split("=", 1) if "=" in x else (x, "") for x in f[7].split(";"))
        t = info.get("SVTYPE", "")
        m = ALT.search(f[4])
        if m or t in ("BND", "TRA"):
            if m: out.append([(f[0], int(f[1]), m.group(1), int(m.group(2)))])
            elif "CHR2" in info: out.append([(f[0], int(f[1]), info["CHR2"], int(info.get("END", info.get("POS2", f[1]))))])
            continue
        ln = abs(int(info.get("SVLEN", "0") or 0))
        end = int(info.get("END", f[1]))
        if t == "INS":
            if ln and ln < MIN: continue
            out.append([(f[0], int(f[1]), f[0], int(f[1]))])
        else:
            if max(ln, end - int(f[1])) < MIN: continue
            out.append([(f[0], int(f[1]), f[0], end)])
    return out


def index(J):
    by = {}
    for j in J:
        by.setdefault(j[0], []).append((j[1], j))
        by.setdefault(j[2], []).append((j[3], j))
    return {c: (sorted(v), [x[0] for x in sorted(v)]) for c, v in by.items()}


def near(a, b):
    return (a[0] == b[0] and abs(a[1] - b[1]) <= TOL and a[2] == b[2] and abs(a[3] - b[3]) <= TOL) or \
           (a[0] == b[2] and abs(a[1] - b[3]) <= TOL and a[2] == b[0] and abs(a[3] - b[1]) <= TOL)


def match(c, I):
    hit = set()
    for (c1, p1, c2, p2) in c:
        d = I.get(c1)
        if not d: continue
        v, ps = d
        for i in range(bisect.bisect_left(ps, p1 - TOL), bisect.bisect_right(ps, p1 + TOL)):
            j = v[i][1]
            if near((c1, p1, c2, p2), j): hit.add(j[4])
    return hit


def score(tpath, cpath, passonly=True):
    sv, J = truth(tpath)
    I = index(J)
    C = calls(cpath, passonly)
    found, tp = set(), 0
    for c in C:
        h = match(c, I)
        if h: tp += 1; found |= h
    P = tp / len(C) if C else 0
    R = len(found) / len(sv) if sv else 0
    F = 2 * P * R / (P + R) if P + R else 0
    return dict(truth=len(sv), calls=len(C), tp_calls=tp, found=len(found), P=P, R=R, F1=F, found_ids=found)


if __name__ == "__main__":
    r = score(sys.argv[1], sys.argv[2], len(sys.argv) < 4)
    print(f"truth {r['truth']} calls {r['calls']} recall {r['R']:.3f} precision {r['P']:.3f} f1 {r['F1']:.3f}")
'''
SCRIPTS["norm"] = r'''import sys, re

ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
name = sys.argv[2]
print("##fileformat=VCFv4.2")
print(f"##source={name}")
print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">')
print('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">')
print('##INFO=<ID=END,Number=1,Type=Integer,Description="end">')
print('##INFO=<ID=MATEID,Number=1,Type=String,Description="mate">')
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
seen, n = set(), 0
for l in open(sys.argv[1]):
    if l[0] == "#":
        continue
    f = l.rstrip("\n").split("\t")
    if f[6] not in ("PASS", "."):
        continue
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    t, c, p = info.get("SVTYPE", "BND"), f[0], int(f[1])
    m = ALT.search(f[4])
    n += 1
    if m or info.get("CHR2", c) != c:
        c2, p2 = (m.group(1), int(m.group(2))) if m else (info["CHR2"], int(info.get("CHR2_POS", info.get("END", p))))
        k = tuple(sorted([(c, p), (c2, p2)]))
        if k in seen:
            continue
        seen.add(k)
        print(f"{c}\t{p}\t{name}_{n}_1\tN\tN[{c2}:{p2}[\t.\tPASS\tSVTYPE=BND;MATEID={name}_{n}_2")
        print(f"{c2}\t{p2}\t{name}_{n}_2\tN\tN[{c}:{p}[\t.\tPASS\tSVTYPE=BND;MATEID={name}_{n}_1")
        continue
    e = int(info.get("END", p))
    try:
        ln = abs(int(info["SVLEN"].split(",")[0]))
    except (KeyError, ValueError):
        ln = abs(e - p)
    if t == "INS":
        print(f"{c}\t{p}\t{name}_{n}\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN={ln};END={p}")
    else:
        t = t if t in ("DEL", "DUP", "INV") else "BND"
        if t == "BND":
            print(f"{c}\t{p}\t{name}_{n}_1\tN\tN[{c}:{e}[\t.\tPASS\tSVTYPE=BND;MATEID={name}_{n}_2")
            print(f"{c}\t{e}\t{name}_{n}_2\tN\tN[{c}:{p}[\t.\tPASS\tSVTYPE=BND;MATEID={name}_{n}_1")
        else:
            print(f"{c}\t{p}\t{name}_{n}\tN\t<{t}>\t.\tPASS\tSVTYPE={t};SVLEN={ln};END={e}")
'''
SCRIPTS["thin"] = r'''import mmap, struct, sys
src, out, stride = sys.argv[1], sys.argv[2], int(sys.argv[3])
n = total = longest = i = 0; off = []
with open(out, "wb") as o, open(src, "rb") as f:
    o.write(b"\0" * 32); pos = 32
    m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ); p, end = 32, len(m)
    while p < end:
        tag = struct.unpack_from("<I", m, p)[0]; ln = tag & 0x7fffffff
        sz = 4 + (ln if tag >> 31 else (ln + 3) // 4)
        if i % stride == 0:
            o.write(m[p:p + sz]); off.append(pos); pos += sz; total += ln; n += 1; longest = max(longest, ln)
        i += 1; p += sz
    off.append(pos)
    o.seek(0); o.write(b"SVB1" + struct.pack("<IQQ", 1, n, total) + b"\0" * 8)
with open(out + ".idx", "wb") as o:
    o.write(b"SVBI" + struct.pack("<IQQ", 1, n, longest) + b"\0" * 8)
    o.write(struct.pack("<%dQ" % len(off), *off))
print(out, n, "reads", total, "bases", "%.2fx" % (total / 3117292070))
'''
SCRIPTS["simjoint"] = r'''import sys, io, os, json, pickle, collections, numpy as np
H = "/home/yuio/janus/bench"
G = {"__name__": "x", "SCRIPTS": SCRIPTS, "run_score": run_score}; sys.argv = sys.argv[:1]
exec(G["SCRIPTS"]["verdict"].split("t0 = time.time()")[0], G)
ALT, np_ = G["ALT"], np
M = pickle.load(open(f"{H}/verdict/model.pkl", "rb")); ms = M["models"]
ns = {}; exec(G["SCRIPTS"]["simscore"].split("if __name__")[0], ns)
MEM = {}; ns["open"] = lambda p, *x: io.StringIO(MEM[p]) if p in MEM else open(p, *x)
TR = f"{H}/runs/sim/truth"; P = os.environ.get("SIMPAIR", f"{H}/runs/sim/settings"); JD = os.environ.get("SIMJOINT", ""); T = ["T%02d" % i for i in range(1, 21)]
jd = lambda s: f"{JD}/{s}" if JD else f"{P}/{s}/joint"
SETS = [("20-20", "20_20"), ("20-10", "20_10"), ("20-5", "20_5"), ("10-10", "10_10"), ("5-5", "5_5")]

sys.path.insert(0, "/home/yuio/janus/src")
from janus import features as JF, classify as JC, joint as JJ

def load(f):
    L, J = [], []
    for l in open(f):
        if l[0] == "#": continue
        g = l.rstrip("\n").split("\t")
        MEM["m"] = "\t".join(g[:6] + ["PASS"] + g[7:]) + "\n"; c = ns["calls"]("m")
        L.append(g); J.append(c[0] if c else None)
    return L, J

out = {}
for name, s in SETS:
    if not os.path.exists(f"{jd(s)}/calls.vcf") or len([t for t in T if os.path.exists(f"{P}/{s}/{t}/t.cnt.hist")]) < 20: print(name, "not ready"); continue
    TI, NT = [], []
    for t in T:
        sv, TJ = ns["truth"](f"{TR}/{t}.tsv"); TI.append(ns["index"](TJ)); NT.append(len(sv))
    pair = []
    for k, t in enumerate(T):
        d = f"{P}/{s}/{t}"; L, J = load(f"{d}/calls.vcf"); _, _, X, K, _, _, _ = JF.records(f"{d}/calls.vcf")
        pr, du, cut = JC.classify(M, X, K, f"{d}/t.cnt.hist"); mt = JF.minimum(X); pp = np.where((pr >= cut) & ~du)[0]
        hits = [ns["match"](J[i], TI[k]) if J[i] else set() for i in pp]; f = set().union(*hits) if hits else set()
        pair.append((sum(1 for h in hits if h) / max(len(pp), 1), len(f) / NT[k], len(pp), mt))
    L, J = load(f"{jd(s)}/calls.vcf"); _, _, X, K, own, SM, SP = JF.records(f"{jd(s)}/calls.vcf"); N = len(L)
    hists = [f"{jd(s)}/{t}.cnt.hist" for t in T]; mins = [int(open(f"{jd(s)}/{t}.minimum").read()) for t in T]
    pr, du, cut, acc1, L1, D1 = JJ.run(M, X, K, own, SM, SP, hists, mins, clade=False)
    NB, grp = JJ.neighbors(K, du)
    def evaluate(lists, adds):
        A = []
        for k in range(20):
            hits = [ns["match"](J[i], TI[k]) if J[i] else set() for i in lists[k]]; f = set().union(*hits) if hits else set()
            A.append((sum(1 for h in hits if h) / max(len(hits), 1), len(f) / NT[k], len(hits), len(adds[k]), sum(1 for h in hits[len(hits) - len(adds[k]):] if h)))
        return A
    clade = lambda lists, mu: JJ.clade_pass(acc1, lists, pr, cut, du, own, grp, 20, mu, M["cut"])
    A1 = evaluate(L1, D1); L2, D2 = clade(L1, JJ.MU); A = evaluate(L2, D2)
    pickle.dump(dict(pr=pr, du=du, own=own, NB=NB, J=J, SM=SM, SP=SP, cut=cut, pair=[(x[0], x[1]) for x in pair]), open(f"{jd(s)}/carry.pkl", "wb"))
    res = lambda A_: dict(P=[x[0] for x in A_], R=[x[1] for x in A_], calls=[x[2] for x in A_], carried=[x[3] for x in A_], carried_true=[x[4] for x in A_])
    out[name] = dict(truth=NT, pair=dict(P=[x[0] for x in pair], R=[x[1] for x in pair], calls=[x[2] for x in pair], minimum=[int(x[3]) for x in pair]), first=res(A1), joint=res(A))
    med = lambda A_: [float(np.median([x[0] for x in A_])), float(np.median([x[1] for x in A_])), float(np.median([2 * x[0] * x[1] / max(x[0] + x[1], 1e-9) for x in A_]))]
    out[name]["mu"] = dict(mu=[0.3, 0.5, 1.0, 1.5, 2.0], PRF=[med(evaluate(*clade(L1, m_))) for m_ in (0.3, 0.5, 1.0, 1.5, 2.0)])
    if name in ("20-20", "20-5"):
        sw = []
        for r in (0.0, 0.1, 0.2, 0.3):
            sw.append([med(evaluate(*JJ.accept(pr, cut, du, own, NB, SM, SP, 20, r, fl)[1:])) for fl in (0.05, 0.1, 0.2, 0.3)])
        out[name]["sweep"] = dict(raise_=[0.0, 0.1, 0.2, 0.3], floor=[0.05, 0.1, 0.2, 0.3], PRF=sw)
    P1, R1, P2, R2 = [np.array(x) for x in (out[name]["pair"]["P"], out[name]["pair"]["R"], out[name]["joint"]["P"], out[name]["joint"]["R"])]
    F1, F2 = 2 * P1 * R1 / np.maximum(P1 + R1, 1e-9), 2 * P2 * R2 / np.maximum(P2 + R2, 1e-9)
    print(f"{name}  one at a time P {np.median(P1):.3f} R {np.median(R1):.3f} F1 {np.median(F1):.3f}  first pass P {med(A1)[0]:.3f} R {med(A1)[1]:.3f} F1 {med(A1)[2]:.3f}"
          f"  joint P {np.median(P2):.3f} R {np.median(R2):.3f} F1 {np.median(F2):.3f}  both better {int(((P2 > P1) & (R2 > R1)).sum())}/20  F1 better {int((F2 > F1).sum())}/20"
          f"  dF1 {np.min(F2 - F1):+.3f}..{np.max(F2 - F1):+.3f}  carriers {sum(x[3] for x in A1)} ({sum(x[4] for x in A1)} true) -> {sum(x[3] for x in A)} ({sum(x[4] for x in A)} true)", flush=True)
json.dump(out, open(f"{H}/runs/sim/scores.json", "w"), indent=1)
'''
SCRIPTS["supp"] = r'''import os, sys, json, itertools, tempfile
os.environ.setdefault("OUTV", "/home/yuio/janus/bench/supp/work")
job = sys.argv[1] if len(sys.argv) > 1 else "all"
G = {"__name__": "supp", "SCRIPTS": SCRIPTS, "run_score": run_score}; sys.argv = ["supp"]
v = G["SCRIPTS"]["verdict"]
exec(v.split("t0 = time.time()")[0], G)
exec("t0 = time.time()" + v.split("t0 = time.time()")[1].split("D = {(p, k)")[0], G)
exec("D = {(p, k): feats(f, p) for p in FULL for k, f in sets(p).items()}", G)
np = G["np"]; D = G["D"]; TRAIN = G["TRAIN"]; OUT = os.environ.get("SUPPOUT", "/home/yuio/janus/bench/supp"); os.makedirs(OUT, exist_ok=True)
from sklearn.metrics import roc_auc_score, average_precision_score

def outcome(ms, cuts, c0, p):
    r = {}
    for (q, k), (L, X, K, y, keep) in D.items():
        if q != p: continue
        pr = G["prob"](ms, X); m = int(round(X[0, 56])); cut = cuts.get(m, c0); d = G["mark"](K, pr)
        with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as o:
            for l, x, u in zip(L, pr, d):
                if x >= cut and not u: g = l.split("\t"); g[6] = "PASS"; o.write("\t".join(g))
        n, tpc, tpt, P, R, F = G["score"](o.name, p); os.unlink(o.name)
        ok = ~d
        auc = roc_auc_score(y[ok], pr[ok]) if 0 < y[ok].sum() < ok.sum() else float("nan")
        ap = average_precision_score(y[ok], pr[ok]) if y[ok].sum() > 0 else float("nan")
        r[k] = dict(P=P, R=R, F1=F, roc_auc=auc, pr_auc=ap, cut=cut, minimum=m)
    return r

def evaluate(name, splits):
    rows = []
    for tr, te in splits:
        ms, cuts, c0, ts, nt = G["fit"](tuple(tr))
        for p in te:
            for k, x in outcome(ms, cuts, c0, p).items(): rows.append(dict(variant=name, train=",".join(tr), pair=p, setting=k, **x))
    return rows

def write(rows, f):
    keys = list(rows[0])
    with open(f, "w") as o:
        o.write("\t".join(keys) + "\n")
        for r in rows: o.write("\t".join(str(round(r[k], 4)) if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")

LOPO = [([q for q in TRAIN if q != p], [p]) for p in TRAIN]
if job in ("auc", "all"):
    write(evaluate("D", LOPO), f"{OUT}/heldout_auc.tsv")
if job in ("splits", "all"):
    sp = [(list(c), [p for p in TRAIN if p not in c]) for c in itertools.combinations(TRAIN, 3)]
    write(evaluate("D", sp), f"{OUT}/splits_3_2.tsv")
if job in ("ablation", "all"):
    base = dict(G["HGB"]); rows = []
    grid = [("max_iter", [100, 600]), ("learning_rate", [0.02, 0.1]), ("max_leaf_nodes", [7, 31]), ("min_samples_leaf", [10, 160]), ("l2_regularization", [0.0, 5.0])]
    for key, vals in grid:
        for val in vals:
            G["HGB"] = dict(base, **{key: val}); rows += evaluate(f"{key}={val}", LOPO)
    G["HGB"] = base
    f5 = G["fitm"]
    G["fitm"] = lambda X, y: [G["HistGradientBoostingClassifier"](**G["HGB"], random_state=0).fit(X, y)]
    rows += evaluate("seeds=1", LOPO)
    G["fitm"] = f5
    write(rows, f"{OUT}/ablation.tsv")
if job in ("grid", "xgb"):
    NT = int(os.environ.get("OMP_NUM_THREADS", "4")); i0, n0 = map(int, os.environ.get("SHARD", "0/1").split("/"))
    if job == "grid":
        AX = [("learning_rate", [0.02, 0.05, 0.1]), ("max_iter", [150, 300, 600]), ("max_leaf_nodes", [7, 15, 31]), ("min_samples_leaf", [10, 40, 160]), ("l2_regularization", [0.0, 1.0, 5.0])]
        fitm = lambda X, y, h: [G["HistGradientBoostingClassifier"](**dict(G["HGB"], **h), random_state=s).fit(X, y) for s in range(5)]
    else:
        import xgboost as xgb
        from sklearn.model_selection import train_test_split
        AX = [("learning_rate", [0.02, 0.05, 0.1]), ("n_estimators", [150, 300, 600]), ("max_depth", [3, 4, 6]), ("min_child_weight", [1, 5, 20]), ("reg_lambda", [0.0, 1.0, 5.0])]
        def fitm(X, y, h):
            ms = []
            for s in range(5):
                a, b, ya, yb = train_test_split(X, y, test_size=0.1, random_state=s, stratify=y)
                ms.append(xgb.XGBClassifier(**h, subsample=0.8, colsample_bytree=0.8, tree_method="hist", early_stopping_rounds=10, n_jobs=NT, random_state=s, verbosity=0).fit(a, ya, eval_set=[(b, yb)], verbose=False))
            return ms
    import glob, collections
    f = f"{OUT}/{job}.{i0}.tsv"; seen = collections.Counter()
    for g in glob.glob(f"{OUT}/{job}.*.tsv"): seen.update(l.split("\t")[0] for l in open(g).read().split("\n")[1:] if l)
    done = {v for v, n in seen.items() if n >= len(D)}
    out = open(f, "a")
    if not out.tell(): out.write("variant\tpair\tsetting\tcalls\tP\tR\tF1\troc_auc\tpr_auc\tcut\tminimum\tinner_f1\n")
    tmp = f"{OUT}/{job}.{i0}.pass.vcf"
    def passed(ms, cuts, c0, p):
        r = {}
        for (q, k), (L, X, K, y, keep) in D.items():
            if q != p: continue
            pr = G["prob"](ms, X); m = int(round(X[0, 56])); cut = cuts.get(m, c0); d = G["mark"](K, pr)
            with open(tmp, "w") as o:
                for l, x, u in zip(L, pr, d):
                    if x >= cut and not u: g = l.split("\t"); g[6] = "PASS"; o.write("\t".join(g))
            n, tpc, tpt, P, R, F = G["score"](tmp, p); ok = ~d
            auc = roc_auc_score(y[ok], pr[ok]) if 0 < y[ok].sum() < ok.sum() else float("nan")
            ap = average_precision_score(y[ok], pr[ok]) if y[ok].sum() > 0 else float("nan")
            r[k] = dict(calls=n, P=P, R=R, F1=F, roc_auc=auc, pr_auc=ap, cut=cut, minimum=m)
        return r
    runs = lambda tr: [k for k in D if k[0] in tr and k[1] != "100p"]
    XY = lambda ks: (np.vstack([D[k][1][D[k][4]] for k in ks]), np.concatenate([D[k][3][D[k][4]] for k in ks]))
    todo = [c for c in itertools.product(*[v for _, v in AX]) if ",".join(f"{a}={v}" for (a, _), v in zip(AX, c)) not in done]
    for c in todo[i0::n0]:
        h = dict(zip([a for a, _ in AX], c)); name = ",".join(f"{a}={v}" for a, v in h.items())
        M = {tr: fitm(*XY(runs(tr)), h) for r in (4, 3) for tr in itertools.combinations(TRAIN, r)}
        for p in TRAIN:
            tr = tuple(q for q in TRAIN if q != p); oof, ys, gs = [], [], []
            for q in tr:
                mi = M[tuple(x for x in tr if x != q)]
                for k in runs((q,)):
                    Xk, yk = D[k][1][D[k][4]], D[k][3][D[k][4]]; oof.append(G["prob"](mi, Xk)); ys.append(yk); gs.append(np.round(Xk[:, 56]))
            o, y, g = np.concatenate(oof), np.concatenate(ys), np.concatenate(gs)
            cuts = {int(m): float(max(G["CUTS"], key=lambda t: G["f1"](o[g == m], y[g == m], t))) for m in set(g)}
            c0 = float(max(G["CUTS"], key=lambda t: G["f1"](o, y, t)))
            inner = float(np.mean([G["f1"](a, b, cuts.get(int(m[0]), c0)) for a, b, m in zip(oof, ys, gs)]))
            for k, x in passed(M[tr], cuts, c0, p).items():
                out.write("\t".join([name, p, k, str(x["calls"])] + [f"{x[v]:.4f}" for v in ("P", "R", "F1", "roc_auc", "pr_auc", "cut")] + [str(x["minimum"]), f"{inner:.4f}"]) + "\n")
        out.flush(); print(job, name, flush=True)
if job == "firstpass":
    orig = dict(D); rows = []
    for name in ("excluded", "kept"):
        O = f"{OUT}/firstpass/{name}"; os.makedirs(O, exist_ok=True)
        for k, v in orig.items(): D[k] = v if name == "excluded" else v[:4] + (np.ones(len(v[3]), bool),)
        for p in TRAIN:
            tr = tuple(q for q in TRAIN if q != p); ms, cuts, c0, ts, nt = G["fit"](tr)
            for k, x in outcome(ms, cuts, c0, p).items(): rows.append(dict(variant=name, pair=p, setting=k, train_records=nt, **x))
            L, X, K, y, keep = D[(p, "full")]; pr = G["prob"](ms, X); cut = cuts.get(int(round(X[0, 56])), c0); d = G["mark"](K, pr)
            with open(f"{O}/{p}.vcf", "w") as o:
                o.write("".join(l for l in open(G["sets"](p)["full"]) if l[0] == "#"))
                for l, r, u in zip(L, pr, d):
                    g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if r >= cut else "LOW"; g[7] += ";VP=%.4f" % r; o.write("\t".join(g))
    D.update(orig)
    write(rows, f"{OUT}/firstpass.tsv")
if job == "features":
    GROUPS = [("Tumor reads, single anchor pair", [1, 2, 3, 4, 5, 35]), ("Normal reads, single anchor pair", [6, 7, 8, 9, 10, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34]),
              ("Tumor reads, anchor grid and paths", [39, 40, 49]), ("Normal reads, anchor grid and paths", [38, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52]),
              ("Graph context", [12, 13, 36, 37, 53, 54, 55]), ("Alignment quality", [14, 15]), ("Reference context and depth", [16, 17, 23]),
              ("Type and size", [18, 19, 20, 21, 22]), ("Minimum and count-stage context", [56, 57, 58, 59, 60])]
    F = json.load(open("/home/yuio/janus/bench/verdict/folds.json")); out = {}
    for p in TRAIN + ["hg008", "colo829"]:
        tr = [q for q in TRAIN if q != p]; ks = [k for k in D if k[0] in tr and k[1] != "100p"]
        ms = G["fitm"](np.vstack([D[k][1][D[k][4]] for k in ks]), np.concatenate([D[k][3][D[k][4]] for k in ks]))
        L, X, K, y, keep = D[(p, "full")]; cut = F[p]["cut"]
        def f1of(Z):
            pr = G["prob"](ms, Z); d = G["mark"](K, pr)
            return float(G["f1"](pr[~d], y[~d], cut)), pr
        base, pr = f1of(X)
        assert np.abs(pr - np.array(F[p]["prob"])).max() < 1e-3, p
        out[p] = dict(base=base, cut=cut, drop={})
        for name, ix in GROUPS:
            v = []
            for r in range(5):
                Z = X.copy(); Z[:, ix] = X[np.random.RandomState(r).permutation(len(X))][:, ix]
                v.append(base - f1of(Z)[0])
            out[p]["drop"][name] = v
        print(p, round(base, 4), {n: round(float(np.mean(v)), 4) for n, v in out[p]["drop"].items()}, flush=True)
    json.dump(dict(groups=GROUPS, pairs=out), open(f"{OUT}/features.json", "w"), indent=1)
print("done", job)
'''
SCRIPTS["labels"] = r'''import os, sys
os.environ.setdefault("OUTV", "/home/yuio/janus/bench/supp/labels11")
os.makedirs(os.environ["OUTV"], exist_ok=True)
G = {"__name__": "supp", "SCRIPTS": SCRIPTS, "run_score": run_score}; sys.argv = ["supp"]
v = G["SCRIPTS"]["verdict"]
head, main = v.split("t0 = time.time()")
head = head.replace('f"{H}/eval_sets/{p}/all_without/None_minda_ensemble.vcf"', 'f"{TR}/castle_without_snape/{p}/None_minda_ensemble.vcf"')
assert "castle_without_snape" in head
exec(head, G)
exec("t0 = time.time()" + main, G)
'''
SCRIPTS["select"] = r'''import sys, re, os
import numpy as np

G = 3100000000
M = (1 << 64) - 1


def load(p, s):
    names, lens = [], []
    for l in open(f"{O}/names/{p}.{s}.tsv"):
        a, b = l.rstrip("\n").split("\t")
        names.append(a)
        lens.append(int(b))
    n = len(names)
    order = list(range(n))
    seed = 0x9e3779b97f4a7c15
    for i in range(n, 1, -1):
        seed ^= (seed << 13) & M
        seed ^= seed >> 7
        seed ^= (seed << 17) & M
        j = seed % i
        order[i - 1], order[j] = order[j], order[i - 1]
    order = np.array(order)
    ln = np.array(lens, dtype=np.int64)[order]
    return names, order, ln, np.cumsum(ln)


def take(cum, ln, skip, want):
    before = cum - ln
    a = int(np.searchsorted(before, skip, side="left"))
    base = int(before[a]) if a < len(cum) else int(cum[-1])
    b = int(np.searchsorted(cum, base + want, side="left"))
    return a, b + 1, int(cum[b] - base)


p, O = sys.argv[1], sys.argv[2]
os.makedirs(f"{O}/lists", exist_ok=True)
D = {s: load(p, s) for s in "nt"}


def put(path, parts):
    with open(path, "w") as o:
        for s, a, b in parts:
            names, order, ln, cum = D[s]
            o.write("\n".join(names[i] for i in order[a:b]) + "\n")


for x in (5, 10, 15, 20, 30):
    for s in "nt":
        names, order, ln, cum = D[s]
        if cum[-1] < x * G:
            continue
        a, b, tot = take(cum, ln, 0, x * G)
        put(f"{O}/lists/{p}.{s}.{x}x.txt", [(s, a, b)])
if p in ("hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"):
    nb = int(D["n"][3][-1])
    nm = min(nb - 24 * G - G // 10, 30 * G)
    a, b, tot = take(D["n"][3], D["n"][2], 0, nm)
    put(f"{O}/lists/{p}.n.purity.txt", [("n", a, b)])
    for q in (20, 40, 60, 80, 100):
        ta, tb, tt = take(D["t"][3], D["t"][2], 0, q * 30 * G // 100)
        parts, n, tot = [("t", ta, tb)], tb - ta, tt
        if q < 100:
            na, nb2, nt = take(D["n"][3], D["n"][2], nm, (100 - q) * 30 * G // 100)
            parts.append(("n", na, nb2))
            n += nb2 - na
            tot += nt
        put(f"{O}/lists/{p}.t.{q}p.tumor_part.txt", parts[:1])
        if q < 100:
            put(f"{O}/lists/{p}.t.{q}p.normal_part.txt", parts[1:])
print(p, "read lists written")
'''

def run_score(t, c):
    out = io.StringIO(); a = sys.argv
    sys.argv = ["score", t, c]
    with contextlib.redirect_stdout(out): exec(SCRIPTS["score"], {"__name__": "__main__"})
    sys.argv = a
    return out.getvalue()
if __name__ == "__main__":
    cmd = sys.argv[1]; sys.argv = [cmd] + sys.argv[2:]
    exec(SCRIPTS[cmd], {"__name__": "__main__", "SCRIPTS": SCRIPTS, "run_score": run_score})
