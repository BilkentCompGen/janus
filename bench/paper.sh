#!/usr/bin/env bash
# usage: bench/paper.sh <experiment> [pairs or tools]
#   ref          LCP core index and minimap2 indexes of the references (hg38 with chrY PAR masked, CHM13 v2), SVDSS index
#   fetch        download the read sets of the paper's pairs into hifi/<pair>
#   align        minimap2 of both read sets of each CASTLE pair, measured (runs/align)
#   full         snape on each pair, measured (runs/full); packing is from the read sets given in hifi/<pair>
#   verdict      held-out classifier over the full runs (verdict/)
#   truth        CASTLE consensus without snape and without Severus from the published calls (hifi/truth)
#   eval         Minda ensemble with snape as a twelfth member (eval/<pair>/main, runs/minda.txt)
#   coverage     snape at 5-30x on the CASTLE pairs, with the classifier trained at the same depth
#   purity       snape on tumor/normal mixtures of HCC1954 and HCC1395 at 30x, with the 30x classifier
#   phase        Clair3 + LongPhase on each normal, haplotag both BAMs, measured (runs/phase)
#   competitors  every other caller on the CASTLE pairs, measured (runs/competitors/<tool>/<pair>)
#   table        F1 of every caller against the consensus without snape and without Severus (runs/competitors_table.tsv)
#   sim          the 20-tumor simulation: truth lift, snape on each tumor against the normal, classifier, scores (runs/sim)
#   figures      the paper's figures and the Overleaf folder (paper/overleaf)
#   guard        memory guard, disk guard and load log (run as services)
#   caller       one caller on one pair (used by competitors): caller <tool> <tumor.bam> <normal.bam> <out dir> <threads>
S=/home/yuio/snape; H=$S/bench; R=$H/runs; K=$H/callers; L=$S/lib/callers; G=3100000000; V=$S/.venv-ml/bin/python
CASTLE="hcc1954 h2009 hcc1937 h1437 hcc1395"; LABELED="$CASTLE hg008 colo829"
TOOLS="severus severus_hp savana nanomonsv sniffles2 sniffles2_mosaic svision-pro cutesv debreak dysgu kled svdss delly"
X=${1:?experiment, see the header of this script}; shift
Q=$R/queue.log; mkdir -p $R
log() { echo "$(date '+%F %T') [$X] $*" >> $Q; }
py() { python3 -c "$PY" "$@"; }
pyml() { $V -c "$PY" "$@"; }
read -r -d '' PY <<'PYEOF'
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

S = "/home/yuio/snape"
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
        t, c, p = info.get("SVTYPE"), f[0], int(f[1])
        m = ALT.search(f[4])
        if m:
            if info.get("MATEID", info.get("MATE_ID")) in mates: continue
            mates.add(f[2])
            out.append((t, c, p, m.group(1), int(m.group(2))))
        else:
            out.append((t, c, p, c, int(info.get("END", p))))
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
H = "/home/yuio/snape/bench"; TR = "/home/yuio/snape/hifi/truth"; S0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
TRUTH = {p: f"{TR}/castle_without_snape/{p}/None_minda_ensemble.vcf" for p in ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]}
TRUTH.update(hg008=f"{TR}/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz",
             colo829=f"{TR}/colo829/truth.vcf")
RUNS = __import__("os").environ.get("RUNS", f"{H}/runs/full"); OUTV = __import__("os").environ.get("OUTV", f"{H}/verdict")
TRAIN = [p for p in ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395", "hg008", "colo829"] if __import__("os").path.exists(f"{RUNS}/{p}/calls.vcf")]
EVAL = []
UNLABELED = [p for p in ["hg009", "h209", "h2126", "hs578t"] if __import__("os").path.exists(f"{RUNS}/{p}/calls.vcf")]
LABELED = TRAIN + EVAL
HGB = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=40, l2_regularization=1.0)
CUTS = np.arange(0.2, 0.81, 0.05)

def records(p):
    L, X, K, R = [], [], [], []
    for l in open(f"{RUNS}/{p}/calls.vcf"):
        if l[0] == "#": continue
        f = l.split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
        K.append((f[0], m.group(1) if m else f[0], int(f[1]), int(m.group(2)) if m else int(f[1]) + int(info["SVLEN"]), info["SVTYPE"], info["O"]))
        R.append((info["SVTYPE"], f[0], int(f[1]), m.group(1), int(m.group(2))) if m else (info["SVTYPE"], f[0], int(f[1]), f[0], int(info["END"])))
        X.append([float(v) for v in info["F"].split(",")]); L.append(l)
    y = None
    if p in TRUTH:
        by = {}
        for t in recs(TRUTH[p]): by.setdefault(t[1], []).append(t); by.setdefault(t[3], []).append(t)
        y = np.array([any(hit(r, t) for t in by.get(r[1], [])) for r in R], float)
    keep = np.array(["\tDUPL\t" not in l for l in L])
    return L, np.array(X), K, y, keep

def mark(K, pr):
    ix = sorted(range(len(K)), key=lambda i: (K[i][:4], i)); d = np.zeros(len(K), bool)
    for a in range(len(ix)):
        x = K[ix[a]]
        for b in ix[a + 1:]:
            y = K[b]
            if y[:2] != x[:2] or y[2] > x[2] + 500: break
            if abs(y[3] - x[3]) <= 500 and (x[4] != y[4] or x[5] == y[5]): d[ix[a] if pr[ix[a]] < pr[b] else b] = True
    return d

def f1(p, y, c):
    return 2 * ((p >= c) & (y > 0)).sum() / max((p >= c).sum() + (y > 0).sum(), 1)

def fit(train):
    X = np.vstack([D[q][1][D[q][4]] for q in train]); y = np.concatenate([D[q][3][D[q][4]] for q in train])
    t = time.time()
    ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(S0, S0 + 5)]
    ts = time.time() - t
    pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)
    cut = max(CUTS, key=lambda c: f1(pt, y, c))
    return ms, cut, ts, len(y)

t0 = time.time()
D = {p: records(p) for p in LABELED + UNLABELED}
log = open(f"{OUTV}/log.txt", "w"); folds = {}
for p in LABELED + UNLABELED:
    tr = [q for q in TRAIN if q != p]
    ms, cut, ts, nt = fit(tr)
    L, X, K = D[p][0], D[p][1], D[p][2]
    t = time.time(); pr = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0); ti = time.time() - t; d = mark(K, pr)
    folds[p] = dict(labeled=p in TRAIN, train=tr,
                    cut=round(float(cut), 2), train_s=round(ts, 3), train_records=nt, infer_s=round(ti, 4), records=len(L),
                    prob=[round(float(v), 4) for v in pr], label=None if D[p][3] is None else [int(v) for v in D[p][3]],
                    dupl=[int(v) for v in d], features=len(X[0]))
    with open(f"{OUTV}/{p}.vcf", "w") as o:
        o.write("".join(l for l in open(f"{RUNS}/{p}/calls.vcf") if l[0] == "#"))
        for l, q, u in zip(L, pr, d):
            f = l.split("\t"); f[6] = "DUPL" if u else "PASS" if q >= cut else "LOW"; f[7] += ";VP=%.4f" % q; o.write("\t".join(f))
    s = f"{p} {'held-out' if p in TRAIN else 'no truth'} PASS {int(((pr >= cut) & ~d).sum())} cut {cut:.2f} DUPL {int(d.sum())} records {len(L)}"
    print(s); log.write(s + "\n")
s = f"verdict ({len(TRAIN)} labeled pairs held out, {len(UNLABELED)} without truth, 5 seeds) {time.time() - t0:.1f} s"
print(s); log.write(s + "\n")
import json; json.dump(folds, open(f"{OUTV}/folds.json", "w"))
'''
SCRIPTS["covverdict"] = r'''# the verdict on the coverage runs, trained at the same depth: for every depth N and labeled pair p, the five-seed trees
# trained on the N x records of the other labeled pairs (labels from their truth sets, as in verdict.py), the cut chosen
# on those training records for the depth's F-beta (params.json verdict_f_beta), applied to p's N x calls; PASS records
# scored against p's truth with score.py's rule. Writes verdict/coverage_matched/<pair>_<N>x.vcf and scores.tsv
import sys, os, json, numpy as np
sys.argv = sys.argv[:1]
exec(SCRIPTS["verdict"].split("t0 = time.time()")[0])
C = f"{H}/runs/coverage"; O = f"{H}/verdict/coverage_matched"; os.makedirs(O, exist_ok=True)
DEPTHS = [5, 10, 15, 20, 30]
TRAIN = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]
TT = {p: recs(TRUTH[p]) for p in TRAIN}
BY = {}
for p in TRAIN:
    by = {}
    for t in TT[p]: by.setdefault(t[1], []).append(t); by.setdefault(t[3], []).append(t)
    BY[p] = by

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

def fb(p, y, c, b):
    tp = ((p >= c) & (y > 0)).sum(); n = (p >= c).sum(); t = (y > 0).sum()
    return (1 + b * b) * tp / max(b * b * t + n, 1)

def score(v, p):
    Cc = recs(v); bc = {}
    for c in Cc: bc.setdefault(c[1], []).append(c); bc.setdefault(c[3], []).append(c)
    tpc = sum(any(hit(c, t) for t in BY[p].get(c[1], [])) for c in Cc)
    tpt = sum(any(hit(t, c) for c in bc.get(t[1], []) + bc.get(t[3], [])) for t in TT[p])
    P = tpc / max(len(Cc), 1); R = tpt / max(len(TT[p]), 1)
    return len(Cc), tpc, tpt, P, R, 2 * P * R / max(P + R, 1e-9)

out = open(f"{O}/scores.tsv", "w"); out.write("pair\tdepth\tbeta\tcut\ttrain_pairs\tpass\ttruth\ttp_calls\ttp_truth\tprecision\trecall\tf1\n")
for n in DEPTHS:
    D = {}
    for p in TRAIN:
        f = f"{C}/{p}_{n}x/calls.vcf"
        if os.path.exists(f): D[p] = crecords(f, p)
    b = float(json.load(open(f"{C}/{next(iter(D))}_{n}x/monitor/params.json"))["verdict_f_beta"])
    for p in D:
        tr = [q for q in D if q != p]
        X = np.vstack([D[q][1][D[q][4]] for q in tr]); y = np.concatenate([D[q][3][D[q][4]] for q in tr])
        ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(5)]
        pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)
        cut = max(CUTS, key=lambda c: fb(pt, y, c, b))
        L, Xc, K = D[p][0], D[p][1], D[p][2]
        pr = np.mean([m.predict_proba(Xc)[:, 1] for m in ms], 0); du = mark(K, pr)
        v = f"{O}/{p}_{n}x.vcf"
        with open(v, "w") as o:
            o.write("".join(l for l in open(f"{C}/{p}_{n}x/calls.vcf") if l[0] == "#"))
            for l, q, u in zip(L, pr, du):
                g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if q >= cut else "LOW"; g[7] += ";VP=%.4f" % q; o.write("\t".join(g))
        k, tpc, tpt, P, R, F = score(v, p)
        s = f"{p}\t{n}\t{b}\t{cut:.2f}\t{','.join(tr)}\t{k}\t{len(TT[p])}\t{tpc}\t{tpt}\t{P:.4f}\t{R:.4f}\t{F:.4f}"
        print(s, flush=True); out.write(s + "\n"); out.flush()
'''
SCRIPTS["purityverdict"] = r'''import os, sys
exec(SCRIPTS["covverdict"].split("out = open(")[0])
U = f"{H}/runs/purity"; O = f"{H}/verdict/purity_matched"; os.makedirs(O, exist_ok=True)
PAIRS = ["hcc1954", "hcc1395"]; LEVELS = [20, 40, 60, 80, 100]
TR30 = {p: crecords(f"{C}/{p}_30x/calls.vcf", p) for p in TRAIN if os.path.exists(f"{C}/{p}_30x/calls.vcf")}
out = open(f"{O}/scores.tsv", "w"); out.write("pair\tpurity\tcut\ttrain_pairs\tpass\ttruth\ttp_calls\ttp_truth\tprecision\trecall\tf1\n")
for p in PAIRS:
    tr = [q for q in TR30 if q != p]
    X = np.vstack([TR30[q][1][TR30[q][4]] for q in tr]); y = np.concatenate([TR30[q][3][TR30[q][4]] for q in tr])
    ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(5)]
    pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)
    cut = max(CUTS, key=lambda c: fb(pt, y, c, 1.0))
    for q in LEVELS:
        f = f"{U}/{p}_{q}p/calls.vcf"
        if not os.path.exists(f): continue
        L, Xc, K, _, _ = crecords(f, p)
        pr = np.mean([m.predict_proba(Xc)[:, 1] for m in ms], 0); du = mark(K, pr)
        v = f"{O}/{p}_{q}p.vcf"
        with open(v, "w") as o:
            o.write("".join(l for l in open(f) if l[0] == "#"))
            for l, r, u in zip(L, pr, du):
                g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if r >= cut else "LOW"; g[7] += ";VP=%.4f" % r; o.write("\t".join(g))
        k, tpc, tpt, P, R, F = score(v, p)
        s = f"{p}\t{q}\t{cut:.2f}\t{','.join(tr)}\t{k}\t{len(TT[p])}\t{tpc}\t{tpt}\t{P:.4f}\t{R:.4f}\t{F:.4f}"
        print(s, flush=True); out.write(s + "\n"); out.flush()
'''
SCRIPTS["table"] = r'''import glob, json, os, subprocess, sys
S = "/home/yuio/snape"; R = f"{S}/bench/runs/competitors"; T = f"{S}/hifi/truth"
PAIRS = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]

def score(truth, vcf):
    o = run_score(truth, vcf).split()
    return dict(zip(o[0::2], o[1::2])) if len(o) >= 10 else {}

rows = []
for p in PAIRS:
    v = f"/tmp/snape_{p}_pass.vcf"
    with open(v, "w") as o:
        o.writelines(l for l in open(f"{S}/bench/verdict/{p}.vcf") if l[0] == "#" or "\tPASS\t" in l)
    j = json.load(open(f"{S}/bench/runs/full/{p}/run.json"))["total"]
    rows.append(("snape", p, 0, j["wall_s"], j["peak_rss_kib"] / 2**20, "clean", v))
for d in sorted(glob.glob(f"{R}/*/*/run.json")):
    j = json.load(open(d)); t = j["total"]; D = os.path.dirname(d)
    if j["params"].get("pair") not in PAIRS: continue
    rows.append((j["params"]["tool"], j["params"]["pair"], t["exit"], t["wall_s"], t["peak_rss_kib"] / 2**20, j["params"].get("run", "clean"), f"{D}/calls.vcf"))
print("tool\tpair\texit\twall_s\tpeak_gb\trun\tcalls\tP\tR\tF1\tP_ws\tR_ws\tF1_ws")
for tool, p, ex, wall, gb, kind, vcf in sorted(rows, key=lambda r: (PAIRS.index(r[1]), r[0] != "snape", r[0])):
    a = b = {}
    if ex == 0 and os.path.exists(vcf):
        a = score(f"{T}/castle_without_snape/{p}/None_minda_ensemble.vcf", vcf)
        b = score(f"{T}/castle_without_severus/{p}/None_minda_ensemble.vcf", vcf)
    print(f"{tool}\t{p}\t{ex}\t{wall:.0f}\t{gb:.1f}\t{kind}\t{a.get('calls', '-')}\t{a.get('precision', '-')}\t{a.get('recall', '-')}\t{a.get('f1', '-')}\t{b.get('precision', '-')}\t{b.get('recall', '-')}\t{b.get('f1', '-')}")

for t, n in [("sniffles2_v22", "sniffles2 2.2")]:
    for p in PAIRS:
        v = f"{R}/{t}/{p}/calls.vcf"
        if not os.path.exists(v): continue
        a = score(f"{T}/castle_without_snape/{p}/None_minda_ensemble.vcf", v); b = score(f"{T}/castle_without_severus/{p}/None_minda_ensemble.vcf", v)
        print(f"{t}\t{p}\t0\t-\t-\tnot timed\t{a.get('calls', '-')}\t{a.get('precision', '-')}\t{a.get('recall', '-')}\t{a.get('f1', '-')}\t{b.get('precision', '-')}\t{b.get('recall', '-')}\t{b.get('f1', '-')}")
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
G = "/home/samples/davp-samples/simulation/phylo/genomes"
OUT = sys.argv[1] if len(sys.argv) > 1 else "/home/yuio/snape/bench/runs/sim/truth"
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
SCRIPTS["simverdict"] = r'''import sys, os, glob, json, numpy as np
sys.argv = sys.argv[:1]
exec(SCRIPTS["covverdict"].split("out = open(")[0])
SIM = "/home/yuio/snape/bench/runs/sim/pair"
D = {p: crecords(f"{C}/{p}_20x/calls.vcf", p) for p in TRAIN}
b = float(json.load(open(f"{C}/hcc1954_20x/monitor/params.json"))["verdict_f_beta"])
X = np.vstack([D[q][1][D[q][4]] for q in TRAIN]); y = np.concatenate([D[q][3][D[q][4]] for q in TRAIN])
ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(5)]
pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)
cut = max(CUTS, key=lambda c: fb(pt, y, c, b))
print(f"trained on the 20x runs of {','.join(TRAIN)}: {len(y)} records, cut {cut:.2f}")
for d in sorted(glob.glob(f"{SIM}/T*")):
    f = f"{d}/calls.vcf"
    if not os.path.exists(f): continue
    L, Xc, K = [], [], []
    for l in open(f):
        if l[0] == "#": continue
        g = l.split("\t"); info = dict(kv.split("=", 1) for kv in g[7].split(";") if "=" in kv); m = ALT.search(g[4])
        K.append((g[0], m.group(1) if m else g[0], int(g[1]), int(m.group(2)) if m else int(g[1]) + int(info["SVLEN"]), info["SVTYPE"], info["O"]))
        Xc.append([float(v) for v in info["F"].split(",")]); L.append(l)
    pr = np.mean([m.predict_proba(np.array(Xc))[:, 1] for m in ms], 0); du = mark(K, pr)
    with open(f"{d}/verdict.vcf", "w") as o:
        o.write("".join(l for l in open(f) if l[0] == "#"))
        for l, q, u in zip(L, pr, du):
            g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if q >= cut else "LOW"; g[7] += ";VP=%.4f" % q; o.write("\t".join(g))
    print(os.path.basename(d), "PASS", int(((pr >= cut) & ~du).sum()), "records", len(L))
'''
SCRIPTS["style"] = r'''import re, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

TEX = "/home/yuio/snape/paper/figures/snape2/system/pieces/style.tex"
_s = open(TEX).read()
_pal = re.search(r"\\providecommand\\palette\{(\w+)\}", _s).group(1)
_th = re.search(r"\\newcommand\\theme%s\s*\{\\settheme((?:\{[0-9A-Fa-f]{6}\}){6})\}" % _pal, _s).group(1)
_v = re.findall(r"\{([0-9A-Fa-f]{6})\}", _th)
C = dict(zip(["normal", "tumor", "ink", "ink2", "muted", "shared"], ["#" + h for h in _v]))
C.update({n: "#" + h for n, h in re.findall(r"\\definecolor\{(\w+)\}\{HTML\}\{([0-9A-Fa-f]{6})\}", _s)})
FS = {n: float(v) for n, v in re.findall(r"\\newcommand\{\\font(\w+)\}\{[^}]*?\\fontsize\{([\d.]+)\}", _s)}
GEO = {n: float(v) for n, v in re.findall(r"\\newcommand\{\\(plot\w+)\}\{([\d.]+)\}", _s)}
MM = 1 / 25.4
FONT = next((f for f in ("Helvetica", "Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans")
             if any(f == x.name for x in font_manager.fontManager.ttflist)), "DejaVu Sans")
plt.rcParams.update({
    "font.family": FONT, "font.size": GEO["plotticks"], "axes.labelsize": GEO["plotlabels"], "axes.titlesize": GEO["plotlabels"],
    "xtick.labelsize": GEO["plotticks"], "ytick.labelsize": GEO["plotticks"], "legend.fontsize": GEO["plotticks"],
    "axes.edgecolor": C["ink2"], "axes.labelcolor": C["ink"], "xtick.color": C["ink2"], "ytick.color": C["ink2"],
    "text.color": C["ink"], "axes.linewidth": GEO["plotaxislw"], "xtick.major.width": GEO["plotaxislw"],
    "ytick.major.width": GEO["plotaxislw"], "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": GEO["plotlw"],
    "lines.markersize": GEO["plotms"], "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 300,
})

from matplotlib.colors import LinearSegmentedColormap
SEQ = LinearSegmentedColormap.from_list("seq", [C["seq%d" % i] for i in range(1, 6)])
DIV = LinearSegmentedColormap.from_list("div", [C["divLow"], C["divMid"], C["divHigh"]])
STAGES = [("pack", "stPack"), ("extract", "stExtract"), ("count", "stCount"), ("graph + prune", "stGraph"),
          ("select .. patch", "stAlign"), ("call", "stCall"), ("verdict", "stVerdict")]
STEPSTAGE = {"pack_n": "pack", "pack_t": "pack", "extract_n": "extract", "extract_t": "extract", "count_n": "count",
             "count_t": "count", "graph": "graph + prune", "prune": "graph + prune", "select": "select .. patch",
             "place": "select .. patch", "local": "select .. patch", "lift": "select .. patch", "clip": "select .. patch",
             "gaps": "select .. patch", "patch": "select .. patch", "call": "call", "verdict": "verdict"}
CALLERCOL = {"Severus": C["clSeverus"], "Sniffles2": C["clSniffles"], "SAVANA": C["clSavana"], "nanomonsv": C["clNanomonsv"],
             "SVision-pro": C["clSvision"], "Delly": C["clDelly"], "minisv": C["clMinisv"], "SomaSV": C["clSomasv"],
             "longcallD": C["clLongcalld"], "SVScope": C["clSvscope"]}
PAIRS = ["hcc1954", "h2009", "hcc1937", "h1437", "hg008", "colo829", "hcc1395", "hg009", "h209", "h2126"]
PAIRNAME = {"hcc1954": "HCC1954", "h2009": "NCI-H2009", "hcc1937": "HCC1937", "h1437": "NCI-H1437", "hg008": "HG008",
            "colo829": "COLO829", "hcc1395": "HCC1395", "hg009": "HG009", "h209": "NCI-H209", "h2126": "NCI-H2126"}
PAIRCOL = {p: C["p" + ("H" + p[1:] if p.startswith("h") and p[1].isdigit() else p.upper())] for p in PAIRS}

def figure(w=None, h=None):
    fig = plt.figure(figsize=((w or GEO["plotw"]) * MM, (h or GEO["ploth"]) * MM))
    return fig, fig.add_axes([0.12, 0.13, 0.84, 0.83])

def save(fig, path):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{path}.{ext}", facecolor="white")
    plt.close(fig)
'''
SCRIPTS["accuracy"] = r'''import csv, json, statistics as st
E = "/home/yuio/snape/bench/eval"
P = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]
NAME = {"A": "Severus HiFi", "B": "SAVANA HiFi", "C": "nanomonsv HiFi", "D": "Sniffles2 HiFi", "E": "snape",
        "F": "Severus ONT", "G": "SAVANA ONT", "H": "nanomonsv ONT", "I": "Sniffles2 ONT", "J": "SvABA", "K": "GRIPSS", "L": "Manta"}


def table(p, f):
    r = list(csv.reader(open(f"{E}/{p}/main/results/{f}"), delimiter="\t"))
    return r[0], {x[0][0]: x[1:] for x in r[1:]}


per = {}
for p in P:
    h, t = table(p, "overall.tsv")
    for k, v in t.items():
        per.setdefault(NAME[k], {})[p] = dict(tp=int(v[0]), fn=int(v[1]), fp=int(v[2]), P=float(v[3]), R=float(v[4]), F1=float(v[5]), calls=int(v[6]), truth=int(v[7]))
out = {"per": per, "median": {}, "rank": {}}
for m in ("F1", "P", "R"):
    med = {c: st.median(per[c][p][m] for p in P) for c in per}
    out["median"][m] = dict(sorted(med.items(), key=lambda x: -x[1]))
out["rank_per_pair"] = {p: 1 + sorted([per[c][p]["F1"] for c in per], reverse=True).index(per["snape"][p]["F1"]) for p in P}
cls = ["BND", "50-100", "100-1k", "1k-10k", ">10k"]
size = {}
for p in P:
    for kind in ("TP", "FN", "FP"):
        h, t = table(p, f"SV_length_{kind}.tsv")
        for k, v in t.items():
            if k not in "ABCDE": continue
            v = [int(x) for x in v]
            vals = [v[0], v[1] + v[2], v[3], v[4], v[5]]
            d = size.setdefault(NAME[k], {c: dict(TP=0, FN=0, FP=0) for c in cls})
            for c, x in zip(cls, vals): d[c][kind] += x
out["size"] = {c: {k: round(2 * d["TP"] / max(1, 2 * d["TP"] + d["FP"] + d["FN"]), 3) for k, d in v.items()} for c, v in size.items()}
out["size_truth"] = {k: size["snape"][k]["TP"] + size["snape"][k]["FN"] for k in cls}
json.dump(out, open("/home/yuio/snape/paper/figures/plots/paper/accuracy.json", "w"), indent=1)
for m in ("F1", "P", "R"):
    print(m, {k: round(v, 3) for k, v in out["median"][m].items()})
print("snape rank per pair", out["rank_per_pair"])
print("size truth", out["size_truth"])
for c, v in out["size"].items(): print(c, v)
print("snape per pair", {p: (per["snape"][p]["calls"], per["snape"][p]["truth"], per["snape"][p]["tp"], per["snape"][p]["fp"], per["snape"][p]["fn"], per["snape"][p]["P"], per["snape"][p]["R"], per["snape"][p]["F1"]) for p in P})
'''
SCRIPTS["figs"] = r'''import sys, json, statistics as st
exec(SCRIPTS["style"])

R = "/home/yuio/snape/bench/runs"
OUT = "/home/yuio/snape/paper/figures/plots/paper"
P = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]
A = json.load(open(f"{OUT}/accuracy.json"))
MAP = "#9a9a9a"


def tag(ax, s):
    ax.text(-0.13, 1.04, s, transform=ax.transAxes, fontweight="bold", fontsize=8, va="bottom")


def bars(ax, metric, label):
    order = list(A["median"][metric].keys())
    for i, c in enumerate(order):
        col = C["tumor"] if c == "snape" else ("#6f8fb0" if "HiFi" in c else ("#b8c6d6" if "ONT" in c else "#d9d9d9"))
        ax.bar(i, A["median"][metric][c], color=col, width=0.7, zorder=1)
        for p in P:
            ax.plot(i, A["per"][c][p][metric], "o", mfc="none", mec=PAIRCOL[p], ms=2.6, mew=0.6, zorder=2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=55, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel(label)


fig = plt.figure(figsize=(180 * MM, 150 * MM))
ax = [fig.add_axes(r) for r in ([0.07, 0.62, 0.40, 0.33], [0.57, 0.62, 0.40, 0.33], [0.07, 0.13, 0.40, 0.33], [0.57, 0.13, 0.40, 0.33])]
bars(ax[0], "F1", "F1"); tag(ax[0], "a")
bars(ax[1], "P", "Precision"); tag(ax[1], "b")
bars(ax[2], "R", "Recall"); tag(ax[2], "c")
cls = ["BND", "50-100", "100-1k", "1k-10k", ">10k"]
lab = ["BND", "50–100 bp", "100 bp–1 kb", "1–10 kb", ">10 kb"]
callers = ["snape", "Severus HiFi", "SAVANA HiFi", "nanomonsv HiFi", "Sniffles2 HiFi"]
cols = [C["tumor"], "#6f8fb0", "#c9a227", "#5a9e6f", "#9a9a9a"]
w = 0.16
for k, (c, col) in enumerate(zip(callers, cols)):
    ax[3].bar([i + (k - 2) * w for i in range(5)], [A["size"][c][x] for x in cls], width=w, color=col, label=c.replace(" HiFi", ""))
ax[3].set_xticks(range(5))
ax[3].set_xticklabels([f"{l}\n({A['size_truth'][x]})" for l, x in zip(lab, cls)])
ax[3].set_ylim(0, 1)
ax[3].set_ylabel("F1")
ax[3].legend(ncol=5, loc="upper left", bbox_to_anchor=(0, 1.13), handlelength=1, columnspacing=0.8)
tag(ax[3], "d")
for p in P:
    ax[0].plot([], [], "o", mfc="none", mec=PAIRCOL[p], label=PAIRNAME[p])
ax[0].legend(ncol=5, loc="upper left", bbox_to_anchor=(0, 1.16), handlelength=0.8, columnspacing=0.6)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/fig_accuracy.{ext}", facecolor="white")
plt.close(fig)

FQ = {}
for p in P:
    b = json.load(open(f"{R}/before_signals/full/{p}/run.json"))
    FQ[p] = {x["step"]: x for x in b["steps"] if x["step"].startswith("pack")}
GROUP = [("Packing", ["pack_n", "pack_t"]), ("Core extraction", ["extract_n", "extract_t"]), ("Counting", ["count_n", "count_t"]),
         ("Graph and pruning", ["graph", "prune"]), ("Selection", ["select"]), ("Placement and alignment", ["place", "lift", "clip", "patch"]),
         ("minimap2", ["local", "gaps"]), ("Calling", ["call"])]
GCOL = ["#c9c3b8", "#a7b8c9", "#7f9cb8", "#5577a0", "#d9a3a8", "#c45c66", MAP, "#4f8a5b"]
T, M = {}, {}
for p in P:
    j = json.load(open(f"{R}/full/{p}/run.json"))
    s = {x["step"]: x for x in j["steps"]}
    s.update(FQ[p])
    T[p] = {g: sum(s[k]["wall_s"] for k in ks) for g, ks in GROUP}
    M[p] = {g: max(s[k]["peak_rss_kib"] for k in ks) / 2**20 for g, ks in GROUP}
json.dump(dict(time=T, mem=M), open(f"{OUT}/cost.json", "w"), indent=1)

ALN = {}
for p in P:
    ALN[p] = sum(json.load(open(f"{R}/align/{p}/{s}/run.json"))["total"]["wall_s"] for s in "nt")
PH = {l.split("\t")[0]: float(l.split("\t")[1]) for l in open(f"{R}/phase/hcc1954/times.tsv") if not l.startswith("step")}
phase = sum(PH.values())
CALL = {}
for t, n in [("severus", "Severus"), ("severus_phased", "Severus (phased)"), ("savana", "SAVANA"), ("nanomonsv", "nanomonsv"),
             ("sniffles2", "Sniffles2"), ("svision-pro", "SVision-pro"), ("cutesv", "cuteSV"), ("debreak", "DeBreak"), ("dysgu", "dysgu")]:
    CALL[n] = json.load(open(f"{R}/competitors/{t}/hcc1954/run.json"))["total"]["wall_s"]
json.dump(dict(align=ALN, phase=PH, callers=CALL), open(f"{OUT}/route.json", "w"), indent=1)

fig = plt.figure(figsize=(180 * MM, 160 * MM))
ax = [fig.add_axes(r) for r in ([0.08, 0.58, 0.40, 0.26], [0.58, 0.58, 0.40, 0.33], [0.08, 0.08, 0.40, 0.30], [0.58, 0.08, 0.40, 0.30])]
for i, p in enumerate(P):
    b = 0
    for (g, _), col in zip(GROUP, GCOL):
        v = T[p][g] / 60
        ax[0].bar(i, v, bottom=b, color=col, width=0.65, label=g if i == 0 else None)
        b += v
ax[0].set_xticks(range(5)); ax[0].set_xticklabels([PAIRNAME[p] for p in P], rotation=30, ha="right")
ax[0].set_ylabel("Wall clock (min)")
ax[0].legend(fontsize=5.5, loc="upper left", bbox_to_anchor=(0.0, 1.40), ncol=3, handlelength=1, columnspacing=0.8)
ax[0].set_ylim(0, 19)
tag(ax[0], "a")
gs = [g for g, _ in GROUP]
for i, (g, col) in enumerate(zip(gs, GCOL)):
    ax[1].bar(i, st.median(M[p][g] for p in P), color=col, width=0.65)
    for p in P:
        ax[1].plot(i, M[p][g], "o", mfc="none", mec=PAIRCOL[p], ms=2.6, mew=0.6)
ax[1].set_xticks(range(len(gs))); ax[1].set_xticklabels(gs, rotation=40, ha="right")
ax[1].set_ylabel("Peak memory (GB)")
tag(ax[1], "b")
for i, p in enumerate(P):
    sn = sum(T[p].values()) / 3600
    ax[2].bar(i - 0.18, sn, width=0.34, color=C["tumor"], label="snape" if i == 0 else None)
    ax[2].bar(i + 0.18, ALN[p] / 3600, width=0.34, color=MAP, hatch=None if p == "hcc1954" else "////", edgecolor="white", lw=0,
              label="minimap2, both read sets" if i == 0 else None)
    ax[2].text(i + 0.18, ALN[p] / 3600 + 0.1, f"{ALN[p] / 60 / (sum(T[p].values()) / 60):.0f}×", ha="center", fontsize=5.5)
ax[2].set_xticks(range(5)); ax[2].set_xticklabels([PAIRNAME[p] for p in P], rotation=30, ha="right")
ax[2].set_ylabel("Wall clock (h)")
ax[2].legend(loc="upper left", fontsize=5.5)
tag(ax[2], "c")
names = list(CALL)
sn = sum(T["hcc1954"].values()) / 3600
ys = ["snape"] + names
ax[3].barh(0, sn, color=C["tumor"])
for k, n in enumerate(names, 1):
    a = ALN["hcc1954"] / 3600
    ax[3].barh(k, a, color=MAP)
    ax[3].barh(k, CALL[n] / 3600, left=a, color="#c9c3b8" if "phased" in n else "#5577a0")
from matplotlib.patches import Patch
ax[3].legend(handles=[Patch(color=MAP, label="minimap2, both read sets"), Patch(color="#c9c3b8", label="phasing + Severus (one run)"),
                      Patch(color="#5577a0", label="caller")], loc="upper right", fontsize=5.5)
ax[3].set_xlim(0, 7.5)
ax[3].set_yticks(range(len(ys))); ax[3].set_yticklabels(ys)
ax[3].invert_yaxis()
ax[3].set_xlabel("Wall clock on HCC1954 (h)")
tag(ax[3], "d")
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/fig_cost.{ext}", facecolor="white")
plt.close(fig)
print(json.dumps(dict(T={p: round(sum(T[p].values())) for p in P}, ALN=ALN, phase=phase, CALL=CALL), indent=0))
'''

SCRIPTS["simjoint"] = r'''import sys, os, json, collections, numpy as np
exec(SCRIPTS["covverdict"].split("out = open(")[0])
O = "/home/yuio/snape/bench/runs/sim"; J = f"{O}/joint"
HI, LO, NT, SEEN, JAC, ALPHA = 0.85, 0.5, 5, 20, 0.5, 0.05
D = {p: crecords(f"{C}/{p}_20x/calls.vcf", p) for p in TRAIN}
b = float(json.load(open(f"{C}/hcc1954_20x/monitor/params.json"))["verdict_f_beta"])
X = np.vstack([D[q][1][D[q][4]] for q in TRAIN]); y = np.concatenate([D[q][3][D[q][4]] for q in TRAIN])
ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(5)]
pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)
cut = max(CUTS, key=lambda c: fb(pt, y, c, b))
H0 = [l for l in open(f"{J}/calls.vcf") if l[0] == "#"]
L, Xc, K, DV, SP, SM = [], [], [], [], [], []
for l in open(f"{J}/calls.vcf"):
    if l[0] == "#": continue
    g = l.rstrip("\n").split("\t"); info = dict(kv.split("=", 1) for kv in g[7].split(";") if "=" in kv); m = ALT.search(g[4])
    K.append((g[0], m.group(1) if m else g[0], int(g[1]), int(m.group(2)) if m else int(g[1]) + int(info["SVLEN"]), info["SVTYPE"], info["O"]))
    Xc.append([float(v) for v in info["F"].split(",")]); L.append(g)
    DV.append([int(x.split(":")[1]) for x in g[10:]])
    SP.append([int(x) for x in info["SPAN"].split(",")][1:]); SM.append([int(x) for x in info["SAME"].split(",")][1:])
DV, SP, SM = np.array(DV), np.array(SP), np.array(SM); NTUM = DV.shape[1]
pr = np.mean([m.predict_proba(np.array(Xc))[:, 1] for m in ms], 0); du = mark(K, pr)
strong = (DV >= 2) | (SM >= 2); ns = strong.sum(1); multi = ns >= 2
acc = ~du & np.where(multi, (pr >= HI) | (((SM >= 2).sum(1) >= NT) & (pr >= LO)), pr >= cut)
pc = collections.Counter(tuple(x) for x in strong[(pr >= 0.75) & ~du & multi])
masks = np.array([p for p, n in pc.items() if n >= SEEN and sum(p) >= 2], bool).reshape(-1, NTUM)
fill = np.zeros_like(strong)
if len(masks):
    inside = (strong[:, None, :] & masks[None, :, :]).sum(2); jac = inside / (ns[:, None] + masks.sum(1)[None, :] - inside + 1e-9)
    best = jac.argmax(1); ok = (jac.max(1) >= JAC) & (inside[np.arange(len(L)), best] >= 2)
    frac = np.clip((SM * strong).sum(1) / np.maximum((SP * strong).sum(1), 1), 0.02, 0.98)
    fill = masks[best] & ok[:, None] & ~strong & ~((SM == 0) & ((1 - frac)[:, None] ** SP < ALPHA))
with open(f"{J}/verdict.vcf", "w") as o:
    o.write("".join(H0))
    for i, g in enumerate(L):
        g = list(g); g[6] = "DUPL" if du[i] else "PASS" if acc[i] else "LOW"; g[7] += ";VP=%.4f" % pr[i]
        g[8] = "GT:DV:CR"
        for k in range(NTUM): g[10 + k] = ("0/1" if strong[i, k] or fill[i, k] else "0/0") + ":%d:%s" % (DV[i, k], "O" if strong[i, k] else "I" if fill[i, k] else ".")
        o.write("\t".join(g) + "\n")
ns_ = {}; exec(SCRIPTS["simscore"].split("if __name__")[0], ns_)
tumors = ["T%02d" % i for i in range(1, NTUM + 1)]
print(f"cut {cut:.2f}; records {len(L)}; accepted {int(acc.sum())}; carrier patterns {len(masks)}")
rows = ["tumor\ttruth\tpair_P\tpair_R\tpair_F1\tjoint_P\tjoint_R\tjoint_F1\tinferred_P\tinferred_R\tinferred_F1"]
for k, t in enumerate(tumors):
    out = [ns_["score"](f"{O}/truth/{t}.tsv", f"{O}/pair/{t}/verdict.vcf")]
    for name, car in (("observed", strong), ("inferred", strong | fill)):
        f2 = f"{J}/{t}.{name}.vcf"
        open(f2, "w").write("".join(H0) + "".join("\t".join(L[i][:6] + ["PASS"] + L[i][7:]) + "\n" for i in np.where(acc & car[:, k])[0]))
        out.append(ns_["score"](f"{O}/truth/{t}.tsv", f2))
    rows.append(t + "\t%d\t" % out[0]["truth"] + "\t".join("%.3f\t%.3f\t%.3f" % (r["P"], r["R"], r["F1"]) for r in out))
open(f"{O}/compare.tsv", "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
'''

def run_score(t, c):
    out = io.StringIO(); a = sys.argv
    sys.argv = ["score", t, c]
    with contextlib.redirect_stdout(out): exec(SCRIPTS["score"], {"__name__": "__main__"})
    sys.argv = a
    return out.getvalue()
cmd = sys.argv[1]; sys.argv = [cmd] + sys.argv[2:]
exec(SCRIPTS[cmd], {"__name__": "__main__", "SCRIPTS": SCRIPTS, "run_score": run_score})
PYEOF
read -r -d '' SVB2FQ <<'CEOF'
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    FILE* f = fopen(argv[1], "rb");
    const char* name = argc > 2 ? argv[2] : "r";
    char hdr[32];
    if (!f || fread(hdr, 1, 32, f) != 32 || memcmp(hdr, "SVB1", 4)) { fprintf(stderr, "not an SVB1 file\n"); return 1; }
    static char out[1 << 20];
    setvbuf(stdout, out, _IOFBF, sizeof out);
    const char B[4] = {'A', 'C', 'G', 'T'};
    size_t cap = 1 << 20;
    uint8_t* p = malloc(cap);
    char* s = malloc(cap * 4 + 1);
    char* q = malloc(cap * 4 + 1);
    uint64_t n = 0;
    uint32_t tag;
    while (fread(&tag, 4, 1, f) == 1) {
        uint32_t L = tag & 0x7fffffffu, raw = tag >> 31, nb = raw ? L : (L + 3) >> 2;
        if (nb > cap || L > cap * 4) {
            cap = (nb > L ? nb : L) * 2;
            p = realloc(p, cap); s = realloc(s, cap * 4 + 1); q = realloc(q, cap * 4 + 1);
        }
        if (fread(p, 1, nb, f) != nb) { fprintf(stderr, "truncated at read %lu\n", (unsigned long)n); return 1; }
        if (raw) memcpy(s, p, L);
        else for (uint32_t i = 0; i < L; i++) s[i] = B[(p[i >> 2] >> (6 - 2 * (i & 3))) & 3];
        char* r = s; uint32_t n = L;
        const char* A = "AATGTACTTCGTTCAGTTACGTATTGCT";
        const char* E = "GCAATACGTAACTGAACGAAGTACATT";
        uint32_t cut = 0;
        for (uint32_t j = 0; j + 21 <= n && j < 100; j++) if (!memcmp(r + j, A + 7, 21)) { cut = j + 21; break; }
        if (!cut) for (uint32_t k = 20; k >= 8; k--) if (k <= n && !memcmp(r, A + 28 - k, k)) { cut = k; break; }
        r += cut; n -= cut;
        uint32_t end = n;
        for (uint32_t j = n > 100 ? n - 100 : 0; j + 21 <= n; j++) if (!memcmp(r + j, E, 21)) { end = j; break; }
        if (end == n) for (uint32_t k = 20; k >= 8; k--) if (k <= n && !memcmp(r + n - k, E, k)) { end = n - k; break; }
        n = end;
        if (n < 100) continue;
        memmove(s, r, n); L = n;
        memset(q, 'I', L);
        printf("@%s_%lu\n", name, (unsigned long)n++);
        fwrite(s, 1, L, stdout); fputs("\n+\n", stdout);
        fwrite(q, 1, L, stdout); fputc('\n', stdout);
    }
    fprintf(stderr, "%lu reads\n", (unsigned long)n);
    return 0;
}
CEOF
src() { [ -s $S/hifi/$1/$2.bam ] && echo $S/hifi/$1/$2.bam || echo $S/hifi/$1/$2.fastq.gz; }
depth() { python3 -c "import json; r=json.load(open('$R/full/$1/run.json'))['reads']; print(int(min(r['n']['bases'], r['t']['bases']) // $G))"; }
minreads() { [ $1 -lt 13 ] && echo '"n2 t2"' || { [ $1 -lt 25 ] && echo '"n3 t2"' || echo '"n3 t3"'; }; }
beta() { [ $1 -le 5 ] && echo 0.5 || { [ $1 -le 10 ] && echo 0.75 || echo 1; }; }
DEV=$(basename $(findmnt -no SOURCE /home) | sed 's/p[0-9]*$//'); HW=$(ls -d /sys/block/$DEV/device/hwmon* 2>/dev/null | head -1)

snape_run() {
    local E=$1 D=$2 p=$3 n=$4 th=16 m="" f=1
    [ -n "$n" ] && { m=$(minreads $n); f=$(beta $n); }
    [ -s $D/run.json ] && return
    mkdir -p $D/work $D/monitor $D/logs; log "[$(basename $D)] start"
    for c in $(src $p n) $(src $p t) $D/work/n.svb $D/work/t.svb; do [ -e $c ] && dd if=$c iflag=nocache count=0 status=none; done
    sync; sleep 30
    printf '{"experiment": "%s", "pair": "%s", "coverage_x": %s, "bases_per_x": %s, "bases_limit": null, "min_reads": %s, "verdict_f_beta": %s, "threads": %s, "intermediates": "deleted"}\n' \
        $E $p "${n:-null}" $G "${m:-\"n3 t3\"}" $f $th > $D/monitor/params.json
    pidstat -h -r -u -d -C 'snape|minimap2' 2 >> $D/monitor/pidstat.txt 2>&1 & m1=$!
    iostat -dxm -t $DEV 5 >> $D/monitor/iostat.txt 2>&1 & m2=$!
    vmstat -t 5 >> $D/monitor/vmstat.txt 2>&1 & m3=$!
    ( while :; do echo "$(date +%s) $(cat $HW/temp1_input) $(cat $HW/temp3_input) $(sensors 2>/dev/null | awk '/Tctl/{print $2}')"; sleep 5; done ) >> $D/monitor/temp.txt & m4=$!
    date +%s.%N > $D/start.txt
    /usr/bin/time -f "total %e s wall %M KB peak %U s user %S s sys %I fs-in-blocks %O fs-out-blocks %x exit" \
        taskset -c 0-$((th - 1)) bash $S/pipeline.sh $D/work $(src $p n) $(src $p t) $th 2> $D/total.txt
    kill $m1 $m2 $m3 $m4 2>/dev/null; wait $m1 $m2 $m3 $m4 2>/dev/null
    [ -s $D/work/calls.vcf ] && [ -e $D/work/done.call ] || { log "[$(basename $D)] FAILED: $(tail -1 $D/total.txt)"; return; }
    py measure $D >> $Q 2>&1 || { log "[$(basename $D)] measure failed, intermediates left in place"; return; }
    mv $D/work/calls.vcf $D/; mv $D/work/log.* $D/logs/
    rm -rf $D/work; log "[$(basename $D)] done"
}

measured() {
    local D=$1 P=$2 mem=$3; shift 3
    [ -s $D/run.json ] && return 0
    mkdir -p $D/out; echo "$P" > $D/params.json
    local sl="sb_$(basename $(dirname $D) | tr -c 'a-z0-9\n' '_')_$(basename $D | tr -c 'a-z0-9\n' '_').slice" c pid
    log "[$(basename $(dirname $D))/$(basename $D)] start"
    CGPARENT=$sl DMEM=$mem systemd-run --user --scope --quiet --slice=$sl -p MemoryMax=$mem -p MemorySwapMax=0 \
        /usr/bin/time -f "%e %M %U %S %x" -o $D/time.txt taskset -c 0-15 "$@" > $D/out/stdout.txt 2>&1 &
    pid=$!; c=/sys/fs/cgroup/$sl
    ( m=0; while kill -0 $pid 2> /dev/null; do x=$(awk '/^anon /{a=$2} /^shmem /{s=$2} END{print a+s}' $c/memory.stat 2> /dev/null); [ "${x:-0}" -gt $m ] && { m=$x; echo $m > $D/anon_peak.txt; }; sleep 2; done ) &
    wait $pid
    [ -d $c ] && { echo "memory_peak $(cat $c/memory.peak)"; [ -s $D/anon_peak.txt ] && echo "anon_peak $(cat $D/anon_peak.txt)"; } > $D/cgroup.txt
    py tool $D >> $Q 2>&1 || log "[$(basename $D)] no measurement"
    for f in calls.vcf log.txt stdout.txt; do [ -e $D/out/$f ] && mv $D/out/$f $D/; done
    rm -rf $D/out
}

aligned() {
    local p=$1 s=$2 B=$3 ref=$4
    minimap2 -ax map-hifi -t 16 -Y --secondary=no -R "@RG\tID:${p}_$s\tSM:${p}_$s" $ref $(src $p $s) 2> $B.minimap2.log \
        | samtools sort -@ 4 -m 2G -T $B.st -o $B.tmp - && mv $B.tmp $B && samtools index -@ 8 $B
}

caller() {
    local tool=$1 T=$(readlink -f $2) N=$(readlink -f $3) P=$5; mkdir -p $4; local O=$(readlink -f $4)
    [ $P -gt 16 ] && P=16
    local REF=${REF:-$S/ref/hg38_parmask.fa} VNTR=${VNTR:-$S/ref/vntr.trf.bed} CPUS=${CPUS:-0-15}
    exec > $O/log.txt 2>&1
    echo "[$(date '+%F %T')] $tool tumor=$T normal=$N out=$O threads=$P"
    dk() {
        local img=$1 m=() d; shift; mkdir -p $O/.tmp
        for d in $S $(printf "%s\n" "$(dirname $T)" "$(dirname $N)" "$O" | grep -v "^$S\(/\|$\)" | sort -u); do m+=(-v "$d:$d"); done
        docker run --rm --log-driver none --cpuset-cpus $CPUS --memory ${DMEM:-40g} --memory-swap ${DMEM:-40g} --oom-score-adj 1000 ${CGPARENT:+--cgroup-parent=$CGPARENT} \
            -e OMP_NUM_THREADS=$P -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$O/.tmp:/tmp" -w "$O" "${m[@]}" "$img" "$@"
    }
    pass() { awk -F'\t' '/^#/ || $7=="PASS" || $7=="."' "$1" > "$2"; }
    sub() { py subtract "$@"; }
    case $tool in
    severus) $L/severus/bin/severus --target-bam $T --control-bam $N --out-dir $O/severus --vntr-bed $VNTR -t $P
        pass $O/severus/somatic_SVs/severus_somatic.vcf $O/calls.vcf ;;
    severus_hp) $L/severus/bin/severus --target-bam $T --control-bam $N --phasing-vcf $(dirname $N)/clair3/phased_merge_output.vcf.gz --out-dir $O/severus --vntr-bed $VNTR -t $P
        pass $O/severus/somatic_SVs/severus_somatic.vcf $O/calls.vcf ;;
    savana) $L/savana/bin/python -c "$PY" savana --tumour $T --normal $N --ref $REF --outdir $O/savana --threads $P --pb --sample s --overwrite
        pass $O/savana/s.classified.somatic.vcf $O/calls.vcf ;;
    nanomonsv) local M=$L/nanomonsv/bin/nanomonsv
        sort -k1,1 -k2,2n $VNTR | bgzip > $O/repeat.bed.gz && tabix -p bed $O/repeat.bed.gz; mkdir -p $O/t $O/n
        $M parse $T $O/t/t --reference_fasta $REF; $M parse $N $O/n/n --reference_fasta $REF
        $M get $O/t/t $T $REF --control_prefix $O/n/n --control_bam $N --processes $P --qv25 --simple_repeat_bed $O/repeat.bed.gz
        pass $O/t/t.nanomonsv.result.vcf $O/calls.vcf ;;
    sniffles2) local F=$L/sniffles/bin/sniffles
        $F --input $N --snf $O/n.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
        $F --input $T --snf $O/t.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
        $F --input $O/n.snf $O/t.snf --vcf $O/joint.vcf --threads $P --reference $REF --allow-overwrite
        awk -F'\t' '/^#/ || ($7=="PASS" && $8 ~ /SUPP_VEC=01(;|$)/)' $O/joint.vcf > $O/calls.vcf ;;
    sniffles2_v22) local F="dk quay.io/biocontainers/sniffles:2.2--pyhdfd78af_0 sniffles"
        $F --input $N --snf $O/n.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
        $F --input $T --snf $O/t.snf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
        $F --input $O/n.snf $O/t.snf --vcf $O/joint.vcf --threads $P --reference $REF --allow-overwrite
        awk -F'\t' '/^#/ || ($7=="PASS" && $8 ~ /SUPP_VEC=01(;|$)/)' $O/joint.vcf > $O/calls.vcf ;;
    sniffles2_mosaic) local F=$L/sniffles/bin/sniffles
        $F --input $T --vcf $O/t.mosaic.vcf --mosaic --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite
        $F --input $N --vcf $O/n.vcf --threads $P --reference $REF --tandem-repeats $VNTR --allow-overwrite --minsupport 1
        sub $O/t.mosaic.vcf $O/n.vcf --dist 1000 > $O/calls.vcf ;;
    svision-pro) local I=snape-bench/svision-pro:2.5-pysam0.21 RR=$K/ref/svision-pro NP=$(( P / 2 )); [ $NP -lt 1 ] && NP=1; mkdir -p $O/svp
        dk $I SVision-pro --target_path $T --base_path $N --genome_path $REF --model_path $RR/src/pre_process/model_liteunet_256_8_16_32_32_32.pth \
           --access_path $RR/src/pre_process/hg38.access.10M.bed --out_path $O/svp --sample_name tumor --preset hifi --detect_mode somatic \
           --process_num $NP --unet_cpu_num 2 --device cpu
        local v=$(ls $O/svp/*.vcf | grep -v somatic_s | head -1)
        dk $I python $RR/extract_op.py --input_vcf $v --extract somatic; pass ${v%.vcf}.somatic_s1.vcf $O/calls.vcf ;;
    cutesv) local I=quay.io/biocontainers/cutesv:2.1.4--pyhdfd78af_1 HI="--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"
        mkdir -p $O/t.work $O/n.work
        dk $I cuteSV $T $REF $O/t.vcf $O/t.work $HI -l 50 -S tumor -t $P; dk $I cuteSV $N $REF $O/n.vcf $O/n.work $HI -l 50 -s 2 -S normal -t $P
        sub $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf ;;
    debreak) local I=quay.io/biocontainers/debreak:1.3--h9ee0642_0
        dk $I debreak --bam $T -o $O/t -t $P --min_size 50; dk $I debreak --bam $N -o $O/n -t $P --min_size 50 --min_support 2
        sub $O/t/debreak.vcf $O/n/debreak.vcf --dist 500 > $O/calls.vcf ;;
    dysgu) local I=quay.io/biocontainers/dysgu:1.8.8--py313h5d164f8_0
        dk $I dysgu call --mode pacbio-revio --diploid False -p $P -x -o $O/t.vcf $REF $O/tmp_t $T
        dk $I dysgu call --mode pacbio-revio -p $P -x -o $O/n.vcf $REF $O/tmp_n $N
        dk $I dysgu filter -p $P --normal-vcf $O/n.vcf -o $O/somatic.vcf $O/t.vcf $N
        bcftools view -f PASS -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="TRA" || abs(INFO/SVLEN)>=50' $O/somatic.vcf -Ov -o $O/calls.vcf ;;
    kled) local I=quay.io/biocontainers/kled:1.2.11--he5eddf3_0
        dk $I kled -R $REF --CCS -t $P $T > $O/t.vcf; dk $I kled -R $REF --CCS -t $P --F2 $N > $O/n.vcf
        sub $O/t.vcf $O/n.vcf --dist 1000 > $O/calls.vcf ;;
    svdss) local I=quay.io/biocontainers/svdss:2.1.1--he17396a_0 FMD=$K/ref/svdss/hg38_parmask.fmd
        one() { dk $I SVDSS smooth --reference $REF --bam $1 --threads $P > $O/$2.smoothed.bam; samtools index -@ $P $O/$2.smoothed.bam
            dk $I SVDSS search --index $FMD --bam $O/$2.smoothed.bam --threads $P > $O/$2.sfs.txt
            dk $I SVDSS call --reference $REF --bam $O/$2.smoothed.bam --sfs $O/$2.sfs.txt --threads $P --min-sv-length 50 --min-cluster-weight $3 > $O/$2.vcf; }
        one $T t 4; one $N n 2; sub $O/t.vcf $O/n.vcf --dist 500 > $O/calls.vcf ;;
    delly) local I=quay.io/biocontainers/delly:2.6.0--h3752d28_0
        sm() { samtools view -H $1 | awk -F'\t' '/^@RG/{for(i=2;i<=NF;i++) if($i~/^SM:/){print substr($i,4); exit}}'; }
        cp ${EXCL:-$K/ref/human.hg38.excl.tsv} $O/excl.tsv; printf "%s\ttumor\n%s\tcontrol\n" "$(sm $T)" "$(sm $N)" > $O/samples.tsv
        dk $I delly lr -y pb -g $REF -x $O/excl.tsv -o $O/joint.bcf $T $N
        dk $I delly filter -f somatic -o $O/somatic.bcf -s $O/samples.tsv $O/joint.bcf
        bcftools view -f PASS $O/somatic.bcf -Ov -o $O/calls.vcf ;;
    *) echo "unknown caller $tool"; exit 2 ;;
    esac
    echo "[$(date '+%F %T')] done: $(grep -vc '^#' $O/calls.vcf) calls in $O/calls.vcf"
}

case $X in
ref)
    cd $S/ref
    $S/build/snape refidx -i hg38_parmask.fa -o hg38_parmask.cidx -t 16
    minimap2 -x map-hifi -d hg38_parmask.mmi hg38_parmask.fa; minimap2 -x sr -d hg38_parmask.sr.mmi hg38_parmask.fa
    $S/build/snape refidx -i chm13v2.fa -o chm13v2.cidx -t 16; minimap2 -x sr -d chm13v2.sr.mmi chm13v2.fa
    mkdir -p $K/ref/svdss; docker run --rm --log-driver none -u "$(id -u):$(id -g)" -v $S:$S quay.io/biocontainers/svdss:2.1.1--he17396a_0 \
        SVDSS index --reference $S/ref/hg38_parmask.fa --index $K/ref/svdss/hg38_parmask.fmd --threads 16 && touch $K/ref/svdss/index.done ;;
fetch)
    bash -s <<'FETCHEOF'
#!/usr/bin/env bash
set -uo pipefail
H=/home/yuio/snape/hifi; RAW=$H/raw; T=$H/truth; TK=$H/tools/sratoolkit/bin; cd "$H"; mkdir -p "$RAW/tmp" "$T"
SDL="https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?location=s3.us-east-1&type=src&acc="
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic
PB=https://downloads.pacbcloud.com/public/revio
DL=4; CV=3; RESERVE=$((150 * 1000 * 1000 * 1000))

JOBS="
sra hcc1954 t SRR28305163
sra hcc1954 n SRR28305160
sra h2009 t SRR28305172
sra h2009 n SRR28305169
sra hcc1937 t SRR28305185
sra hcc1937 n SRR28305182
sra h1437 t SRR28305179
sra h1437 n SRR28305175
giab hg008 t $GIAB/HG008/Liss_lab/BCM_Revio_20240313/HG008-T_ubams $GIAB/HG008/Liss_lab/BCM_Revio_20240313/checksum.md5 m84101_240226_212100_s2.hifi_reads.bam m84101_240226_232018_s3.hifi_reads.bam
giab hg008 n $GIAB/HG008/Liss_lab/BCM_Revio_20240313/HG008-N-D_ubams $GIAB/HG008/Liss_lab/BCM_Revio_20240313/checksum.md5 m84059_240304_183205_s3.hifi_reads.bam m84059_240304_203144_s4.hifi_reads.bam
pacb hcc1395 t $PB/2023Q2/HCC1395/HCC1395/m84039_230414_235240_s2.hifi_reads.default.bam $PB/2023Q2/HCC1395/HCC1395/m84039_230415_002321_s3.hifi_reads.default.bam
pacb hcc1395 n $PB/2023Q2/HCC1395/HCC1395-BL/m84039_230412_183900_s2.hifi_reads.bc2002.bam $PB/2023Q2/HCC1395/HCC1395-BL/m84039_230415_005427_s4.hifi_reads.bc2002.bam
pacb colo829 t $PB/2023Q2/COLO829/COLO829/m84039_230312_025934_s1.hifi_reads.bc2026.bam $PB/2023Q2/COLO829/COLO829/m84039_230328_000836_s3.hifi_reads.bc2025.bam
pacb colo829 n $PB/2023Q2/COLO829/COLO829-BL/m84039_230327_230708_s1.hifi_reads.bc2007.bam $PB/2023Q2/COLO829/COLO829-BL/m84039_230327_233730_s2.hifi_reads.bc2007.bam
"

free() { df --output=avail -B1 "$H" | tail -1; }
log() { echo "$(date '+%F %T') $*"; }
slot() {
    local k=$1 n=$2 i fd; [ $k = dl ] && fd=8 || fd=7
    while :; do for i in $(seq 0 $((n - 1))); do eval "exec $fd> $RAW/tmp/$k.$i.lock"; flock -n $fd && return 0; done; sleep 20; done
}
unslot() { [ $1 = dl ] && { flock -u 8; exec 8>&-; } || { flock -u 7; exec 7>&-; }; }

etag_ok() {
    python3 - "$1" "$2" <<'EOF'
import hashlib, math, os, sys
f, tag = sys.argv[1], sys.argv[2].strip('"'); want, n = tag.split('-'); n = int(n); size = os.path.getsize(f)
for mib in (8, 16, 32, 64, 128, 256, 512):
    ps = mib << 20
    if math.ceil(size / ps) != n: continue
    d = b''
    with open(f, 'rb') as h:
        while True:
            b = h.read(ps)
            if not b: break
            d += hashlib.md5(b).digest()
    sys.exit(0 if hashlib.md5(d).hexdigest() == want else 1)
sys.exit(2)
EOF
}

get() {
    local f=$1 url=$2 bytes=$3 md5=$4 tries=0 ck=()
    [ -s "$f" ] && [ ! -e "$f.aria2" ] && return 0
    [ -n "$md5" ] && ck=(--checksum=md5="$md5")
    [ -e "$f.aria2" ] || while [ "$(free)" -lt $(( bytes + RESERVE )) ]; do log "[$f] waiting for disk"; sleep 600; done
    until aria2c -x16 -s16 -k20M --file-allocation=falloc --console-log-level=warn --summary-interval=60 \
          --show-console-readout=false --max-tries=0 --retry-wait=20 --auto-file-renaming=false "${ck[@]}" \
          -d "$(dirname "$f")" -o "$(basename "$f")" "$url" >> "$RAW/tmp/$(basename "$f").dl.log" 2>&1; do
        tries=$((tries + 1)); log "[$f] aria2c failed ($tries), resuming from the control file"
        [ $tries -ge 30 ] && { log "[FATAL] $f: 30 attempts"; return 1; }; sleep 30
    done
}

bams_to_fastq() {
    local p=$1 s=$2 g=$H/$1/$2.fastq.gz n=0 got b; shift 2
    slot cv $CV
    while [ "$(free)" -lt $(( $(du -cB1 "$@" | tail -1 | cut -f1) * 12 / 10 + RESERVE / 3 )) ]; do log "[$p $s] waiting for disk to convert"; sleep 600; done
    log "[$p $s] converting"
    for b in "$@"; do samtools fastq -@ 4 "$b" 2>> "$RAW/$p/$s.convert.log"; done | pigz -p 8 > "$g.tmp"
    unslot cv
    for b in "$@"; do n=$((n + $(samtools view -c -@ 4 "$b"))); done
    got=$(( $(pigz -dc "$g.tmp" | wc -l) / 4 ))
    if [ "$n" -gt 0 ] && [ "$n" = "$got" ]; then
        mv "$g.tmp" "$g"; rm -f -- "$@" "${@/%/.ok}"; log "[$p $s] done: $got reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $p $s: $n reads in the BAMs, $got in the fastq.gz; BAMs kept"; return 1; fi
}

pacb() {
    local p=$1 s=$2 u f tag bytes files=(); shift 2
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"
    for u in "$@"; do
        f=$RAW/$p/$s.$(basename "$u"); files+=("$f"); [ -e "$f.ok" ] && continue
        tag=$(curl -sSI --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="etag:"{print $2}')
        bytes=$(curl -sSI --retry 5 "$u" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
        slot dl $DL; get "$f" "$u" "$bytes" "" || { unslot dl; return 1; }; unslot dl
        etag_ok "$f" "$tag" && { touch "$f.ok"; log "[$f] downloaded, ETag ok"; } || { log "[FATAL] $f: ETag $tag does not match"; return 1; }
    done
    bams_to_fastq "$p" "$s" "${files[@]}"
}

giab() {
    local p=$1 s=$2 dir=$3 sums m f md5 bytes files=(); sums=$(curl -sS --retry 5 "$4"); shift 4
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"
    for m in "$@"; do
        f=$RAW/$p/$s.$m; files+=("$f")
        md5=$(printf '%s\n' "$sums" | grep -F "$m" | grep -v "\.pbi" | awk '{print $1}' | head -1)
        bytes=$(curl -sSI --retry 5 "$dir/$m" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
        [ -n "$md5" ] && [ -n "$bytes" ] || { log "[FATAL] $p $s $m: no md5/size from GIAB"; return 1; }
        slot dl $DL; get "$f" "$dir/$m" "$bytes" "$md5" || { unslot dl; return 1; }; unslot dl
        log "[$f] downloaded, md5 ok"
    done
    bams_to_fastq "$p" "$s" "${files[@]}"
}

sra() {
    local p=$1 s=$2 acc=$3 f g=$H/$1/$2.fastq.gz meta url md5 bytes r w
    [ -s "$g" ] && { log "[$p $s] present"; return 0; }
    mkdir -p "$H/$p" "$RAW/$p"; f=$RAW/$p/$s.$acc.sra
    meta=$(curl -sS --retry 5 "$SDL$acc")
    url=$(printf '%s' "$meta" | grep -o '"link": "[^"]*"' | head -1 | cut -d'"' -f4)
    md5=$(printf '%s' "$meta" | grep -o '"md5": "[^"]*"' | head -1 | cut -d'"' -f4)
    bytes=$(printf '%s' "$meta" | grep -o '"size": [0-9]*' | head -1 | grep -o '[0-9]*$')
    [ -n "$url" ] && [ -n "$md5" ] && [ -n "$bytes" ] || { log "[FATAL] $acc: no link/md5/size from NCBI"; return 1; }
    slot dl $DL; get "$f" "$url" "$bytes" "$md5" || { unslot dl; return 1; }; unslot dl
    log "[$f] downloaded, md5 ok"
    slot cv $CV
    while [ "$(free)" -lt $(( bytes * 13 / 10 + RESERVE / 3 )) ]; do log "[$p $s] waiting for disk to convert"; sleep 600; done
    log "[$p $s] converting"
    "$TK/fastq-dump" --stdout --split-spot "$f" 2> "$RAW/$p/$s.convert.log" | pigz -p 8 > "$g.tmp"
    unslot cv
    r=$(grep -o '^Read [0-9]* spots' "$RAW/$p/$s.convert.log" | grep -o '[0-9]*')
    w=$(grep -o '^Written [0-9]* spots' "$RAW/$p/$s.convert.log" | grep -o '[0-9]*')
    if [ -n "$r" ] && [ "$r" = "$w" ] && pigz -t "$g.tmp"; then
        mv "$g.tmp" "$g"; rm -f -- "$f"; log "[$p $s] done: $r reads, $(du -BG --apparent-size "$g" | cut -f1)"
    else log "[FATAL] $p $s: conversion check failed (read $r, written $w); $f kept"; return 1; fi
}

sras() {
    local p=$1 s=$2 i=0 a parts=(); shift 2
    [ -s "$H/$p/$s.fastq.gz" ] && { log "[$p $s] present"; return 0; }
    for a in "$@"; do i=$((i + 1)); sra "$p" "${s}_part$i" "$a" || return 1; parts+=("$H/$p/${s}_part$i.fastq.gz"); done
    cat "${parts[@]}" > "$H/$p/$s.fastq.gz.tmp" && pigz -t "$H/$p/$s.fastq.gz.tmp" && mv "$H/$p/$s.fastq.gz.tmp" "$H/$p/$s.fastq.gz" || { log "[FATAL] $p $s: joining failed"; return 1; }
    mkdir -p /home/yuio/trash/fetch-parts && mv "${parts[@]}" /home/yuio/trash/fetch-parts/ && log "[$p $s] done: ${#parts[@]} runs joined, parts in ~/trash/fetch-parts"
}
zen() {
    mkdir -p "$T/$2"
    curl -s "https://zenodo.org/api/records/$1" | python3 -c "
import json,sys
for f in json.load(sys.stdin)['files']: print(f['key'], f['checksum'].split(':')[1], f['links']['self'])" |
    while read -r k md5 url; do
        [ -s "$T/$2/$k" ] && [ "$(md5sum "$T/$2/$k" | cut -d' ' -f1)" = "$md5" ] && continue
        aria2c -x16 -s16 -k1M --file-allocation=none --console-log-level=warn --summary-interval=0 --max-tries=0 \
               --retry-wait=10 --auto-file-renaming=false --allow-overwrite=true -d "$T/$2" -o "$k" "$url" > /dev/null
        [ "$(md5sum "$T/$2/$k" | cut -d' ' -f1)" = "$md5" ] && log "[truth $2] $k md5 ok" || log "[FATAL] truth $2 $k md5 mismatch"
    done
}
web() {
    local d=$1 u; shift; mkdir -p "$T/$d"
    for u in "$@"; do [ -s "$T/$d/$(basename "$u")" ] || curl -sSfL --retry 5 -o "$T/$d/$(basename "$u")" "$u" || log "[FATAL] truth $d $u"; done
}
truth() {
    zen 5019338 colo829
    zen 21231495 prcgap
    web hcc1395 $(for n in 3 4 5 15; do echo "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-022-02816-6/MediaObjects/13059_2022_2816_MOESM${n}_ESM.xlsx"; done)
        echo "$GIAB/HG009/NIST/HG009-T_bulk/20250715p16/analysis/GRCh38/BCM_Revio_HG009T-p16_20260313/HG009-T_vs_HG009-N-WT-p4/structural_variants/$f"; done)
    log "[truth] done"
}

{
log "########## start ${*:-all} ##########"
[ $# -eq 0 ] && truth &
while read -r fn p rest; do
    [ -z "$fn" ] && continue
    [ $# -gt 0 ] && [[ " $* " != *" $p "* ]] && continue
    $fn $p $rest &
    sleep 3
done <<< "$JOBS"
wait
log "########## end ##########"
} 2>&1 | tee -a "$H/fetch.log"
FETCHEOF
    ;;
align)
    for p in ${*:-$CASTLE}; do for s in n t; do
        measured $R/align/$p/$s "{\"experiment\": \"align\", \"tool\": \"minimap2+samtools\", \"pair\": \"$p\", \"sample\": \"$s\", \"bam\": \"$S/hifi/$p/$s.bam\"}" 60G \
            bash -c "$(declare -f src aligned); S=$S aligned $p $s $S/hifi/$p/$s.bam $S/ref/hg38_parmask.mmi"
    done; done ;;
full)
    for p in ${*:-$LABELED}; do snape_run full $R/full/$p $p ""; done ;;
verdict)
    ( cd $H && OMP_NUM_THREADS=16 taskset -c 0-15 $V -c "$PY" verdict > $H/verdict/stdout.txt 2>&1 ); cat $H/verdict/log.txt ;;
truth)
    C=$S/hifi/castle-calls/variant_calls_and_benchmarks; T=$S/hifi/truth; M=$HOME/.local/bin/minda
    for p in ${*:-$CASTLE}; do
        P=${p^^}; D=$C/$P/$P; SA=$P; GR=${P}_subsample.bam.gripss.filtered_vaf.vcf
        [ $p = hcc1395 ] && { SA=H1395; GR=H1395_realigned.bam.gripss.filtered_vaf.vcf; }
        SH=$D/savana_$SA.haplotagged.classified.somatic_HIFI_vaf.vcf; SO=$D/savana_$SA.haplotagged.classified.somatic_ONT_vaf.vcf
        SR="$D/svaba_merged.vcf $D/$GR $D/manta_somaticSV_vaf.vcf"; mkdir -p $T/castle_without_snape $T/castle_without_severus
        $M ensemble --vcfs $D/severus_HIFI_vaf.vcf $SH $D/tumor.nanomonsv.result_HIFI_vaf.vcf $D/sniffles2_HIFI_vaf.vcf \
            $D/severus_ont_vaf.vcf $SO $D/tumor.nanomonsv.result_ONT_vaf.vcf $D/sniffles2_ONT_vaf.vcf $SR \
            --conditions "[[[caller_names[:4], caller_names[4:8], caller_names[8:]], '>=', 2],'&', [caller_names,'>=', 4]]" \
            --min_size 50 --tolerance 500 --out_dir $T/castle_without_snape/$p > $T/castle_without_snape/$p.log 2>&1
        $M ensemble --vcfs $SH $D/tumor.nanomonsv.result_HIFI_vaf.vcf $D/sniffles2_HIFI_vaf.vcf \
            $SO $D/tumor.nanomonsv.result_ONT_vaf.vcf $D/sniffles2_ONT_vaf.vcf $SR \
            --conditions "[[[caller_names[:3], caller_names[3:6], caller_names[6:9]], '>=', 2],'&', [caller_names,'>=', 4]]" \
            --min_size 50 --tolerance 500 --out_dir $T/castle_without_severus/$p > $T/castle_without_severus/$p.log 2>&1
    done ;;
eval)
    C=$S/hifi/castle-calls/variant_calls_and_benchmarks; M=$HOME/.local/bin/minda; : > $R/minda.txt
    for p in ${*:-$CASTLE hg008}; do
        E=$H/eval/$p; O=$E/main; mkdir -p $E; grep -E '^#|PASS' $H/verdict/$p.vcf > $O.vcf
        if [ $p = hg008 ]; then
            $M truthset --base $S/hifi/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz --vcfs $O.vcf --min_size 50 --tolerance 500 --out_dir $O > $O.log 2>&1
            awk -F'\t' 'NR==2{printf "%s\tsnape\tA_nist\tF1 %.4f\tR %.4f\tP %.4f\n","'$p'",$7,$6,$5}' $O/results/overall.tsv >> $R/minda.txt
        else
            P=${p^^}; D=$C/$P/$P; SA=$P; GR=${P}_subsample.bam.gripss.filtered_vaf.vcf
            [ $p = hcc1395 ] && { SA=H1395; GR=H1395_realigned.bam.gripss.filtered_vaf.vcf; }
            $M ensemble --vcfs $D/severus_HIFI_vaf.vcf $D/savana_$SA.haplotagged.classified.somatic_HIFI_vaf.vcf $D/tumor.nanomonsv.result_HIFI_vaf.vcf $D/sniffles2_HIFI_vaf.vcf $O.vcf \
                $D/severus_ont_vaf.vcf $D/savana_$SA.haplotagged.classified.somatic_ONT_vaf.vcf $D/tumor.nanomonsv.result_ONT_vaf.vcf $D/sniffles2_ONT_vaf.vcf \
                $D/svaba_merged.vcf $D/$GR $D/manta_somaticSV_vaf.vcf \
                --conditions "[[[caller_names[:5], caller_names[5:9], caller_names[9:]], '>=', 2],'&', [caller_names,'>=', 4]]" --min_size 50 --tolerance 500 --out_dir $O > $O.log 2>&1
            awk -F'\t' '$1~/^E_/{printf "%s\tsnape\tA_fair\tF1 %.4f\tR %.4f\tP %.4f\n","'$p'",$7,$6,$5}' $O/results/overall.tsv >> $R/minda.txt
        fi
    done; cat $R/minda.txt ;;
coverage)
    E=$R/coverage
    for p in ${*:-$CASTLE}; do
        M=$E/mix/$p; mkdir -p $M
        for s in n t; do [ -s $M/$s.svb.idx ] || taskset -c 0-15 $S/build/snape pack -i $(src $p $s) -o $M/$s.svb -t 16 -s 1 > $M/pack_$s.log 2>&1; done
        for n in 5 10 15 20 30; do
            [ $n -le $(depth $p) ] || continue
            D=$E/${p}_${n}x; [ -s $D/run.json ] && continue; mkdir -p $D/work
            for s in n t; do py mix $D/work/$s.svb $M/$s.svb:0:$((n * G)) >> $Q 2>&1 || continue 2; done
            touch $D/work/done.pack_n $D/work/done.pack_t
            snape_run coverage $D $p $n
        done
        rm -rf $M
    done
    ( cd $H && OMP_NUM_THREADS=16 $V -c "$PY" covverdict > $H/verdict/coverage_stdout.txt 2>&1 ) ;;
purity)
    E=$R/purity
    for p in ${*:-hcc1954 hcc1395}; do
        M=$E/mix/$p; mkdir -p $M
        for s in n t; do [ -s $M/$s.svb.idx ] || taskset -c 0-15 $S/build/snape pack -i $(src $p $s) -o $M/$s.svb -t 16 -s 1 > $M/pack_$s.log 2>&1; done
        nb=$(python3 -c "import struct; print(struct.unpack_from('<Q', open('$M/n.svb', 'rb').read(32), 16)[0])")
        nm=$((nb - 24 * G - G / 10)); [ $nm -gt $((30 * G)) ] && nm=$((30 * G))
        [ -s $M/nm.svb.idx ] || py mix $M/nm.svb $M/n.svb:0:$nm >> $Q 2>&1 || continue
        for q in 20 40 60 80 100; do
            [ -s $M/t$q.svb.idx ] || py mix $M/t$q.svb $M/t.svb:0:$((q * 30 * G / 100)) $([ $q -lt 100 ] && echo $M/n.svb:$nm:$(((100 - q) * 30 * G / 100))) >> $Q 2>&1 || continue
            D=$E/${p}_${q}p; mkdir -p $D/work
            for f in svb svb.idx; do ln -sf $M/nm.$f $D/work/n.$f; ln -sf $M/t$q.$f $D/work/t.$f; done
            touch $D/work/done.pack_n $D/work/done.pack_t
            snape_run purity $D $p ""
        done
        rm -rf $M
    done
    ( cd $H && OMP_NUM_THREADS=16 $V -c "$PY" purityverdict > $H/verdict/purity_stdout.txt 2>&1 ) ;;
phase)
    IMG=hkubal/clair3:v2.0.3; RF=$S/ref/hg38_parmask.fa
    for p in ${*:-$CASTLE}; do
        W=$S/hifi/$p; O=$R/phase/$p; T=$O/times.tsv; mkdir -p $O
        [ -s $T ] || echo -e "step\twall_s\tcpu_s\tpeak_gib\texit" > $T
        step() {
            local name=$1 sl=sb_phase_${p}_$1.slice; shift
            grep -q "^$name	" $T && return 0
            CG=$sl /usr/bin/time -f "%e %U %S %M %x" -o $O/$name.time "$@" > $O/$name.log 2>&1
            local c=/sys/fs/cgroup/$sl cpu=0 peak=0
            [ -d $c ] && { cpu=$(awk '/^usage_usec/{print $2/1e6}' $c/cpu.stat); peak=$(awk '/^anon /{print $2}' $c/memory.stat); }
            read wall us sy rss ex < $O/$name.time
            awk -v n=$name -v w=$wall -v c=$cpu -v u=$us -v y=$sy -v m=$peak -v r=$rss -v x=$ex 'BEGIN{r *= 1024; if (r > m) m = r; printf "%s\t%s\t%.1f\t%.2f\t%s\n", n, w, c + u + y, m / 1073741824, x}' >> $T
        }
        dkc() { docker run --rm --log-driver none --cpuset-cpus 0-15 --cgroup-parent=$CG -u "$(id -u):$(id -g)" -v $S:$S $IMG "$@"; }
        export -f dkc; export S IMG
        step clair3 bash -c "dkc /opt/bin/run_clair3.sh --bam_fn=$W/n.bam --ref_fn=$RF --threads=16 --platform=hifi --model_path=/opt/models/hifi_revio --output=$W/clair3 --enable_phasing --longphase_for_phasing"
        for s in n t; do
            [ -s $W/$s.tagged ] && continue
            step haplotag_$s bash -c "dkc /opt/bin/longphase haplotag -s $W/clair3/phased_merge_output.vcf.gz -b $W/$s.bam -r $RF -t 16 --tagSupplementary -o $W/$s.hp"
            samtools view -@ 8 -b -o $W/$s.unmapped.bam $W/$s.bam '*'
            samtools cat -o $W/$s.hpu.bam $W/$s.hp.bam $W/$s.unmapped.bam && mv $W/$s.hpu.bam $W/$s.hp.bam
            step index_$s taskset -c 0-15 samtools index -@ 16 $W/$s.hp.bam
            a=$(samtools idxstats $W/$s.bam | awk '{m+=$3; u+=$4} END{print m, u}'); b=$(samtools idxstats $W/$s.hp.bam | awk '{m+=$3; u+=$4} END{print m, u}')
            [ "$a" = "$b" ] || { log "$p $s: read counts differ ($a vs $b); original kept"; continue; }
            mv $W/$s.hp.bam $W/$s.bam; mv $W/$s.hp.bam.bai $W/$s.bam.bai; rm -f $W/$s.unmapped.bam; echo "$a" > $W/$s.tagged
        done
    done ;;
competitors)
    for t in ${1:-$TOOLS}; do for p in ${2:-$CASTLE}; do
        mem=60G; th=16; [ $t = savana ] && { mem=150G; th=8; }; [ $t = severus_hp ] && [ ! -s $S/hifi/$p/n.tagged ] && continue
        n=$t; [ $t = severus_hp ] && n=severus_phased
        measured $R/competitors/$n/$p "{\"experiment\": \"competitor\", \"tool\": \"$n\", \"pair\": \"$p\", \"threads\": $th}" $mem \
            bash $H/paper.sh caller $t $S/hifi/$p/t.bam $S/hifi/$p/n.bam $R/competitors/$n/$p/out $th
    done; done
    for p in ${2:-$CASTLE}; do [ -s $R/competitors/sniffles2_v22/$p/calls.vcf ] || { bash $H/paper.sh caller sniffles2_v22 $S/hifi/$p/t.bam $S/hifi/$p/n.bam $R/competitors/sniffles2_v22/$p/out 16; mv $R/competitors/sniffles2_v22/$p/out/calls.vcf $R/competitors/sniffles2_v22/$p/; }; done ;;
table)
    ( cd $H && py table > $R/competitors_table.tsv ); column -t $R/competitors_table.tsv ;;
sim)
    SP=/home/samples/davp-samples/simulation/phylo; O=$R/sim; FQ=$O/fq; mkdir -p $FQ
    [ -x $O/svb2fq ] || { echo "$SVB2FQ" > $O/svb2fq.c; gcc -O2 -o $O/svb2fq $O/svb2fq.c; }
    tofq() { [ -s $FQ/$1.fastq.gz ] || { $O/svb2fq $SP/reads/$1/$1.svb $1 2> $FQ/$1.log | pigz -p 16 -3 > $FQ/$1.fastq.gz.part && mv $FQ/$1.fastq.gz.part $FQ/$1.fastq.gz; }; }
    [ -s $O/truth/tree.json ] || py simtruth $O/truth
    tofq germline
    for t in $(seq -f "T%02g" 1 20); do
        d=$O/pair/$t; [ -s $d/calls.vcf ] && continue
        tofq $t; mkdir -p $d; log "snape $t start"
        REF=chm13v2 /usr/bin/time -f "%e %M %U %S %x" -o $d/time.txt taskset -c 0-15 bash $S/pipeline.sh $d/work $FQ/germline.fastq.gz $FQ/$t.fastq.gz 16 > $d/pipeline.log 2>&1
        cp $d/work/steps.txt $d/; mkdir -p $d/logs; cp $d/work/log.* $d/logs/; [ -s $d/work/calls.vcf ] && cp $d/work/calls.vcf $d/
        rm -rf $d/work $FQ/$t.fastq.gz; log "snape $t done: $(cat $d/time.txt)"
    done
    ( cd $H && $V -c "$PY" simverdict > $O/verdict.txt 2>&1 )
    for t in $(seq -f "T%02g" 1 20); do echo -e "$t\t$(py simscore $O/truth/$t.tsv $O/pair/$t/verdict.vcf)"; done > $O/scores.tsv; cat $O/scores.tsv ;;
simjoint)
    SP=/home/samples/davp-samples/simulation/phylo; O=$R/sim; FQ=$O/fq; J=$O/joint; B=${SNAPE:-$S/build/snape}; RF=$S/ref/chm13v2; mkdir -p $J/work $FQ
    [ -x $O/svb2fq ] || { echo "$SVB2FQ" > $O/svb2fq.c; gcc -O2 -o $O/svb2fq $O/svb2fq.c; }
    SM="germline $(seq -f 'T%02g' 1 20 | tr '\n' ' ')"
    for s in $SM; do [ -s $FQ/$s.fastq.gz ] || { $O/svb2fq $SP/reads/$s/$s.svb $s 2> $FQ/$s.log | pigz -p 16 -3 > $FQ/$s.fastq.gz.part && mv $FQ/$s.fastq.gz.part $FQ/$s.fastq.gz; }; done
    jp() {
        cd $J/work
        st() { n=$1; shift; [ -e done.$n ] && return; /usr/bin/time -f "$n %e s %M KB %U user %S sys" "$@" > log.$n 2>&1 || exit 1; tail -1 log.$n >> steps.txt; touch done.$n; }
        local svb=$(for s in $SM; do printf "%s.svb," $s; done | sed 's/,$//') cnt=$(for s in $SM; do printf "%s.cnt," $s; done | sed 's/,$//')
        for s in $SM; do st pack_$s $B pack -i $FQ/$s.fastq.gz -o $s.svb -t 16; done
        for s in $SM; do st extract_$s $B extract -i $s.svb -o $s.svb.exc -t 16; done
        for s in $SM; do st count_$s $B count -i $s.svb -o $s.cnt -t 16 -m 2; done
        st graph $B graph -i $cnt -o graph.gph -m 2
        st prune $B prune -i graph.gph -r $svb -o prune -t 16 -m 2
        st select $B select -i graph.gph -d prune.sites -r $svb -o ev -t 16
        st place $B place -i ev.fa -r $RF.cidx -d $RF.fa -o ev -t 16
        st local minimap2 -cx map-hifi --secondary=no -t 16 -o local.paf ev.win.fa ev.on.fa
        st lift $B lift -i local.paf -d $RF.fa -o aln.paf
        st clip $B clip -i ev -r aln.paf -o ev.gaps.fa
        st gaps minimap2 -cx sr --secondary=no -t 16 -o ev.gaps.raw.paf $RF.sr.mmi ev.gaps.fa
        st patch $B patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf
        st call $B call -i ev -r $svb -d $RF.fa -g graph.gph -o calls.vcf -v $RF.trf.bed -t 16 -m 2
    }
    if [ ! -s $J/calls.vcf ]; then
        log "joint start ($B)"
        /usr/bin/time -f "%e %M %U %S %x" -o $J/time.txt taskset -c 0-15 bash -c "$(declare -f jp); B=$B J=$J FQ=$FQ RF=$RF SM='$SM' jp" > $J/run.log 2>&1
        cp $J/work/steps.txt $J/; mkdir -p $J/logs; cp $J/work/log.* $J/logs/; [ -s $J/work/calls.vcf ] && cp $J/work/calls.vcf $J/
        log "joint done: $(cat $J/time.txt)"
    fi
    ( cd $H && $V -c "$PY" simjoint > $J/scores.txt 2>&1 ); cat $J/scores.txt ;;
figures)
    F=$S/paper/figures/plots/paper; mkdir -p $F $S/paper/overleaf/figures
    ( cd $F && py accuracy && py figs )
    cp $S/paper/figures/snape2/system/pieces/20_system.pdf $S/paper/overleaf/figures/fig_system.pdf
    cp $F/fig_accuracy.pdf $F/fig_cost.pdf $S/paper/overleaf/figures/
    docker run --rm --log-driver none -v $S/paper/overleaf:/w -w /w -u "$(id -u):$(id -g)" -e HOME=/tmp latex-tikz \
        sh -c 'pdflatex -interaction=nonstopmode main.tex > /dev/null; bibtex main > /dev/null; pdflatex -interaction=nonstopmode main.tex > /dev/null; pdflatex -interaction=nonstopmode main.tex | grep -E "^!|undefined"' ;;
guard)
    case $1 in
    memory) while :; do a=$(awk '/MemAvailable/{print int($2 / 1048576)}' /proc/meminfo)
            [ "$a" -lt 14 ] && { echo "$(date '+%F %T') memguard: $a GB available; killing containers" >> $R/night/INCIDENTS.md; c=$(docker ps -q); [ -n "$c" ] && docker kill $c > /dev/null; sleep 90; }
            sleep 4; done ;;
    disk) held=0; while :; do a=$(df --output=avail -BG /home | tail -1 | tr -dc 0-9)
            if [ $held = 0 ] && [ $a -lt 60 ]; then c=$(docker ps -q); [ -n "$c" ] && docker pause $c > /dev/null; echo "$(date '+%F %T') DISK GUARD: $a GB free; containers paused" >> $R/night/INCIDENTS.md; held=1
            elif [ $held = 1 ] && [ $a -ge 120 ]; then [ -n "$c" ] && docker unpause $c > /dev/null; echo "$(date '+%F %T') DISK GUARD: $a GB free, resumed" >> $R/night/INCIDENTS.md; held=0; fi
            sleep 60; done ;;
    load) O=$R/load.tsv; [ -s $O ] || echo -e "time\tload1\tall_cores" > $O
          cpu() { awk '/^cpu /{print $2 + $3 + $4 + $7 + $8}' /proc/stat; }; a=$(cpu); t=$(date +%s.%N); h=$(getconf CLK_TCK)
          while :; do sleep 30; b=$(cpu); u=$(date +%s.%N)
              awk -v a=$a -v b=$b -v t=$t -v u=$u -v h=$h -v l="$(cut -d' ' -f1 /proc/loadavg)" -v d="$(date '+%F %T')" 'BEGIN{printf "%s\t%s\t%.1f\n", d, l, (b - a) / h / (u - t)}' >> $O
              a=$b; t=$u; done ;;
    esac ;;
caller)
    caller "$@" ;;
*)
    echo "unknown experiment $X; see the header of $0" >&2; exit 1 ;;
esac
