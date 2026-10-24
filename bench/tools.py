#!/usr/bin/env python3
"""Small tools of the benchmark: measurement of runs, packed-read manipulation, VCF normalisation and scoring, and the
simulation truth. Usage: tools.py <command> [arguments]; see the functions at the end."""
import bisect, collections, datetime, glob, gzip, json, mmap, os, platform, re, struct, subprocess, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCH = f"{ROOT}/bench"
HAPLOID = 3117292070  # CHM13v2.0, the reference of the simulation
ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
TOL = 500  # breakpoint tolerance of every comparison
MIN = 50  # minimum SV size


# ---------------------------------------------------------------- measurement of a janus run (run.json, steps.tsv, files.tsv)
def measure(D):
    """Measurements of a janus run made by `janus run` under D/work with the monitors of bench/env.sh in D/monitor."""
    D = os.path.abspath(D); W = f"{D}/work"; M = f"{D}/monitor"
    for f in ("start.txt", "total.txt", "note.txt"):
        if os.path.exists(f"{D}/{f}"): os.replace(f"{D}/{f}", f"{M}/{f}")
    P = json.load(open(f"{M}/params.json")) if os.path.exists(f"{M}/params.json") else {}

    def iostat():
        out, t, first = [], None, True
        for l in open(f"{M}/iostat.txt"):
            f = l.split()
            if len(f) == 3 and f[2] in ("AM", "PM") and "/" in f[0]:
                t = datetime.datetime.strptime(" ".join(f), "%m/%d/%Y %I:%M:%S %p").timestamp()
            elif f and f[0].startswith("nvme") and t is not None:
                if first: first = False; continue
                out.append((t, float(f[2]), float(f[8]), float(f[-1])))
        return out

    def temps():
        out = []
        for l in open(f"{M}/temp.txt"):
            f = l.split()
            if len(f) < 3: continue
            out.append((int(f[0]), int(f[1]) / 1000, int(f[2]) / 1000, float(f[3].strip("+°C")) if len(f) > 3 else None))
        return out

    IO, TP = iostat(), temps()

    def window(a, b):
        io_ = [r for r in IO if a < r[0] <= b + 5]; tp = [r for r in TP if a <= r[0] <= b + 5]
        if not io_: rd = wr = ut = None
        else: rd = round(sum(r[1] for r in io_) * 5 / 1000, 2); wr = round(sum(r[2] for r in io_) * 5 / 1000, 2); ut = round(sum(r[3] for r in io_) / len(io_), 1)
        mx = lambda i: max((r[i] for r in tp if r[i] is not None), default=None)
        return rd, wr, ut, mx(1), mx(2), mx(3)

    def report(n):
        f = f"{W}/log.{n}"
        return [l.strip() for l in open(f, errors="replace") if l.startswith("[") and not l.startswith("[M::")] if os.path.exists(f) else []

    steps = []
    for l in open(f"{W}/steps.tsv").read().split("\n")[1:]:
        if not l: continue
        n, wall, user, sy, rss = l.split("\t")
        end = os.path.getmtime(f"{W}/done.{n}"); a = end - float(wall)
        s = dict(step=n, start=round(a, 2), end=round(end, 2), wall_s=float(wall), peak_rss_kib=int(rss), user_s=float(user), sys_s=float(sy))
        s["cpu_s"] = round(s["user_s"] + s["sys_s"], 2)
        s["read_gb"], s["written_gb"], s["util_pct"], s["drive_c"], s["controller_c"], s["cpu_c"] = window(a, end)
        s["report"] = report(n)
        steps.append(s)
    prev = steps[-1]["end"]
    if os.path.exists(f"{M}/total.txt"):
        start = float(open(f"{M}/start.txt").read()); tot = open(f"{M}/total.txt").read().split()
        total = dict(wall_s=float(tot[1]), peak_rss_kib=int(tot[4]), user_s=float(tot[7]), sys_s=float(tot[10]))
    else:
        start = steps[0]["start"]
        total = dict(wall_s=round(sum(s["wall_s"] for s in steps), 2), peak_rss_kib=max(s["peak_rss_kib"] for s in steps), user_s=round(sum(s["user_s"] for s in steps), 2), sys_s=round(sum(s["sys_s"] for s in steps), 2))
    total["cpu_s"] = round(total["user_s"] + total["sys_s"], 2)
    total["step_wall_s"] = round(sum(s["wall_s"] for s in steps), 2)
    total["peak_step"] = max(steps, key=lambda s: s["peak_rss_kib"])["step"]
    total["read_gb"], total["written_gb"], total["util_pct"], total["drive_c"], total["controller_c"], total["cpu_c"] = window(start, prev)
    files = sorted(((os.path.getsize(p), os.path.relpath(p, W)) for p in glob.glob(f"{W}/*") if os.path.isfile(p) and not os.path.islink(p)), reverse=True)
    KEEP = {"calls.vcf", "janus.vcf", "steps.tsv"}
    files = [dict(file=f, bytes=b, kept=f in KEEP or f.startswith("log.")) for b, f in files if not f.startswith("done.")]
    total["intermediate_gb"] = round(sum(f["bytes"] for f in files) / 1e9, 2)
    cols = ["step", "start", "end", "wall_s", "peak_rss_kib", "user_s", "sys_s", "cpu_s", "read_gb", "written_gb", "util_pct", "drive_c", "controller_c", "cpu_c"]
    with open(f"{D}/steps.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for s in steps + [dict(total, step="total", start=round(start, 2), end=round(prev, 2))]: o.write("\t".join("" if s.get(c) is None else str(s[c]) for c in cols) + "\n")
    with open(f"{D}/files.tsv", "w") as o:
        o.write("file\tbytes\tkept\n")
        for f in files: o.write(f"{f['file']}\t{f['bytes']}\t{'yes' if f['kept'] else 'deleted'}\n")
    reads = {}
    for name, path in P.get("samples", {}).items():
        reads[name] = dict(file=os.path.realpath(path), bytes=os.path.getsize(path))
        for st in steps:
            if st["step"] == f"pack_{name}":
                for l in st["report"]:
                    m = re.match(r"\[pack\] (\d+) reads, (\d+) bases", l)
                    if m: reads[name].update(reads=int(m.group(1)), bases=int(m.group(2)))
    cpu = next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")), "")
    mem = int(next(l.split()[1] for l in open("/proc/meminfo") if l.startswith("MemTotal")))
    commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    note = open(f"{M}/note.txt").read().strip() if os.path.exists(f"{M}/note.txt") else ""
    json.dump(dict(pair=P.get("pair", os.path.basename(D)), run=os.path.basename(D), params=P, commit=commit, note=note, started=datetime.datetime.fromtimestamp(start).isoformat(timespec="seconds"),
                   threads=P.get("threads", 16), cores="taskset 0-%d" % (P.get("threads", 16) - 1), machine=dict(cpu=cpu, mem_kib=mem, kernel=platform.release()), reads=reads, total=total, steps=steps, files=files),
              open(f"{D}/run.json", "w"), indent=1)
    print(f"{P.get('pair', D)}: {total['wall_s']} s wall, {total['cpu_s']} CPU s, peak {total['peak_rss_kib'] / 2**20:.1f} GiB ({total['peak_step']}), intermediates {total['intermediate_gb']} GB")


def tool(D):
    """Measurements of another caller's run (bench/env.sh measured()): time.txt and the cgroup of its container."""
    D = os.path.abspath(D); P = json.load(open(f"{D}/params.json"))
    wall, rss, us, sy, ex = open(f"{D}/time.txt").read().split()[-5:]
    t = dict(wall_s=float(wall), host_peak_kib=int(rss), host_user_s=float(us), host_sys_s=float(sy), exit=int(ex))
    cg = {}
    if os.path.exists(f"{D}/cgroup.txt"):
        for l in open(f"{D}/cgroup.txt"):
            k, v = l.split(); cg[k] = int(v)
    t["docker_user_s"] = cg.get("user_usec", 0) / 1e6; t["docker_sys_s"] = cg.get("system_usec", 0) / 1e6; t["docker_peak_kib"] = cg.get("memory_peak", 0) // 1024
    t["user_s"] = round(t["host_user_s"] + t["docker_user_s"], 2); t["sys_s"] = round(t["host_sys_s"] + t["docker_sys_s"], 2); t["cpu_s"] = round(t["user_s"] + t["sys_s"], 2)
    t["docker_anon_peak_kib"] = cg.get("anon_peak", 0) // 1024
    t["peak_rss_kib"] = max(t["host_peak_kib"], t["docker_anon_peak_kib"] if "anon_peak" in cg else t["docker_peak_kib"])
    O = f"{D}/out"; files = []
    for root, _, fs in os.walk(O):
        for f in fs:
            p = os.path.join(root, f)
            if os.path.isfile(p) and not os.path.islink(p): files.append(dict(file=os.path.relpath(p, O), bytes=os.path.getsize(p), kept=os.path.relpath(p, O) in ("calls.vcf", "log.txt")))
    for k in ("bam", "bam_n", "bam_t"):
        if k in P and os.path.exists(P[k]): files.append(dict(file=os.path.basename(P[k]), bytes=os.path.getsize(P[k]), kept=False))
    files.sort(key=lambda f: -f["bytes"])
    calls = sum(1 for l in open(f"{O}/calls.vcf") if not l.startswith("#")) if os.path.exists(f"{O}/calls.vcf") else None
    json.dump(dict(params=P, total=t, calls=calls, files=files), open(f"{D}/run.json", "w"), indent=1)
    with open(f"{D}/files.tsv", "w") as o:
        o.write("file\tbytes\tkept\n")
        for f in files: o.write(f"{f['file']}\t{f['bytes']}\t{'yes' if f['kept'] else 'deleted'}\n")
    print(f"{P.get('tool')} {P.get('pair')} {P.get('coverage_x') or 'full'}: exit {t['exit']}, {t['wall_s']:.0f} s wall, {t['cpu_s']:.0f} CPU s, peak {t['peak_rss_kib'] / 2**20:.1f} GiB, calls {calls}")


# ---------------------------------------------------------------- packed reads (SVB1 files of janus pack)
def _svb_write(out, parts):
    """Write the records of `parts` (iterables of (bytes, length) records) as out and out.idx."""
    n = total = longest = 0; off = []
    with open(out, "wb") as o:
        o.write(b"\0" * 32); pos = 32
        for rec, ln in parts:
            o.write(rec); off.append(pos); pos += len(rec); total += ln; n += 1; longest = max(longest, ln)
        off.append(pos)
        o.seek(0); o.write(b"SVB1" + struct.pack("<IQQ", 1, n, total) + b"\0" * 8)
    with open(out + ".idx", "wb") as o:
        o.write(b"SVBI" + struct.pack("<IQQ", 1, n, longest) + b"\0" * 8)
        o.write(struct.pack("<%dQ" % len(off), *off))
    print(out, n, "reads", total, "bases", "%.2fx" % (total / HAPLOID))


def _svb_records(src, skip=0, take=None):
    """Records of a packed file after skipping `skip` bases, up to `take` bases (the whole file by default)."""
    with open(src, "rb") as f:
        m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        p, seen, got, end = 32, 0, 0, len(m)
        while p < end and (take is None or got < take):
            tag = struct.unpack_from("<I", m, p)[0]; ln = tag & 0x7fffffff
            sz = 4 + (ln if tag >> 31 else (ln + 3) // 4)
            if seen >= skip: yield m[p:p + sz], ln; got += ln
            seen += ln; p += sz
        if take is not None and got < take: sys.exit(f"{src}: only {got} of {take} bases after skipping {skip}")


def mix(out, *parts):
    """mix OUT SRC:SKIP:TAKE ... : the first TAKE bases of SRC after skipping SKIP bases, from each part in turn."""
    def gen():
        for a in parts:
            src, skip, take = a.split(":"); yield from _svb_records(src, int(skip), int(take))
    _svb_write(out, gen())


def thin(src, out, stride):
    """Every STRIDE-th read of SRC."""
    _svb_write(out, (r for i, r in enumerate(_svb_records(src)) if i % int(stride) == 0))


def svbcat(out, *ins):
    """The reads of several packed files, one after another."""
    def gen():
        for f in ins: yield from _svb_records(f)
    _svb_write(out, gen())


# ---------------------------------------------------------------- VCF records and comparisons
def recs(path):
    """PASS records of a VCF as (type, chrom, pos, chrom2, pos2), one per breakend pair."""
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
            mates.add(f[2]); out.append((t, c, p, m.group(1), int(m.group(2))))
        else: out.append((t, c, p, c, int(float(info.get("END", p)))))
    return out


def hit(a, b):
    for x in ((b[1], b[2], b[3], b[4]), (b[3], b[4], b[1], b[2])):
        if a[1] == x[0] and a[3] == x[2] and abs(a[2] - x[1]) <= TOL and abs(a[4] - x[3]) <= TOL: return True
    if a[0] == "INS" or b[0] == "INS": return a[1] == b[1] and abs(a[2] - b[2]) <= TOL
    return False


def score(truth, calls):
    """Recall, precision and F1 of a call set against a truth set, breakpoints within TOL."""
    T, C = recs(truth), recs(calls)
    tp_t = sum(any(hit(t, c) for c in C) for t in T); tp_c = sum(any(hit(c, t) for t in T) for c in C)
    r, p = tp_t / len(T), (tp_c / len(C) if C else 0)
    print("truth %d calls %d recall %.3f precision %.3f f1 %.3f" % (len(T), len(C), r, p, 2 * r * p / (r + p) if r + p else 0))


def norm(path, name):
    """A caller's PASS records in the form Minda reads: symbolic DEL/DUP/INV/INS, everything else as breakend pairs."""
    print("##fileformat=VCFv4.2"); print(f"##source={name}")
    print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">'); print('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">')
    print('##INFO=<ID=END,Number=1,Type=Integer,Description="end">'); print('##INFO=<ID=MATEID,Number=1,Type=String,Description="mate">')
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    seen, n = set(), 0
    for l in open(path):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t")
        if f[6] not in ("PASS", "."): continue
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        t, c, p = info.get("SVTYPE", "BND"), f[0], int(f[1]); m = ALT.search(f[4]); n += 1
        if m or info.get("CHR2", c) != c:
            c2, p2 = (m.group(1), int(m.group(2))) if m else (info["CHR2"], int(info.get("CHR2_POS", info.get("END", p))))
            k = tuple(sorted([(c, p), (c2, p2)]))
            if k in seen: continue
            seen.add(k)
            print(f"{c}\t{p}\t{name}_{n}_1\tN\tN[{c2}:{p2}[\t.\tPASS\tSVTYPE=BND;MATEID={name}_{n}_2"); print(f"{c2}\t{p2}\t{name}_{n}_2\tN\tN[{c}:{p}[\t.\tPASS\tSVTYPE=BND;MATEID={name}_{n}_1")
            continue
        e = int(info.get("END", p))
        try: ln = abs(int(info["SVLEN"].split(",")[0]))
        except (KeyError, ValueError): ln = abs(e - p)
        if t == "INS": print(f"{c}\t{p}\t{name}_{n}\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN={ln};END={p}")
        else:
            t = t if t in ("DEL", "DUP", "INV") else "BND"
            if t == "BND": print(f"{c}\t{p}\t{name}_{n}_1\tN\tN[{c}:{e}[\t.\tPASS\tSVTYPE=BND;MATEID={name}_{n}_2"); print(f"{c}\t{e}\t{name}_{n}_2\tN\tN[{c}:{p}[\t.\tPASS\tSVTYPE=BND;MATEID={name}_{n}_1")
            else: print(f"{c}\t{p}\t{name}_{n}\tN\t<{t}>\t.\tPASS\tSVTYPE={t};SVLEN={ln};END={e}")


def subtract(tumor, normal, dist=1000, min_len=50, somatic_tag=False):
    """Tumor-only PASS SVs: records of the tumor VCF without a record of the same type within `dist` in the normal VCF."""
    dist, min_len = int(dist), int(min_len)

    def info(f):
        return dict(kv.partition("=")[::2] for kv in f[7].split(";"))

    def svtype(f, d): return d["SVTYPE"] if "SVTYPE" in d else ("INS" if len(f[4]) > len(f[3]) else "DEL")

    def is_sv(f, d):
        if f[4].startswith("<") or "[" in f[4] or "]" in f[4]: return True
        if "SVLEN" in d:
            try: return abs(int(d["SVLEN"].split(",")[0])) >= min_len
            except ValueError: pass
        return abs(len(f[4]) - len(f[3])) >= min_len

    nm = collections.defaultdict(list)
    for l in open(normal):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); d = info(f)
        if is_sv(f, d): nm[(f[0], svtype(f, d))].append(int(f[1]))
    for l in open(tumor):
        if l[0] == "#": sys.stdout.write(l); continue
        f = l.rstrip("\n").split("\t"); d = info(f)
        if f[6] not in ("PASS", ".") or not is_sv(f, d): continue
        if somatic_tag and "SOMATIC" not in d: continue
        if any(abs(p - int(f[1])) <= dist for p in nm[(f[0], svtype(f, d))]): continue
        sys.stdout.write(l)


def savana(*args):
    """Run SAVANA with its worker processes reaped on timeout (its multiprocessing join hangs otherwise)."""
    import multiprocessing.process as mp
    from savana.savana import main
    join = mp.BaseProcess.join

    def reap(self, timeout=None):
        join(self, timeout)
        if timeout is not None and self.is_alive(): self.kill(); join(self)
    mp.BaseProcess.join = reap
    sys.argv = ["savana"] + list(args)
    sys.exit(main())


# ---------------------------------------------------------------- janus records: the joint-run grouping used by the simulation scoring
NEAR = 500


def janus_records(path):
    """Records of a janus VCF: keys (chrom, chrom2, pos, pos2), FILTER, the tumor index of the record and its CARRIERS."""
    K, F, own, carriers, names = [], [], [], [], []
    for l in open(path):
        if l[0] == "#":
            if l.startswith("#CHROM"): names = l.rstrip("\n").split("\t")[10:]
            continue
        f = l.rstrip("\n").split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
        K.append((f[0], m.group(1) if m else f[0], int(f[1]), int(m.group(2)) if m else int(f[1]) + int(info["SVLEN"])))
        F.append(f[6]); gt = [x.split(":")[0] for x in f[10:]]; own.append(gt.index("0/1") if "0/1" in gt else 0)
        carriers.append([names.index(c) for c in info["CARRIERS"].split(",")] if "CARRIERS" in info else [])
    return K, F, own, carriers


def groups(K, skip):
    """The SV groups of janus records: records of different samples within NEAR of each other (transitive closure), skipping `skip` records."""
    key = lambda q, da=0, db=0: (q[0], q[1], q[2] // NEAR + da, q[3] // NEAR + db)
    grid = collections.defaultdict(list)
    for i, k in enumerate(K):
        if not skip[i]: grid[key(k)].append(i)
    root = list(range(len(K)))

    def find(i):
        while root[i] != i: root[i] = root[root[i]]; i = root[i]
        return i
    for i, k in enumerate(K):
        if skip[i]: continue
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for j in grid.get(key(k, da, db), ()):
                    if abs(K[j][2] - k[2]) <= NEAR and abs(K[j][3] - k[3]) <= NEAR: root[find(i)] = find(j)
    return [find(i) for i in range(len(K))]


# ---------------------------------------------------------------- the simulation: truth junctions, scoring, lost copies
def simtruth(G, OUT):
    """Lift every simulated SV of the tree to reference junctions: OUT/<tumor>.tsv (one line per junction) and OUT/tree.json."""
    plan = json.load(open(f"{G}/tumor/plan.json")); PAR = plan["parent"]
    base = lambda c: "_".join(c.split("_")[:2])

    def load(path, segs):
        for l in open(path):
            f = l.rstrip("\n").split("\t"); m = re.search(r"cg:Z:(\d+)M$", f[12]) if len(f) > 12 else None
            if not m: continue
            segs.setdefault(base(f[0]), []).append((int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8])))
    PAF = {}

    def paf(node):
        if node not in PAF:
            segs = {}
            if node == "germline":
                for c in os.listdir(f"{G}/germline"):
                    for h in ("A", "B"):
                        p = f"{G}/germline/{c}/sim.hap{h}.paf"
                        if os.path.exists(p): load(p, segs)
            else: load(f"{G}/tumor/{node}/sim.paf", segs)
            for k in segs: segs[k].sort()
            PAF[node] = {k: (v, [x[0] for x in v]) for k, v in segs.items()}
        return PAF[node]

    def step(node, c, p):
        d = paf(node).get(base(c))
        if not d: return None
        v, starts = d; i = bisect.bisect_right(starts, p) - 1
        if i < 0: i, p = 0, v[0][0]
        qs, qe, st, t, ts, te = v[i]; snap = False
        if p >= qe:
            if i + 1 < len(v) and v[i + 1][0] - p < p - (qe - 1): qs, qe, st, t, ts, te = v[i + 1]; p = qs
            else: p = qe - 1
            snap = True
        return (t, ts + (p - qs) if st == "+" else te - 1 - (p - qs), snap)

    def lift(node, c, p):
        snapped = False
        while True:
            r = step(node, c, p)
            if r is None: return None
            c, p, s = r; snapped |= s
            if node == "germline": return c, p, snapped
            node = PAR[node] or "germline"

    def chain(n):
        out = []
        while n: out.append(n); n = PAR[n]
        return out
    rows = [l.rstrip("\n").split("\t") for l in open(f"{G}/tumor/mutations.tsv")]
    rows = [r for r in rows[1:] if r[5] == "SV"]
    lifted, failed, snaps = {}, 0, 0
    for r in rows:
        mid, birth, c, pos, end, svt, svlen = r[0], r[1], r[2], int(r[3]), int(r[4]), r[6], int(r[7])
        info = dict(x.split("=", 1) for x in r[10].split(";") if "=" in x); start = PAR[birth] or "germline"
        a, b = lift(start, c, pos), lift(start, c, end)
        if not a or not b: failed += 1; continue
        if svt in ("dDUP", "rTRA"):
            tg = lift(start, info["TARGET_CHROM"], int(info["TARGET"]))
            if not tg: failed += 1; continue
            J = [(tg[0], tg[1], a[0], a[1]), (b[0], b[1], tg[0], tg[1])]; snaps += tg[2]
        elif svt == "INS": J = [(a[0], a[1], a[0], a[1])]
        else: J = [(a[0], a[1], b[0], b[1])]
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
            full = [x for x in L if x[2] > 0]; empty = [x for x in L if x[2] == 0]; J = []
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
                    for k, j in enumerate(d["j"]): o.write(f"{mid}\t{d['type']}\t{d['len']}\t{k}\t{j[0]}\t{j[1]}\t{j[2]}\t{j[3]}\n")
    json.dump(dict(parent=PAR, leaves=leaves), open(f"{OUT}/tree.json", "w"))
    print(f"SVs {len(rows)}, lifted {len(lifted)} ({len(lifted) / len(rows):.3%}), failed {failed}, snapped {snaps}")
    print(" ".join(f"{t} {len(set(l.split(chr(9))[0] for l in open(f'{OUT}/{t}.tsv')))}" for t in leaves))


def truth(path):
    """The simulated SVs of a tumor: {id: type} and the junction list (chrom, pos, chrom2, pos2, id)."""
    sv, J = {}, []
    for l in open(path):
        mid, typ, ln, k, c1, p1, c2, p2 = l.rstrip("\n").split("\t")
        if typ not in ("dDUP", "rTRA") and int(ln) < MIN: continue
        sv[mid] = typ; J.append((c1, int(p1), c2, int(p2), mid))
    return sv, J


def calls(path, passonly=True):
    """The junctions of a call set, one list per record."""
    out = []
    for l in (open(path) if isinstance(path, str) else path):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t")
        if passonly and f[6] not in ("PASS", "."): continue
        info = dict(x.split("=", 1) if "=" in x else (x, "") for x in f[7].split(";")); t = info.get("SVTYPE", ""); m = ALT.search(f[4])
        if m or t in ("BND", "TRA"):
            if m: out.append([(f[0], int(f[1]), m.group(1), int(m.group(2)))])
            elif "CHR2" in info: out.append([(f[0], int(f[1]), info["CHR2"], int(info.get("END", info.get("POS2", f[1]))))])
            continue
        ln = abs(int(info.get("SVLEN", "0") or 0)); end = int(info.get("END", f[1]))
        if t == "INS":
            if ln and ln < MIN: continue
            out.append([(f[0], int(f[1]), f[0], int(f[1]))])
        else:
            if max(ln, end - int(f[1])) < MIN: continue
            out.append([(f[0], int(f[1]), f[0], end)])
    return out


def junction(line):
    """The junction of one janus record regardless of its FILTER (None when it is too small)."""
    g = line.rstrip("\n").split("\t"); c = calls(["\t".join(g[:6] + ["PASS"] + g[7:]) + "\n"])
    return c[0] if c else None


def index(J):
    by = {}
    for j in J: by.setdefault(j[0], []).append((j[1], j)); by.setdefault(j[2], []).append((j[3], j))
    return {c: (sorted(v), [x[0] for x in sorted(v)]) for c, v in by.items()}


def near(a, b):
    return (a[0] == b[0] and abs(a[1] - b[1]) <= TOL and a[2] == b[2] and abs(a[3] - b[3]) <= TOL) or (a[0] == b[2] and abs(a[1] - b[3]) <= TOL and a[2] == b[0] and abs(a[3] - b[1]) <= TOL)


def match(c, I):
    found = set()
    for (c1, p1, c2, p2) in c:
        d = I.get(c1)
        if not d: continue
        v, ps = d
        for i in range(bisect.bisect_left(ps, p1 - TOL), bisect.bisect_right(ps, p1 + TOL)):
            j = v[i][1]
            if near((c1, p1, c2, p2), j): found.add(j[4])
    return found


def simscore(tpath, cpath, passonly=True):
    """Precision, recall and F1 of a call set against a tumor's simulated SVs."""
    sv, J = truth(tpath); I = index(J); C = calls(cpath, passonly); found, tp = set(), 0
    for c in C:
        h = match(c, I)
        if h: tp += 1; found |= h
    P = tp / len(C) if C else 0; R = len(found) / len(sv) if sv else 0
    return dict(truth=len(sv), calls=len(C), tp_calls=tp, found=len(found), P=P, R=R, F1=2 * P * R / (P + R) if P + R else 0, found_ids=found)


def lost(G, OUT, fai):
    """OUT/<tumor>.lost: the simulated SVs of each tumor that lie on a chromosome copy the tumor lost (copy number 0 of the haplotype)."""
    FAI = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
    hap = {f[0]: f[2].split("_")[1] for f in (l.split("\t") for l in open(f"{G}/tumor/mutations.tsv")) if f[5] == "SV"}
    for t in sorted(json.load(open(f"{OUT}/tree.json"))["leaves"]):
        seg = collections.defaultdict(list)
        for l in open(f"{G}/tumor/{t}/cn_segments.bed"):
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t"); seg[f[3]].append((float(f[4]), float(f[5]), int(f[9]), int(f[10])))
        sv, J = truth(f"{OUT}/{t}.tsv"); first = {}
        for j in J: first.setdefault(j[4], j)
        gone = set()
        for m in sv:
            j = first[m]; c = j[2] if sv[m] in ("dDUP", "rTRA") else j[0]; p = j[3] if sv[m] in ("dDUP", "rTRA") else j[1]; fr = p / FAI[c]; cn = None
            for s0, s1, ca, cb in seg.get(c, []):
                if s0 <= fr < s1: cn = ca if hap[m] == "hapA" else cb; break
            if cn == 0: gone.add(m)
        open(f"{OUT}/{t}.lost", "w").write("\n".join(sorted(gone)) + "\n")
        print(t, "truth", len(sv), "on a lost copy", len(gone), "%.1f%%" % (100 * len(gone) / len(sv)))


# ---------------------------------------------------------------- read lists of the subsampled BAMs for the other callers (coverage and purity)
G_HG38 = 3100000000


def select(p, O, castle):
    """Read names of the depth and purity subsets of pair p (BAM read names in O/names/<p>.<s>.tsv, lists into O/lists), with the
    same shuffled order and base budgets as the packed-read mixes of janus."""
    import numpy as np
    M = (1 << 64) - 1

    def load(s):
        names, lens = [], []
        for l in open(f"{O}/names/{p}.{s}.tsv"):
            a, b = l.rstrip("\n").split("\t"); names.append(a); lens.append(int(b))
        n = len(names); order = list(range(n)); seed = 0x9e3779b97f4a7c15
        for i in range(n, 1, -1):
            seed ^= (seed << 13) & M; seed ^= seed >> 7; seed ^= (seed << 17) & M; j = seed % i
            order[i - 1], order[j] = order[j], order[i - 1]
        order = np.array(order); ln = np.array(lens, dtype=np.int64)[order]
        return names, order, ln, np.cumsum(ln)

    def take(cum, ln, skip, want):
        before = cum - ln; a = int(np.searchsorted(before, skip, side="left")); base = int(before[a]) if a < len(cum) else int(cum[-1])
        b = int(np.searchsorted(cum, base + want, side="left"))
        return a, b + 1, int(cum[b] - base)
    os.makedirs(f"{O}/lists", exist_ok=True); D = {s: load(s) for s in "nt"}

    def put(path, parts):
        with open(path, "w") as o:
            for s, a, b in parts:
                names, order, ln, cum = D[s]; o.write("\n".join(names[i] for i in order[a:b]) + "\n")
    for x in (5, 10, 15, 20, 30):
        for s in "nt":
            names, order, ln, cum = D[s]
            if cum[-1] < x * G_HG38: continue
            a, b, tot = take(cum, ln, 0, x * G_HG38); put(f"{O}/lists/{p}.{s}.{x}x.txt", [(s, a, b)])
    if p in castle.split(","):
        nb = int(D["n"][3][-1]); nm = min(nb - 24 * G_HG38 - G_HG38 // 10, 30 * G_HG38)
        a, b, tot = take(D["n"][3], D["n"][2], 0, nm); put(f"{O}/lists/{p}.n.purity.txt", [("n", a, b)])
        for q in (20, 40, 60, 80, 100):
            ta, tb, tt = take(D["t"][3], D["t"][2], 0, q * 30 * G_HG38 // 100); parts = [("t", ta, tb)]
            if q < 100:
                na, nb2, nt = take(D["n"][3], D["n"][2], nm, (100 - q) * 30 * G_HG38 // 100); parts.append(("n", na, nb2))
            put(f"{O}/lists/{p}.t.{q}p.tumor_part.txt", parts[:1])
            if q < 100: put(f"{O}/lists/{p}.t.{q}p.normal_part.txt", parts[1:])
    print(p, "read lists written")


# ---------------------------------------------------------------- scoring of the simulation runs (simulation.sh score)
SIMSETS = [("20-20", "20_20"), ("20-10", "20_10"), ("20-5", "20_5"), ("10-10", "10_10"), ("5-5", "5_5")]
T20 = ["T%02d" % i for i in range(1, 21)]


def _minimum(vcf):
    """The minimum read support of a run, from the smallest tumor read count of its records (log2 count + 1 is feature 1)."""
    lo = min(float(l.split("\t")[7].split("F=")[1].split(";")[0].split(",")[1]) for l in open(vcf) if l[0] != "#")
    return int(round(2 ** lo - 1))


def _junctions(vcf):
    return [junction(l) for l in open(vcf) if l[0] != "#"]


def _tumor_lists(vcf, T):
    """Per tumor: the indices of its PASS records followed by the records it carries without a record of its own."""
    K, F, own, carriers = janus_records(vcf); lists, adds = [], []
    for k in range(T):
        mine = [i for i in range(len(K)) if F[i] == "PASS" and own[i] == k]; add = [i for i in range(len(K)) if k in carriers[i]]
        lists.append(mine + add); adds.append(add)
    return lists, adds


def simjoint(R=f"{BENCH}/runs/sim"):
    """runs/sim/scores.json: per depth setting, the precision and recall of every tumor from the pairwise runs, the first
    pass and the clade pass of the joint run, the clade pass at every penalty mu, and (20/20, 20/5) the first pass under
    every acceptance rule of the supplement; from the classified VCFs written by simulation.sh score."""
    TI, NT = [], []
    for t in T20: sv, J = truth(f"{R}/truth/{t}.tsv"); TI.append(index(J)); NT.append(len(sv))
    out = {}
    for name, s in SIMSETS:
        d = f"{R}/settings/{s}"
        if not all(os.path.exists(f"{d}/{t}/janus.vcf") for t in T20) or not os.path.exists(f"{d}/joint/janus.vcf"): print(name, "not ready"); continue
        pair = []
        for k, t in enumerate(T20):
            v = f"{d}/{t}/janus.vcf"; J = _junctions(v); F = [l.split("\t")[6] for l in open(v) if l[0] != "#"]
            pp = [i for i in range(len(J)) if F[i] == "PASS"]; hits = [match(J[i], TI[k]) if J[i] else set() for i in pp]; f = set().union(*hits) if hits else set()
            pair.append((sum(1 for h in hits if h) / max(len(pp), 1), len(f) / NT[k], len(pp), _minimum(f"{d}/{t}/calls.vcf")))
        JJ = _junctions(f"{d}/joint/calls.vcf")

        def evaluate(vcf):
            lists, adds = _tumor_lists(vcf, 20); A = []
            for k in range(20):
                hits = [match(JJ[i], TI[k]) if JJ[i] else set() for i in lists[k]]; f = set().union(*hits) if hits else set()
                A.append((sum(1 for h in hits if h) / max(len(hits), 1), len(f) / NT[k], len(hits), len(adds[k]), sum(1 for h in hits[len(hits) - len(adds[k]):] if h)))
            return A
        res = lambda A: dict(P=[x[0] for x in A], R=[x[1] for x in A], calls=[x[2] for x in A], carried=[x[3] for x in A], carried_true=[x[4] for x in A])
        med = lambda A: [float(np.median([x[0] for x in A])), float(np.median([x[1] for x in A])), float(np.median([2 * x[0] * x[1] / max(x[0] + x[1], 1e-9) for x in A]))]
        A1 = evaluate(f"{d}/joint/first.vcf"); A = evaluate(f"{d}/joint/janus.vcf")
        out[name] = dict(truth=NT, pair=dict(P=[x[0] for x in pair], R=[x[1] for x in pair], calls=[x[2] for x in pair], minimum=[int(x[3]) for x in pair]), first=res(A1), joint=res(A))
        MU = [0.3, 0.5, 1.0, 1.5, 2.0]; out[name]["mu"] = dict(mu=MU, PRF=[med(evaluate(f"{d}/joint/mu_{m}.vcf")) for m in MU])
        if name in ("20-20", "20-5"):
            RA, FL = [0.0, 0.1, 0.2, 0.3], [0.05, 0.1, 0.2, 0.3]
            out[name]["sweep"] = dict(raise_=RA, floor=FL, PRF=[[med(evaluate(f"{d}/joint/rule_{r}_{f}.vcf")) for f in FL] for r in RA])
        P1, R1, P2, R2 = [np.array(x) for x in (out[name]["pair"]["P"], out[name]["pair"]["R"], out[name]["joint"]["P"], out[name]["joint"]["R"])]
        F1, F2 = 2 * P1 * R1 / np.maximum(P1 + R1, 1e-9), 2 * P2 * R2 / np.maximum(P2 + R2, 1e-9)
        print(f"{name}  one at a time P {np.median(P1):.3f} R {np.median(R1):.3f} F1 {np.median(F1):.3f}  first pass P {med(A1)[0]:.3f} R {med(A1)[1]:.3f} F1 {med(A1)[2]:.3f}"
              f"  joint P {np.median(P2):.3f} R {np.median(R2):.3f} F1 {np.median(F2):.3f}  both better {int(((P2 > P1) & (R2 > R1)).sum())}/20  F1 better {int((F2 > F1).sum())}/20"
              f"  dF1 {np.min(F2 - F1):+.3f}..{np.max(F2 - F1):+.3f}  carriers {sum(x[3] for x in A1)} ({sum(x[4] for x in A1)} true) -> {sum(x[3] for x in A)} ({sum(x[4] for x in A)} true)", flush=True)
    json.dump(out, open(f"{R}/scores.json", "w"), indent=1)


def pooledscore(s="20_5", R=f"{BENCH}/runs/sim"):
    """runs/sim/pooled/scores.json: the 20 tumors pooled into one sample, the 20 pairwise runs and the joint run, scored against
    the union of the tumors' SVs, with the recall by the number of tumors carrying the SV."""
    sv, J, NC = {}, {}, {}
    for t in T20:
        a, b = truth(f"{R}/truth/{t}.tsv"); sv.update(a)
        for mid in a: NC[mid] = NC.get(mid, 0) + 1
        for j in b: J[j] = j
    I = index(list(J.values())); CLS = [("1", 1, 1), ("2-5", 2, 5), ("6-19", 6, 19), ("20", 20, 20)]

    def score_(calls):
        found, tp = set(), 0
        for c in calls:
            h = match(c, I)
            if h: tp += 1; found |= h
        P = tp / max(len(calls), 1); Rr = len(found) / len(sv)
        cls = {lab: (lambda ids: [len(ids), len(ids & found) / max(len(ids), 1)])({m for m, n in NC.items() if lo <= n <= hi}) for lab, lo, hi in CLS}
        return dict(calls=len(calls), tp=tp, found=len(found), P=P, R=Rr, F1=2 * P * Rr / max(P + Rr, 1e-9), classes=cls)

    def passing(vcf):
        K, F, own, carriers = janus_records(vcf); C = _junctions(vcf)
        return [(k, c) for k, c, f in zip(K, C, F) if f == "PASS" and c]

    def merged(KC):
        g = groups([k for k, _ in KC], [False] * len(KC)); first = {}
        for i, r in enumerate(g): first.setdefault(r, i)
        return [KC[i][1] for i in first.values()]
    out = {"truth": len(sv), "pooled": score_(merged(passing(f"{R}/pooled/{s}/janus.vcf")))}
    pairs = []
    for t in T20: pairs += passing(f"{R}/settings/{s}/{t}/janus.vcf")
    out["pairwise"] = score_(merged(pairs)); out["joint"] = score_(merged(passing(f"{R}/settings/{s}/joint/janus.vcf")))
    out["time"] = [float(x) for x in open(f"{R}/pooled/{s}/time.txt").read().split()[:2]]
    A = json.load(open(f"{R}/pooled/scores.json")) if os.path.exists(f"{R}/pooled/scores.json") else {}
    A[s.replace("_", "-")] = out; json.dump(A, open(f"{R}/pooled/scores.json", "w"), indent=1)
    for n in ("pooled", "pairwise", "joint"):
        o = out[n]; print(f"{s} {n:9s} calls {o['calls']} tp {o['tp']} found {o['found']} P {o['P']:.3f} R {o['R']:.3f} F1 {o['F1']:.3f}  " + " ".join(f"{k}:{v[1]:.3f}" for k, v in o["classes"].items()))


# ---------------------------------------------------------------- the count-stage threshold sweep (sweep.sh)
def sweep(mode, *args):
    import hashlib, itertools
    R_ = f"{BENCH}/runs"; W = f"{R_}/sweep"; CASTLE = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]
    N0 = [5, 10, 15, 20, 25]; NT = [20, 30, 40, 50, 60]; GRID = list(itertools.product(N0, NT)) + [(70, 40), (10, 5), (70, 5)]

    def exps(p):
        e = [("full", p)] + [("coverage", f"{p}_{n}x") for n in (5, 10, 15, 20, 30)] + [("purity", f"{p}_{q}p") for q in (20, 40, 60, 80, 100)]
        return [(k, n) for k, n in e if os.path.exists(f"{R_}/{k}/{n}/t.cnt.hist")]

    def hist_(d): return {l.split()[0]: np.array(l.split()[1:], float) for l in open(f"{d}/t.cnt.hist")}

    def median(h):
        n = h["shared_normal"] + h["normal_only"]; tot = n[3:].sum(); run = 0; m = 3
        while m + 1 < len(n) and 2 * (run + n[m]) < tot: run += n[m]; m += 1
        return m

    def p95(h):
        t = h["tumor_only"]; tot = t[3:].sum(); run = 0; p = 3
        while p + 1 < len(t) and 20 * (run + t[p]) < 19 * tot: run += t[p]; p += 1
        return p

    def pair(d): return open(f"{d}/n.cnt", "rb").read(27)[26], open(f"{d}/t.cnt", "rb").read(27)[26]

    def runs(k, n):
        out = {}
        for d in [f"{R_}/{k}/{n}"] + sorted(glob.glob(f"{W}/runs/{k}/{n}/*/")):
            d = d.rstrip("/")
            if all(os.path.exists(f"{d}/{f}") for f in ("calls.vcf", "n.cnt", "t.cnt", "t.cnt.hist")): out.setdefault(pair(d), d)
        return out

    def choice(k, n, a, b):
        rs = runs(k, n); hs = {}
        for (m0, mt), d in rs.items(): hs.setdefault(m0, hist_(d))
        m0 = 2 if median(next(iter(hs.values()))) < a else 3
        if m0 not in hs: return m0, 0, None
        mt = 3 if p95(hs[m0]) >= b else 2
        return m0, mt, rs.get((m0, mt))
    if mode == "plan":
        for k, n in exps(args[0]):
            want = {}
            for a, b in GRID:
                m0, mt, d = choice(k, n, a, b)
                if not d: want.setdefault((m0, mt), (a, b))
            for (m0, mt), (a, b) in sorted(want.items()): print(k, n, a, b, m0, mt or "?")
    elif mode == "check":
        k, n, d = args[:3]; m0, mt = pair(d)
        ref = [x for (a, b), x in runs(k, n).items() if a == m0 and x != d]
        same = [all((hist_(d)[g] == hist_(x)[g]).all() for g in hist_(d)) for x in ref]
        print(f"{k}/{n} {os.path.basename(d)}: minimum read support {m0}/{mt}" + (f", count spectra {'identical to' if all(same) else 'DIFFERENT FROM'} the earlier run with m0 {m0}" if ref else ""))
    elif mode == "configs":
        os.makedirs(W, exist_ok=True); keys = {}; P = open(f"{W}/pairs.tsv", "w"); P.write("theta0\tthetat\texperiment\tm0\tmt\trun\n")
        for a, b in GRID:
            sel = {}
            for p in CASTLE:
                for k, n in exps(p):
                    m0, mt, d = choice(k, n, a, b); sel[f"{k}/{n}"] = d; P.write(f"{a}\t{b}\t{k}/{n}\t{m0}\t{mt}\t{d}\n")
            if not all(sel.values()): print(f"tau {a}/{b}: missing runs for {[e for e, d in sel.items() if not d]}", file=sys.stderr); continue
            key = hashlib.sha1("".join(f"{e}={d};" for e, d in sorted(sel.items())).encode()).hexdigest()[:10]; keys[(a, b)] = key
            for e, d in sel.items():
                l = f"{W}/conf/{key}/{e}"; os.makedirs(os.path.dirname(l), exist_ok=True)
                if not os.path.islink(l): os.symlink(d, l)
        P.close()
        with open(f"{W}/combos.tsv", "w") as o:
            o.write("theta0\ttheta_t\tconfig\n"); o.write("".join(f"{a}\t{b}\t{k}\n" for (a, b), k in keys.items()))
        print("\n".join(sorted(set(keys.values()))))
    elif mode == "rescore":
        C = args[0]; S = {}
        WITH = ["Severus", "SAVANA", "nanomonsv", "Sniffles2", "dysgu", "SVision-pro", "cuteSV", "DeBreak", "kled", "SVDSS", "janus", "Severus ONT", "SAVANA ONT", "nanomonsv ONT", "Sniffles2 ONT", "SvABA", "GRIPSS", "Manta"]
        LOO = {"janus": "janus", "severus_phased": "Severus", "savana_phased": "SAVANA", "nanomonsv": "nanomonsv", "dysgu": "dysgu", "svision-pro": "SVision-pro", "sniffles2": "Sniffles2", "cutesv": "cuteSV", "debreak": "DeBreak", "kled": "kled", "svdss": "SVDSS"}
        TOOL = {"janus": "janus", "severus": "Severus", "savana": "SAVANA", "nanomonsv": "nanomonsv"}
        SETS = ["with", "without", "loo", "30x", "20x", "15x", "10x", "5x", "100p", "80p", "60p", "40p", "20p"]
        rows = lambda d: [l.split("\t") for l in open(f"{d}/results/overall.tsv").read().split("\n")[1:] if l]
        for p in CASTLE:
            for s, names in (("with", WITH), ("without", [n for n in WITH if n != "janus"])):
                r = rows(f"{C}/eval_sets/{p}/all_{s}"); assert len(r) == len(names), (p, s)
                for n, x in zip(names, r): S[(s, n, p)] = [float(v) for v in x[4:7]]
            S[("without", "janus", p)] = [float(v) for v in rows(f"{C}/eval_sets/{p}/all_without_janus")[0][4:7]]
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


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    f = globals()[cmd]
    if cmd == "subtract":
        kw = {}
        for i, a in enumerate(args):
            if a == "--dist": kw["dist"] = args[i + 1]
            if a == "--min-len": kw["min_len"] = args[i + 1]
            if a == "--somatic-tag": kw["somatic_tag"] = True
        f(args[0], args[1], **kw)
    elif cmd == "simscore":
        r = simscore(args[0], args[1], len(args) < 3)
        print(f"truth {r['truth']} calls {r['calls']} recall {r['R']:.3f} precision {r['P']:.3f} f1 {r['F1']:.3f}")
    else: f(*args)
