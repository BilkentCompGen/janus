#!/usr/bin/env python3
# measure.py <run dir>: turns one pair's run (bench/run.sh) into structured measurements, written next to calls.vcf:
#   steps.tsv  per step (run.json also has each step's own report lines, e.g. reads and bases from pack): start/end (unix s), wall s, peak RSS KiB, user s, sys s, CPU s (user+sys), drive read/written GB,
#              mean drive utilization %, max drive composite / controller temperature and CPU Tctl (C)
#   files.tsv  every file the run produced, with its size in bytes and whether it is kept (calls.vcf, logs) or moved to trash
#   run.json   the pair, how it was run (note.txt, if any), machine, code commit, totals (wall, peak RSS, user, sys, CPU s, drive I/O, temperatures), the steps
#              and the file sizes, all in one file
# Sources: work/steps.txt (pipeline.sh's /usr/bin/time per step), work/done.<step> modification times (step ends),
# start.txt / total.txt (the runner's /usr/bin/time of the whole pipeline), monitor/{iostat,vmstat,temp}.txt.
import sys, os, re, json, glob, datetime, subprocess, platform

D = os.path.abspath(sys.argv[1]); W = os.path.join(D, "work"); M = os.path.join(D, "monitor")
KEEP = {"calls.vcf"}

def iostat():
    """[(unix s, read MB/s, write MB/s, util %)] of nvme0n1, one per 5 s report (the first, since-boot report dropped)"""
    out, t, first = [], None, True
    for l in open(os.path.join(M, "iostat.txt")):
        f = l.split()
        if len(f) == 3 and f[2] in ("AM", "PM") and "/" in f[0]:
            t = datetime.datetime.strptime(" ".join(f), "%m/%d/%Y %I:%M:%S %p").timestamp()
        elif f and f[0] == "nvme0n1" and t is not None:
            if first: first = False; continue
            out.append((t, float(f[2]), float(f[8]), float(f[-1])))
    return out

def temps():
    """[(unix s, drive composite C, drive controller C, CPU Tctl C or None)]"""
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
    """snape's own lines in the step's log ("[pack] 11337735 reads, 195114129352 bases", "[call] ... records, ... PASS")"""
    for f in (os.path.join(W, "log." + n), os.path.join(D, "logs", "log." + n)):
        if os.path.exists(f):
            return [l.strip() for l in open(f, errors="replace") if l.startswith("[") and not l.startswith("[M::")]
    return []

steps = []
for l in open(os.path.join(W, "steps.txt")):
    f = l.split()                                   # <step> <wall> s <KB> KB <user> user <sys> sys
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

if os.path.exists(os.path.join(D, "total.txt")):
    start = float(open(os.path.join(D, "start.txt")).read())
    tot = open(os.path.join(D, "total.txt")).read().split()  # total <wall> s wall <KB> KB peak <user> s user <sys> s sys ...
    total = dict(wall_s=float(tot[1]), peak_rss_kib=int(tot[4]), user_s=float(tot[7]), sys_s=float(tot[10]))
else:                                                # a run made step by step (paused between steps): totals of the steps
    start = steps[0]["start"]
    total = dict(wall_s=round(sum(s["wall_s"] for s in steps), 2), peak_rss_kib=max(s["peak_rss_kib"] for s in steps),
                 user_s=round(sum(s["user_s"] for s in steps), 2), sys_s=round(sum(s["sys_s"] for s in steps), 2))
total["cpu_s"] = round(total["user_s"] + total["sys_s"], 2)
total["step_wall_s"] = round(sum(s["wall_s"] for s in steps), 2)
total["peak_step"] = max(steps, key=lambda s: s["peak_rss_kib"])["step"]
total["read_gb"], total["written_gb"], total["util_pct"], total["drive_c"], total["controller_c"], total["cpu_c"] = window(start, prev)

files = sorted(((os.path.getsize(p), os.path.relpath(p, W)) for p in glob.glob(os.path.join(W, "*")) if os.path.isfile(p)), reverse=True)
files = [dict(file=f, bytes=b, kept=f in KEEP or f.startswith("log.")) for b, f in files if not f.startswith("done.")]
total["intermediate_gb"] = round(sum(f["bytes"] for f in files) / 1e9, 2)

cols = ["step", "start", "end", "wall_s", "peak_rss_kib", "user_s", "sys_s", "cpu_s", "read_gb", "written_gb", "util_pct",
        "drive_c", "controller_c", "cpu_c"]
with open(os.path.join(D, "steps.tsv"), "w") as o:
    o.write("\t".join(cols) + "\n")
    for s in steps + [dict(total, step="total", start=round(start, 2), end=round(prev, 2))]:
        o.write("\t".join("" if s.get(c) is None else str(s[c]) for c in cols) + "\n")
with open(os.path.join(D, "files.tsv"), "w") as o:
    o.write("file\tbytes\tkept\n")
    for f in files: o.write(f"{f['file']}\t{f['bytes']}\t{'yes' if f['kept'] else 'trash'}\n")

S = "/home/yuio/snape"; pair = os.path.basename(D)
reads = {s: dict(file=os.path.realpath(f"{S}/hifi/{pair}/{s}.fastq.gz"), bytes=os.path.getsize(f"{S}/hifi/{pair}/{s}.fastq.gz")) for s in "nt"}
for st in steps:                                     # read and base counts from pack's report
    if st["step"] in ("pack_n", "pack_t"):
        for l in st["report"]:
            m = re.match(r"\[pack\] (\d+) reads, (\d+) bases", l)
            if m: reads[st["step"][-1]].update(reads=int(m.group(1)), bases=int(m.group(2)))
cpu = next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")), "")
mem = int(next(l.split()[1] for l in open("/proc/meminfo") if l.startswith("MemTotal")))
commit = subprocess.run(["git", "-C", S, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
note = open(os.path.join(D, "note.txt")).read().strip() if os.path.exists(os.path.join(D, "note.txt")) else ""
json.dump(dict(pair=pair, commit=commit, note=note, started=datetime.datetime.fromtimestamp(start).isoformat(timespec="seconds"),
               threads=16, cores="taskset 0-15", machine=dict(cpu=cpu, mem_kib=mem, kernel=platform.release()),
               reads=reads, total=total, steps=steps, files=files),
          open(os.path.join(D, "run.json"), "w"), indent=1)
print(f"{pair}: {total['wall_s']} s wall, {total['cpu_s']} CPU s, peak {total['peak_rss_kib'] / 2**20:.1f} GiB ({total['peak_step']}), "
      f"intermediates {total['intermediate_gb']} GB")
