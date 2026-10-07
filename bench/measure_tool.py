#!/usr/bin/env python3
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
