#!/usr/bin/env python3
# compare_groups.py: today's snape (bench/runs/full) against the per-distance-group select (bench/runs/groups) on the 4
# CASTLE pairs + HG008. Both get the same verdict (bench/verdict.py's code: boosted trees, 5 seeds, each pair held out once,
# trained on the other four, records pooled, cut by best F1 on the training records) and the same scoring (score.py's
# matching at each pair's own held-out cut, against the CASTLE consensus without snape / NIST for HG008). Writes
# runs/groups/compare/<version>/<pair>.vcf (the verdict VCFs, for minda) and compare.json, and prints the table.
import json, os, sys, time, numpy as np
H = "/home/yuio/snape/bench"
src = open(f"{H}/verdict.py").read().split("t0 = time.time()")[0]
PAIRS = ["hcc1954", "h2009", "hcc1937", "h1437", "hg008"]
OUT = f"{H}/runs/groups/compare"
res = {}
for ver in ("full", "groups"):
    g = {}
    exec(src.replace('f"{H}/runs/full/{p}/calls.vcf"', 'f"{H}/runs/' + ver + '/{p}/calls.vcf"'), g)
    D = {p: g["records"](p) for p in PAIRS}
    g["D"] = D
    os.makedirs(f"{OUT}/{ver}", exist_ok=True)
    res[ver] = {}
    for p in PAIRS:
        ms, cut, ts, nt = g["fit"]([q for q in PAIRS if q != p])
        L, X, K = D[p][0], D[p][1], D[p][2]
        pr = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0); d = g["mark"](K, pr)
        with open(f"{OUT}/{ver}/{p}.vcf", "w") as o:
            o.write("".join(l for l in open(f"{H}/runs/{ver}/{p}/calls.vcf") if l[0] == "#"))
            for l, q, u in zip(L, pr, d):
                f = l.split("\t"); f[6] = "DUPL" if u else "PASS" if q >= cut else "LOW"; f[7] += ";VP=%.4f" % q; o.write("\t".join(f))
        # event-level scoring: PASS records against the truth set (score.py's recs / hit)
        T = g["recs"](g["TRUTH"][p]); C = g["recs"](f"{OUT}/{ver}/{p}.vcf")
        bt, bc = {}, {}
        for x in T: bt.setdefault(x[1], []).append(x); bt.setdefault(x[3], []).append(x)
        for x in C: bc.setdefault(x[1], []).append(x); bc.setdefault(x[3], []).append(x)
        P = sum(any(g["hit"](x, y) for y in bt.get(x[1], [])) for x in C) / max(len(C), 1)
        R = sum(any(g["hit"](x, y) for y in bc.get(x[1], [])) for x in T) / len(T)
        rj = json.load(open(f"{H}/runs/{ver}/{p}/run.json"))
        st = {s["step"]: s for s in rj["steps"]}
        res[ver][p] = dict(records=len(L), pass_=len(C), cut=round(float(cut), 2), precision=round(P, 4), recall=round(R, 4),
                           f1=round(2 * P * R / (P + R), 4) if P + R else 0, wall_s=rj["total"]["wall_s"], cpu_s=rj["total"]["cpu_s"],
                           peak_gb=round(rj["total"]["peak_rss_kib"] * 1024 / 1e9, 2), peak_step=rj["total"]["peak_step"],
                           select_s=st["select"]["wall_s"], local_s=st["local"]["wall_s"], call_s=st["call"]["wall_s"],
                           select_report=st["select"]["report"], call_report=st["call"]["report"],
                           intermediate_gb=rj["total"]["intermediate_gb"])
json.dump(res, open(f"{OUT}/compare.json", "w"), indent=1)
print("%-8s %-7s %7s %6s %5s %6s %6s %6s %8s %8s %7s %7s %7s" % ("pair", "version", "records", "PASS", "cut", "P", "R", "F1", "wall min", "CPU h", "peak GB", "select s", "local s"))
for p in PAIRS:
    for ver in ("full", "groups"):
        r = res[ver][p]
        print("%-8s %-7s %7d %6d %5.2f %6.3f %6.3f %6.3f %8.1f %8.2f %7.1f %7.1f %7.1f" % (p, ver, r["records"], r["pass_"], r["cut"], r["precision"], r["recall"], r["f1"], r["wall_s"] / 60, r["cpu_s"] / 3600, r["peak_gb"], r["select_s"], r["local_s"]))
for ver in ("full", "groups"):
    print(ver, "mean F1 %.4f" % np.mean([res[ver][p]["f1"] for p in PAIRS]))
