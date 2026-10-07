import json, os
exec(open("/home/yuio/snape/bench/cov_verdict_matched.py").read().split("out = open(")[0])
R = f"{H}/runs"; X = f"{R}/mmfull"; CAS = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]; EXTRA = ["hg008", "colo829"]
for p in EXTRA:
    TT[p] = recs(TRUTH[p]); by = {}
    for t in TT[p]: by.setdefault(t[1], []).append(t); by.setdefault(t[3], []).append(t)
    BY[p] = by

def lopo(name, S, beta, out):
    D = {p: crecords(f, p) for p, f in S.items() if os.path.exists(f)}
    os.makedirs(out, exist_ok=True); tot = []; models = {}
    for p in D:
        tr = [q for q in D if q != p]
        Xt = np.vstack([D[q][1][D[q][4]] for q in tr]); y = np.concatenate([D[q][3][D[q][4]] for q in tr])
        ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(Xt, y) for s in range(5)]
        pt = np.mean([m.predict_proba(Xt)[:, 1] for m in ms], 0); cut = max(CUTS, key=lambda c: fb(pt, y, c, beta))
        models[p] = (ms, cut)
        tot.append(apply(name, p, S[p], ms, cut, out, D[p]))
    print(f"{name}\tmean ({len(tot)} pairs)\t\t\t\t\t\t{np.mean(tot):.4f}", flush=True)
    return models

def apply(name, p, f, ms, cut, out, d=None):
    L, Xc, K = (d or crecords(f, p))[:3]
    os.makedirs(out, exist_ok=True)
    pr = np.mean([m.predict_proba(Xc)[:, 1] for m in ms], 0); du = mark(K, pr); v = f"{out}/{p}.vcf"
    with open(v, "w") as o:
        o.write("".join(l for l in open(f) if l[0] == "#"))
        for l, r, u in zip(L, pr, du):
            g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if r >= cut else "LOW"; g[7] += ";VP=%.4f" % r; o.write("\t".join(g))
    k, tpc, tpt, P, Rr, F = score(v, p)
    print(f"{name}\t{p}\t{len(L)}\t{cut:.2f}\t{k}\t{P:.4f}\t{Rr:.4f}\t{F:.4f}", flush=True)
    return F

print("setting\tpair\trecords\tcut\tpass\tprecision\trecall\tf1")
cur = {p: f"{X}/{p}/full_m3/calls.vcf" for p in CAS}; cur.update({p: f"{R}/full/{p}/calls.vcf" for p in EXTRA})
lopo("full depth, current (3/3)", cur, 1.0, f"{X}/verdict/full_m3")
lopo("full depth, tumor 2 / normal 3", {p: f"{X}/{p}/full_m2/calls.vcf" for p in CAS + EXTRA}, 1.0, f"{X}/verdict/full_m2")
m = lopo("10x, current (2/2)", {p: f"{X}/{p}/x10_n2t2/calls.vcf" for p in CAS}, 0.75, f"{X}/verdict/x10_n2t2")
lopo("10x, tumor 2 / normal 3", {p: f"{X}/{p}/x10_n3t2/calls.vcf" for p in CAS}, 0.75, f"{X}/verdict/x10_n3t2")
for p in ["hcc1954", "hcc1395"]:
    for tag, f in (("10x, tumor 1 (record >= 2), verdict of 2/2", f"{X}/{p}/x10_n2t1/calls.vcf"), ("10x, tumor 1 (record >= 1), verdict of 2/2", f"{X}/{p}/x10_n2t1/calls_m1.vcf")):
        if os.path.exists(f) and p in m: apply(tag, p, f, m[p][0], m[p][1], f"{X}/verdict/x10_n2t1_{f[-6:-4]}")

print("\ntime and memory from the tumor count step to the calls (16 threads, alone)")
print("pair\tsetting\tseconds\tpeak_gb\tsites\tpieces_bp")
for p in CAS + EXTRA:
    for v in ("full_m3", "full_m2"):
        d = f"{X}/{p}/{v}"
        if not os.path.exists(f"{d}/steps.txt"): continue
        st = [l.split() for l in open(f"{d}/steps.txt")]
        sel = [l for l in open(f"{d}/log.select") if l.startswith("[select]")][-1]
        print(f"{p}\t{v}\t{sum(float(x[1]) for x in st):.0f}\t{max(int(x[3]) for x in st) / 1048576:.1f}\t{sel.split()[1]}\t{sel.split()[-2]}")
    j = json.load(open(f"{R}/full/{p}/run.json"))["total"]
    print(f"{p}\tofficial full run (all steps)\t{j['wall_s']:.0f}\t{j['peak_rss_kib'] / 1048576:.1f}")
