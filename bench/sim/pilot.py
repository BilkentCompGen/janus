import glob, os, sys
exec(open("/home/yuio/snape/bench/cov_verdict_matched.py").read().split("out = open(")[0])
M = f"{H}/runs/sim"
def sim(n):
    X, y = [], []
    for d in sorted(glob.glob(f"{M}/*/x{n}/calls.vcf")):
        T = recs(os.path.dirname(os.path.dirname(d)) + "/sim.truth.vcf"); by = {}
        for t in T: by.setdefault(t[1], []).append(t); by.setdefault(t[3], []).append(t)
        BY["_sim"] = by
        L, Xs, K, ys, keep = crecords(d, "_sim")
        X.append(Xs[keep]); y.append(ys[keep])
        print(f"sim {d.split('/')[-3]} {n}x: {len(L)} records, {int(ys.sum())} true", flush=True)
    return np.vstack(X), np.concatenate(y)
real = {}
for l in list(open(f"{H}/verdict/coverage_matched/scores.tsv"))[1:]:
    f = l.split("\t"); real[(f[0], int(f[1]))] = (f[9], f[10], f[11].strip())
print("pair\tdepth\tsim_cut\tsim_pass\tsim_P\tsim_R\tsim_F1\t|\treal_P\treal_R\treal_F1")
for n in [int(a) for a in os.environ["SIMDEPTHS"].split()]:
    X, y = sim(n)
    sh = float(os.environ.get("SIMPOS", "0")); w = np.ones(len(y))
    if sh > 0: w[y > 0] = sh * (y == 0).sum() / ((1 - sh) * (y > 0).sum())
    ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y, sample_weight=w) for s in range(5)]
    pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)
    def wf(c): tp = (w * ((pt >= c) & (y > 0))).sum(); return 2 * tp / max((w * (pt >= c)).sum() + (w * (y > 0)).sum(), 1e-9)
    cut = max(CUTS, key=wf)
    for p in TRAIN:
        f = f"{C}/{p}_{n}x/calls.vcf"
        if not os.path.exists(f): continue
        L, Xc, K, _, _ = crecords(f, p)
        pr = np.mean([m.predict_proba(Xc)[:, 1] for m in ms], 0); du = mark(K, pr); v = f"/tmp/simv_{p}_{n}.vcf"
        with open(v, "w") as o:
            o.write("".join(l for l in open(f) if l[0] == "#"))
            for l, r, u in zip(L, pr, du):
                g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if r >= cut else "LOW"; o.write("\t".join(g))
        k, tpc, tpt, P, R, F = score(v, p); a = real.get((p, n), ("-", "-", "-"))
        print(f"{p}\t{n}\t{cut:.2f}\t{k}\t{P:.3f}\t{R:.3f}\t{F:.3f}\t|\t{a[0]}\t{a[1]}\t{a[2]}", flush=True)
