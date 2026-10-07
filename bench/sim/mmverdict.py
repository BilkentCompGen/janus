import os, sys
exec(open("/home/yuio/snape/bench/cov_verdict_matched.py").read().split("out = open(")[0])
p = os.environ["PAIR"]; n = 20
D = {q: crecords(f"{C}/{q}_{n}x/calls.vcf", q) for q in TRAIN if q != p}
X = np.vstack([D[q][1][D[q][4]] for q in D]); y = np.concatenate([D[q][3][D[q][4]] for q in D])
ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(5)]
pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0); cut = max(CUTS, key=lambda c: fb(pt, y, c, 1.0))
for name, f in [("current run", f"{C}/{p}_{n}x/calls.vcf")] + [(a.split("@")[0], a.split("@")[1]) for a in os.environ["VCFS"].split(";")]:
    L, Xc, K, _, _ = crecords(f, p)
    pr = np.mean([m.predict_proba(Xc)[:, 1] for m in ms], 0); du = mark(K, pr); v = f"/tmp/mmv_{p}.vcf"
    with open(v, "w") as o:
        o.write("".join(l for l in open(f) if l[0] == "#"))
        for l, r, u in zip(L, pr, du):
            g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if r >= cut else "LOW"; o.write("\t".join(g))
    k, tpc, tpt, P, R, F = score(v, p)
    print(f"{p} {name:32s} cut {cut:.2f} PASS {k:5d} P {P:.3f} R {R:.3f} F1 {F:.3f}", flush=True)
