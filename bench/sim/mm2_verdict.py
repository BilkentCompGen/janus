import os
exec(open("/home/yuio/snape/bench/cov_verdict_matched.py").read().split("out = open(")[0])
R = f"{H}/runs"
SETS = {"current": {}, "MM=2": {}}
for p in ["hcc1954", "hcc1395"]:
    SETS["current"][p] = f"{R}/coverage/{p}_20x/calls.vcf"; SETS["MM=2"][p] = f"{R}/trace/{p}_20x_mm2/calls_mv.vcf"
for p in ["h2009", "hcc1937", "h1437"]:
    SETS["current"][p] = f"{R}/mm2/{p}_20x/m3/calls.vcf"; SETS["MM=2"][p] = f"{R}/mm2/{p}_20x/m2/calls.vcf"
print("setting\tpair\tcut\tpass\tprecision\trecall\tf1")
for name, S in SETS.items():
    D = {p: crecords(f, p) for p, f in S.items() if os.path.exists(f)}
    tot = []
    for p in D:
        tr = [q for q in D if q != p]
        X = np.vstack([D[q][1][D[q][4]] for q in tr]); y = np.concatenate([D[q][3][D[q][4]] for q in tr])
        ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(5)]
        pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0); cut = max(CUTS, key=lambda c: fb(pt, y, c, 1.0))
        L, Xc, K = D[p][0], D[p][1], D[p][2]
        pr = np.mean([m.predict_proba(Xc)[:, 1] for m in ms], 0); du = mark(K, pr); v = f"/tmp/mm2v_{p}.vcf"
        with open(v, "w") as o:
            o.write("".join(l for l in open(S[p]) if l[0] == "#"))
            for l, r, u in zip(L, pr, du):
                g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if r >= cut else "LOW"; o.write("\t".join(g))
        k, tpc, tpt, P, Rr, F = score(v, p); tot.append(F)
        print(f"{name}\t{p}\t{cut:.2f}\t{k}\t{P:.3f}\t{Rr:.3f}\t{F:.3f}", flush=True)
    print(f"{name}\tmean\t\t\t\t\t{np.mean(tot):.3f}", flush=True)
