import os, sys
exec(open("/home/yuio/snape/bench/cov_verdict_matched.py").read().split("out = open(")[0])
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
