# the verdict on the coverage runs, trained at the same depth: for every depth N and labeled pair p, the five-seed trees
# trained on the N x records of the other labeled pairs (labels from their truth sets, as in verdict.py), the cut chosen
# on those training records for the depth's F-beta (params.json verdict_f_beta), applied to p's N x calls; PASS records
# scored against p's truth with score.py's rule. Writes verdict/coverage_matched/<pair>_<N>x.vcf and scores.tsv
import sys, os, json, numpy as np
sys.argv = sys.argv[:1]
exec(open("/home/yuio/snape/bench/verdict.py").read().split("t0 = time.time()")[0])
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
