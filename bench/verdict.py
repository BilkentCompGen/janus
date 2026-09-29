# verdict.py <run dir> [first seed]: held-out verdict of record. For each pair, boosted trees (5 seeds averaged) are fitted on the
# other four pairs' calls.vcf records (labels: match to the no-snape consensus, NIST for hg008), the cut is chosen on those pairs,
# every record of the held-out pair is scored, and duplicates within 500 bp are marked by the probability. Writes <pair>/verdict.vcf.
import sys, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
exec(open("/home/yuio/snape/bench/score.py").read().split("T, C = recs")[0])
W, S0 = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 0
PAIRS = ["hcc1954", "hcc1937", "hg008", "h1437", "h2009"]
NIST = "/home/yuio/snape/hifi/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz"

def records(p):
    L, X, K, R = [], [], [], []
    for l in open(f"{W}/{p}/calls.vcf"):
        if l[0] == "#": continue
        f = l.split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
        K.append((f[0], m.group(1) if m else f[0], int(f[1]), int(m.group(2)) if m else int(f[1]) + int(info["SVLEN"]), info["SVTYPE"], info["O"]))
        R.append((info["SVTYPE"], f[0], int(f[1]), m.group(1), int(m.group(2))) if m else (info["SVTYPE"], f[0], int(f[1]), f[0], int(info["END"])))
        X.append([float(v) for v in info["F"].split(",")]); L.append(l)
    T = recs(NIST if p == "hg008" else f"/home/yuio/snape/hifi/truth/castle_without_snape/{p}/None_minda_ensemble.vcf")
    by = {}
    for t in T: by.setdefault(t[1], []).append(t); by.setdefault(t[3], []).append(t)
    y = np.array([any(hit(r, t) for t in by.get(r[1], [])) for r in R], float)
    keep = np.array(["\tDUPL\t" not in l for l in L])
    return L, np.array(X), K, y, keep

def mark(K, pr):
    ix = sorted(range(len(K)), key=lambda i: (K[i][:4], i)); d = np.zeros(len(K), bool)
    for a in range(len(ix)):
        x = K[ix[a]]
        for b in ix[a + 1:]:
            y = K[b]
            if y[:2] != x[:2] or y[2] > x[2] + 500: break
            if abs(y[3] - x[3]) <= 500 and (x[4] != y[4] or x[5] == y[5]): d[ix[a] if pr[ix[a]] < pr[b] else b] = True
    return d

D = {p: records(p) for p in PAIRS}
for p in D:
    Xt = np.vstack([D[q][1][D[q][4]] for q in D if q != p]); yt = np.concatenate([D[q][3][D[q][4]] for q in D if q != p])
    ms = [HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=40, l2_regularization=1.0, random_state=s).fit(Xt, yt) for s in range(S0, S0 + 5)]
    pt = np.mean([m.predict_proba(Xt)[:, 1] for m in ms], 0); cut = max(np.arange(0.2, 0.81, 0.05), key=lambda c: 2 * ((pt >= c) & (yt > 0)).sum() / max((pt >= c).sum() + (yt > 0).sum(), 1))
    L, X, K = D[p][0], D[p][1], D[p][2]; pr = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0); d = mark(K, pr)
    with open(f"{W}/{p}/verdict.vcf", "w") as o:
        o.write("".join(l for l in open(f"{W}/{p}/calls.vcf") if l[0] == "#"))
        for l, q, u in zip(L, pr, d): f = l.split("\t"); f[6] = "DUPL" if u else "PASS" if q >= cut else "LOW"; o.write("\t".join(f))
    print(p, "held-out PASS", int(((pr >= cut) & ~d).sum()), "cut %.2f" % cut, "DUPL", int(d.sum()))
