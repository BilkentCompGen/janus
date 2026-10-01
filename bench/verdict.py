# verdict.py [first seed]: the verdict, a learned component (boosted trees), always applied held out.
# Labeled pairs (a truth set exists): for each one, trees (5 seeds averaged) are fitted on the records of the other labeled
# pairs, every pair weighted equally (a record weighs 1 / its pair's record count), the cut on the probability is the one
# with the best mean F1 over those pairs, and the held-out pair is scored. Unlabeled pairs are scored by trees fitted on
# all labeled pairs, with the cut chosen the same way. Labels: a record is true when it matches a truth event (score.py).
# Duplicates within 500 bp are marked by the probability. Reads bench/runs/full/<pair>/calls.vcf, writes
# bench/verdict/<pair>.vcf (FILTER PASS / LOW / DUPL, INFO VP = the held-out probability) and bench/verdict/log.txt.
import sys, time, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
exec(open("/home/yuio/snape/bench/score.py").read().split("T, C = recs")[0])
H = "/home/yuio/snape/bench"; TR = "/home/yuio/snape/hifi/truth"; S0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
TRUTH = {p: f"{TR}/castle_without_snape/{p}/None_minda_ensemble.vcf" for p in ["hcc1954", "h2009", "hcc1937", "h1437"]}
TRUTH.update(hg008=f"{TR}/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz",
             hcc1395=f"{TR}/hcc1395/truth.vcf", colo829=f"{TR}/colo829/truth.vcf")
LABELED = list(TRUTH)
UNLABELED = ["hg009", "h209", "h2126"]
HGB = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=40, l2_regularization=1.0)
CUTS = np.arange(0.2, 0.81, 0.05)

def records(p):
    L, X, K, R = [], [], [], []
    for l in open(f"{H}/runs/full/{p}/calls.vcf"):
        if l[0] == "#": continue
        f = l.split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
        K.append((f[0], m.group(1) if m else f[0], int(f[1]), int(m.group(2)) if m else int(f[1]) + int(info["SVLEN"]), info["SVTYPE"], info["O"]))
        R.append((info["SVTYPE"], f[0], int(f[1]), m.group(1), int(m.group(2))) if m else (info["SVTYPE"], f[0], int(f[1]), f[0], int(info["END"])))
        X.append([float(v) for v in info["F"].split(",")]); L.append(l)
    y = None
    if p in TRUTH:
        by = {}
        for t in recs(TRUTH[p]): by.setdefault(t[1], []).append(t); by.setdefault(t[3], []).append(t)
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

def f1(p, y, c):
    return 2 * ((p >= c) & (y > 0)).sum() / max((p >= c).sum() + (y > 0).sum(), 1)

def fit(train):
    """trees fitted on the labeled pairs in train (equal pair weights) and the cut with the best mean F1 over them"""
    X = np.vstack([D[q][1][D[q][4]] for q in train]); y = np.concatenate([D[q][3][D[q][4]] for q in train])
    w = np.concatenate([np.full(D[q][4].sum(), 1 / D[q][4].sum()) for q in train])
    ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y, sample_weight=w) for s in range(S0, S0 + 5)]
    pq = {q: np.mean([m.predict_proba(D[q][1][D[q][4]])[:, 1] for m in ms], 0) for q in train}
    cut = max(CUTS, key=lambda c: np.mean([f1(pq[q], D[q][3][D[q][4]], c) for q in train]))
    return ms, cut

t0 = time.time()
D = {p: records(p) for p in LABELED + UNLABELED}
log = open(f"{H}/verdict/log.txt", "w")
for p in LABELED + UNLABELED:
    ms, cut = fit([q for q in LABELED if q != p])
    L, X, K = D[p][0], D[p][1], D[p][2]; pr = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0); d = mark(K, pr)
    with open(f"{H}/verdict/{p}.vcf", "w") as o:
        o.write("".join(l for l in open(f"{H}/runs/full/{p}/calls.vcf") if l[0] == "#"))
        for l, q, u in zip(L, pr, d):
            f = l.split("\t"); f[6] = "DUPL" if u else "PASS" if q >= cut else "LOW"; f[7] += ";VP=%.4f" % q; o.write("\t".join(f))
    s = f"{p} {'held-out' if p in TRUTH else 'unlabeled'} PASS {int(((pr >= cut) & ~d).sum())} cut {cut:.2f} DUPL {int(d.sum())} records {len(L)}"
    print(s); log.write(s + "\n")
s = f"verdict ({len(LABELED)} labeled pairs held out + {len(UNLABELED)} unlabeled, 5 seeds) {time.time() - t0:.1f} s"
print(s); log.write(s + "\n")
