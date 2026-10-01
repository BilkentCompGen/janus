# verdict.py [first seed]: the verdict, a learned component (boosted trees), always applied held out.
# Training pairs (CASTLE x4, HG008): for each one, trees (5 seeds averaged) are fitted on the records of the other training
# pairs (all records pooled), the cut on the probability is the one with the best F1 over those records, and the held-out
# pair is scored. Evaluation pairs (HCC1395, COLO829: truth sets that are incomplete, so they are scored but never trained
# on) and unlabeled pairs are scored by trees fitted on all five training pairs, with the cut chosen the same way.
# Labels: a record is true when it matches a truth event (score.py). bench/verdict/folds.json keeps, per scored pair, the
# training pairs, the cut, the training and inference times, and every record's probability and label.
# Duplicates within 500 bp are marked by the probability. Reads bench/runs/full/<pair>/calls.vcf, writes
# bench/verdict/<pair>.vcf (FILTER PASS / LOW / DUPL, INFO VP = the held-out probability) and bench/verdict/log.txt.
import sys, time, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
exec(open("/home/yuio/snape/bench/score.py").read().split("T, C = recs")[0])
H = "/home/yuio/snape/bench"; TR = "/home/yuio/snape/hifi/truth"; S0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
TRUTH = {p: f"{TR}/castle_without_snape/{p}/None_minda_ensemble.vcf" for p in ["hcc1954", "h2009", "hcc1937", "h1437"]}
TRUTH.update(hg008=f"{TR}/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz",
             hcc1395=f"{TR}/hcc1395/truth.vcf", colo829=f"{TR}/colo829/truth.vcf")
TRAIN = ["hcc1954", "h2009", "hcc1937", "h1437", "hg008"]
EVAL = ["hcc1395", "colo829"]
UNLABELED = ["hg009", "h209", "h2126"]
LABELED = TRAIN + EVAL
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
    """trees fitted on the records of the pairs in train (pooled) and the cut with the best F1 over those records"""
    X = np.vstack([D[q][1][D[q][4]] for q in train]); y = np.concatenate([D[q][3][D[q][4]] for q in train])
    t = time.time()
    ms = [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(S0, S0 + 5)]
    ts = time.time() - t
    pt = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)
    cut = max(CUTS, key=lambda c: f1(pt, y, c))
    return ms, cut, ts, len(y)

t0 = time.time()
D = {p: records(p) for p in LABELED + UNLABELED}
log = open(f"{H}/verdict/log.txt", "w"); folds = {}
for p in LABELED + UNLABELED:
    tr = [q for q in TRAIN if q != p]
    ms, cut, ts, nt = fit(tr)
    L, X, K = D[p][0], D[p][1], D[p][2]
    t = time.time(); pr = np.mean([m.predict_proba(X)[:, 1] for m in ms], 0); ti = time.time() - t; d = mark(K, pr)
    folds[p] = dict(role="held-out training pair" if p in TRAIN else "evaluation pair" if p in EVAL else "unlabeled", train=tr,
                    cut=round(float(cut), 2), train_s=round(ts, 3), train_records=nt, infer_s=round(ti, 4), records=len(L),
                    prob=[round(float(v), 4) for v in pr], label=None if D[p][3] is None else [int(v) for v in D[p][3]],
                    dupl=[int(v) for v in d], features=len(X[0]))
    with open(f"{H}/verdict/{p}.vcf", "w") as o:
        o.write("".join(l for l in open(f"{H}/runs/full/{p}/calls.vcf") if l[0] == "#"))
        for l, q, u in zip(L, pr, d):
            f = l.split("\t"); f[6] = "DUPL" if u else "PASS" if q >= cut else "LOW"; f[7] += ";VP=%.4f" % q; o.write("\t".join(f))
    s = f"{p} {'held-out' if p in TRAIN else 'evaluation' if p in EVAL else 'unlabeled'} PASS {int(((pr >= cut) & ~d).sum())} cut {cut:.2f} DUPL {int(d.sum())} records {len(L)}"
    print(s); log.write(s + "\n")
s = f"verdict ({len(TRAIN)} training pairs held out, {len(EVAL)} evaluation pairs, {len(UNLABELED)} unlabeled, 5 seeds) {time.time() - t0:.1f} s"
print(s); log.write(s + "\n")
import json; json.dump(folds, open(f"{H}/verdict/folds.json", "w"))
