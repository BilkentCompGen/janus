"""Classification of the records of a janus VCF with gradient-boosted trees."""
import argparse, pickle
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from . import features as F

HGB = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=40, l2_regularization=1.0)
CUTS = np.arange(0.2, 0.81, 0.05)
SEEDS = 5


def fit_models(X, y):
    return [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X, y) for s in range(SEEDS)]


def predict(ms, X):
    return np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)


def f1(p, y, c):
    return 2 * ((p >= c) & (y > 0)).sum() / max((p >= c).sum() + (y > 0).sum(), 1)


def train(sets, fitm=None):
    """Train on a list of (group, X, y) sets, in this order. The probability threshold of each minimum read support is chosen on
    held-out predictions: each group is predicted by models trained on the other groups."""
    import time
    fitm = fitm or fit_models
    X = np.vstack([x for _, x, _ in sets]); y = np.concatenate([v for _, _, v in sets]); g = np.round(X[:, 56])
    t = time.time(); ms = fitm(X, y); ts = time.time() - t
    oof = np.zeros(len(y)); o = 0; sl = []
    for _, x, _ in sets: sl.append(slice(o, o + len(x))); o += len(x)
    for q in dict.fromkeys(k for k, _, _ in sets):
        mi = fitm(np.vstack([x for k, x, _ in sets if k != q]), np.concatenate([v for k, _, v in sets if k != q]))
        for (k, x, _), s in zip(sets, sl):
            if k == q: oof[s] = predict(mi, x)
    cuts = {int(m): float(max(CUTS, key=lambda c: f1(oof[g == m], y[g == m], c))) for m in set(g)}
    return dict(models=ms, cuts=cuts, cut=float(max(CUTS, key=lambda c: f1(oof, y, c))), train=list(dict.fromkeys(k for k, _, _ in sets)), train_s=ts, n=len(y))


def threshold(model, mt):
    return model["cuts"].get(int(mt), model["cut"])


def classify(model, X, K, hist, mt=None):
    """Probabilities, duplicate marks, and the threshold for the records of one run."""
    mt = F.minimum(X) if mt is None else mt
    pr = predict(model["models"], F.features(X, hist, mt)); return pr, F.duplicates(K, pr), threshold(model, mt)


def write(out, H, L, pr, du, cut):
    with open(out, "w") as o:
        o.write("".join(H))
        for l, r, u in zip(L, pr, du):
            g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if r >= cut else "LOW"; g[7] += ";VP=%.4f" % r; o.write("\t".join(g))


def main():
    a = argparse.ArgumentParser(description="Classify the records of a janus VCF: PASS, LOW, or DUPL, with the probability VP.")
    a.add_argument("-i", required=True, help="calls.vcf of janus call"); a.add_argument("-o", required=True, help="output VCF")
    a.add_argument("-m", required=True, help="model (models/janus.pkl)"); a.add_argument("-H", required=True, help="t.cnt.hist of the tumor sample")
    a.add_argument("-s", type=int, help="minimum read support of the run (default: from the records)")
    a = a.parse_args()
    H, L, X, K, _, _, _ = F.records(a.i); model = pickle.load(open(a.m, "rb"))
    pr, du, cut = classify(model, X, K, a.H, a.s); write(a.o, H, L, pr, du, cut)
    print(f"{a.o}: {int(((pr >= cut) & ~du).sum())} PASS, {int(((pr < cut) & ~du).sum())} LOW, {int(du.sum())} DUPL of {len(L)} records, threshold {cut:.2f}")


if __name__ == "__main__":
    main()
