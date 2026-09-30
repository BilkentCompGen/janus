# ml.py <tag> [pairs] (XGBoost + MLP for inspiration only, never the method): XGBoost and an MLP on the per-record measurements of calls.<tag>.vcf, leave-one-pair-out, labels as
# bench/verdict.py derives them (match within 500 bp to the no-snape CASTLE consensus, NIST for hg008).
# FOR INSPIRATION ONLY (the user allowed ML/DL to observe, not as the method): prints feature medians for true/false records,
# XGBoost gain importances and held-out record-level scores, and writes calls.<tag>.xgb.vcf / calls.<tag>.mlp.vcf whose PASS is
# the held-out model's decision (cut chosen on the training pairs) so the fair harness can score them, labelled ML/DL.
import sys, math, collections as C, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.argv = [sys.argv[0], sys.argv[1]] + sys.argv[2:]
sys.path.insert(0, "/home/yuio/snape/bench"); import rerule as RR
tag = sys.argv[1]
PAIRS = sys.argv[2].split(",") if len(sys.argv) > 2 else ["hcc1954", "hcc1937", "hg008", "h1437", "h2009"]
NAMES = ["INS", "DEL", "DUP", "INV", "BND", "log_len", "mq", "pieces", "tdp", "ndp", "gcount_t", "span_n", "span_t", "same_n", "same_t", "near_n", "near_t", "diff_n", "diff_t",
         "seq_t", "same_n_frac", "vaf_span", "vaf_depth", "nmin2_ge50", "nmin2_log", "nstable", "tstable", "nclu", "nsite", "top", "inside_site", "vntr", "site_len", "no_measure"]
def feats(r):
    return [r["ty"] == "INS", r["ty"] == "DEL", r["ty"] == "DUP", r["ty"] == "INV", r["ty"] == "BND", math.log2(r["L"] + 1), r["mq"], r["pc"], r["tdp"], r["ndp"], r["gc"],
            len(r["dn"]), len(r["dt"]), r["same_n"], r["same_t"], r["near_n"], sum(20 <= abs(x) < 50 for x in r["dt"]), sum(abs(x) >= 50 for x in r["dn"]), sum(abs(x) >= 50 for x in r["dt"]),
            r["sq"][1], r["same_n"] / max(len(r["dn"]), 1), r["vaf"], r["vafd"], r["nmin2"] >= 50, math.log2(min(r["nmin2"], 10**6) + 1), r["nstable"], r["tstable"], r["nclu"], r["nsite"], r["top"],
            r["a"] <= r["q"] <= r["b"], r["vn"], r["b"] - r["a"], not r["dn"] and not r["dt"]]
D = {}
for p in PAIRS:
    R = RR.parse(p)
    T = RR.recs(RR.NIST if p == "hg008" else f"/home/yuio/snape/hifi/truth/castle_without_snape/{p}/None_minda_ensemble.vcf")
    by = {}
    for t in T: by.setdefault(t[1], []).append(t); by.setdefault(t[3], []).append(t)
    y = np.array([any(RR.hit(r["key"], t) for t in by.get(r["key"][1], [])) for r in R], float)
    D[p] = (R, np.array([feats(r) for r in R], float), y)
print("records per pair (n, true):", {p: (len(D[p][2]), int(D[p][2].sum())) for p in D})
Xa = np.vstack([D[p][1] for p in D]); ya = np.concatenate([D[p][2] for p in D])
print("\nfeature medians (true | false) over all pairs' records:")
for i, n in enumerate(NAMES):
    t, f = Xa[ya > 0, i], Xa[ya == 0, i]
    print("  %-12s true med %9.3f q25 %9.3f q75 %9.3f | false med %9.3f q25 %9.3f q75 %9.3f" % (n, np.median(t), np.quantile(t, .25), np.quantile(t, .75), np.median(f), np.quantile(f, .25), np.quantile(f, .75)))
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import xgboost as xgb
def f1(y, pr, c): tp = ((pr >= c) & (y > 0)).sum(); return 2 * tp / max((pr >= c).sum() + (y > 0).sum(), 1)
imp = C.Counter()
for name in ["xgb", "mlp"]:
    print("\n==", name.upper(), "leave-one-pair-out (record level; LABELLED ML/DL, for inspiration, not a result)")
    for p in PAIRS:
        Xt = np.vstack([D[q][1] for q in D if q != p]); yt = np.concatenate([D[q][2] for q in D if q != p])
        if name == "xgb":
            m = xgb.XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, min_child_weight=5, reg_lambda=1.0, n_jobs=8, verbosity=0).fit(Xt, yt)
            for k, v in zip(NAMES, m.feature_importances_): imp[k] += v / len(PAIRS)
        else:
            m = make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(64, 32), alpha=1e-3, max_iter=400, random_state=0)).fit(Xt, yt)
        pt = m.predict_proba(Xt)[:, 1]; cut = max(np.arange(0.2, 0.81, 0.05), key=lambda c: f1(yt, pt, c))
        R, X, y = D[p]; pr = m.predict_proba(X)[:, 1]
        keep = [q >= cut for q in pr]; dup = RR.dedupe(R, keep)
        print("  %-8s held-out record F1 %.3f at cut %.2f ; PASS %d of %d records, true %d" % (p, f1(y, pr, cut), cut, int(sum(keep)), len(y), int(y.sum())))
        with open(f"{RR.CB}/{p}/calls.{tag}.{name}.vcf", "w") as o:
            o.write("".join(l for l in open(f"{RR.CB}/{p}/calls.{tag}.vcf") if l[0] == "#"))
            for r, k, d in zip(R, keep, dup):
                f = list(r["f"]); f[6] = "DUPL" if k and d else "PASS" if k else "LOW"; o.write("\t".join(f) + "\n")
    if name == "xgb":
        print("  XGBoost gain importance (mean over folds):")
        for k, v in imp.most_common(): print("    %-12s %.3f" % (k, v))
