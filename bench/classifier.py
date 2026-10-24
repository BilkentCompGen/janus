#!/usr/bin/env python3
"""The classifier of janus: training on the records of the CASTLE pairs, the leave-one-pair-out classification of every run
(verdict), the exported models, and the supplementary experiments on the classifier.

    classifier.py verdict [--runs DIR --coverage DIR --purity DIR --out DIR --labels DIR --models DIR]
    classifier.py supp <auc|splits|ablation|features|grid|firstpass> [--out DIR] [--shard i/n]

A pair's records are labelled by the consensus of the other call sets (eval_sets/<pair>/all_without, or the published
11-call-set consensus with --labels); the model of a pair is trained on the other four CASTLE pairs (full depth, every
depth and every purity below 100%); HG008 and COLO829 are classified by the model trained on all five, which is also the
model shipped in models/."""
import glob, itertools, json, os, subprocess, sys, time
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tools import ALT, BENCH, ROOT, hit, recs

HGB = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=40, l2_regularization=1.0)
CUTS = np.arange(0.2, 0.81, 0.05)
SEEDS = 5
DUPT = 500
CASTLE = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]
EXTERNAL = {"hg008": f"{ROOT}/datasets/hifi/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz", "colo829": f"{ROOT}/datasets/hifi/truth/colo829/truth.vcf"}
DEPTHS = [5, 10, 15, 20, 30]
PUR = [20, 40, 60, 80, 100]
JANUS = f"{ROOT}/build/janus"


# ---------------------------------------------------------------- features (the python counterpart of src/model.cpp, used for training)
def hist(path):
    return {l.split()[0]: np.array(l.split()[1:], float) for l in open(path)}


def pct(h, f, lo):
    c = np.cumsum(h[lo:]) / max(h[lo:].sum(), 1)
    return float(lo + np.searchsorted(c, f))


def context(path, mt):
    h = hist(path); mn = int(np.nonzero(h["shared_normal"])[0][0])
    tm = pct(h["shared_tumor"], 0.5, mt); nm = pct(h["shared_normal"], 0.5, mn)
    return tm, nm, [np.log2(tm), np.log2(nm), pct(h["tumor_only"], 0.05, mt) / tm, pct(h["tumor_only"], 0.95, mt) / tm]


def minimum(X):
    return int(round(2 ** X[:, 1].min() - 1))


def features(X, path, mt):
    """The record features with the log2 counts normalised by the sample medians, the minimum read support and the four context values."""
    tm, nm, c4 = context(path, mt); Z = X.copy()
    for j in (1, 3, 39): Z[:, j] -= np.log2(tm + 1)
    for j in (26, 27, 33, 34, 38, 41, 47, 51): Z[:, j] -= np.log2(nm + 1)
    return np.hstack([Z, np.full((len(X), 1), mt, float), np.tile(c4, (len(X), 1))])


def records(path):
    """Lines, raw features, breakpoint keys and comparison tuples of a janus VCF."""
    L, X, K, R = [], [], [], []
    for l in open(path):
        if l[0] == "#": continue
        f = l.split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
        K.append((f[0], m.group(1) if m else f[0], int(f[1]), int(m.group(2)) if m else int(f[1]) + int(info["SVLEN"]), info["SVTYPE"], info["O"]))
        R.append((info["SVTYPE"], f[0], int(f[1]), m.group(1), int(m.group(2))) if m else (info["SVTYPE"], f[0], int(f[1]), f[0], int(info["END"])))
        X.append([float(v) for v in info["F"].split(",")]); L.append(l)
    return L, np.array(X), K, R


def mark(K, pr):
    """The record with the lower probability in each pair of duplicate records (as `janus classify` does)."""
    ix = sorted(range(len(K)), key=lambda i: (K[i][:4], i)); d = np.zeros(len(K), bool)
    for a in range(len(ix)):
        x = K[ix[a]]
        for b in ix[a + 1:]:
            y = K[b]
            if y[:2] != x[:2] or y[2] > x[2] + DUPT: break
            if abs(y[3] - x[3]) <= DUPT and (x[4] != y[4] or x[5] == y[5]): d[ix[a] if pr[ix[a]] < pr[b] else b] = True
    return d


# ---------------------------------------------------------------- training
def fit_models(X, y, hgb=None):
    return [HistGradientBoostingClassifier(**(hgb or HGB), random_state=s).fit(X, y) for s in range(SEEDS)]


def predict(ms, X):
    return np.mean([m.predict_proba(X)[:, 1] for m in ms], 0)


def f1(p, y, c):
    return 2 * ((p >= c) & (y > 0)).sum() / max((p >= c).sum() + (y > 0).sum(), 1)


def train(sets, fitm=fit_models):
    """Train on a list of (pair, X, y) sets. The threshold of each minimum read support is chosen on held-out predictions:
    each pair is predicted by models trained on the other pairs."""
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


def export(model, path):
    """The text model read by `janus classify`: a header, then every tree of every seed, one node per line
    (value feature threshold missing_left left right; leaves have feature -1). Doubles are written exactly (repr)."""
    with open(path, "w") as o:
        ms = model["models"]
        o.write("janus-model 1\n")
        o.write(f"features {ms[0].n_features_in_} seeds {len(ms)} cut {model['cut']!r} cuts " + " ".join(f"{m}:{c!r}" for m, c in sorted(model["cuts"].items())) + " train " + ",".join(model["train"]) + "\n")
        for i, m in enumerate(ms):
            o.write(f"seed {i} baseline {float(m._baseline_prediction[0, 0])!r} trees {len(m._predictors)}\n")
            for p in m._predictors:
                nodes = p[0].nodes
                o.write(f"tree {len(nodes)}\n")
                for n in nodes:
                    if n["is_leaf"]: o.write(f"{float(n['value'])!r} -1 0 0 0 0\n")
                    else: o.write(f"{float(n['value'])!r} {int(n['feature_idx'])} {float(n['num_threshold'])!r} {int(n['missing_go_to_left'])} {int(n['left'])} {int(n['right'])}\n")


# ---------------------------------------------------------------- the runs and their labels
class Data:
    def __init__(self, runs, cov, pur, labels=None):
        self.runs, self.cov, self.pur = runs, cov, pur
        self.truth = {p: f"{labels}/{p}/None_minda_ensemble.vcf" if labels else f"{BENCH}/eval_sets/{p}/all_without/None_minda_ensemble.vcf" for p in CASTLE}
        self.truth.update(EXTERNAL)
        self.pairs = [p for p in self.truth if os.path.exists(f"{runs}/{p}/calls.vcf")]
        self.TT, self.BY = {}, {}
        for p in self.pairs:
            self.TT[p] = recs(self.truth[p]); self.BY[p] = {}
            for t in self.TT[p]: self.BY[p].setdefault(t[1], []).append(t); self.BY[p].setdefault(t[3], []).append(t)
        self.D = {(p, k): self.load(f, p) for p in self.pairs for k, f in self.sets(p).items()}

    def sets(self, p):
        s = {"full": f"{self.runs}/{p}/calls.vcf"}
        s.update({f"{n}x": f"{self.cov}/{p}_{n}x/calls.vcf" for n in DEPTHS}); s.update({f"{q}p": f"{self.pur}/{p}_{q}p/calls.vcf" for q in PUR})
        return {k: f for k, f in s.items() if os.path.exists(f)}

    def load(self, path, p):
        """Lines, model inputs, keys, labels (a record matches a consensus SV) and the records that are not DUPL in calls.vcf."""
        L, X, K, R = records(path)
        y = np.array([any(hit(r, t) for t in self.BY[p].get(r[1], [])) for r in R], float)
        keep = np.array(["\tDUPL\t" not in l for l in L])
        return L, features(X, f"{os.path.dirname(path)}/t.cnt.hist", minimum(X)), K, y, keep

    def score(self, v, p):
        """Precision, recall and F1 of the PASS records of a classified VCF against the pair's consensus."""
        C = recs(v); bc = {}
        for c in C: bc.setdefault(c[1], []).append(c); bc.setdefault(c[3], []).append(c)
        tpc = sum(any(hit(c, t) for t in self.BY[p].get(c[1], [])) for c in C)
        tpt = sum(any(hit(t, c) for c in bc.get(t[1], []) + bc.get(t[3], [])) for t in self.TT[p])
        P = tpc / max(len(C), 1); R = tpt / max(len(self.TT[p]), 1)
        return len(C), tpc, tpt, P, R, 2 * P * R / max(P + R, 1e-9)

    def training(self, pairs):
        """The training sets of a model: every run of the pairs except the 100% purity runs (they repeat the 30x runs)."""
        return [(k[0], self.D[k][1][self.D[k][4]], self.D[k][3][self.D[k][4]]) for k in self.D if k[0] in pairs and k[1] != "100p"]


def classified(model_path, calls, hist, out):
    """`janus classify` on one run; returns the probabilities and duplicate marks it wrote, and the time it took."""
    t = time.time()
    subprocess.run([JANUS, "classify", "--calls", calls, "--hist", hist, "--model", model_path, "--out", out], check=True, capture_output=True)
    t = time.time() - t
    pr, du = [], []
    for l in open(out):
        if l[0] == "#": continue
        f = l.split("\t"); pr.append(float(f[7].split("VP=")[1].split(";")[0])); du.append(f[6] == "DUPL")
    return np.array(pr), np.array(du, bool), t


def verdict(runs=f"{BENCH}/runs/full", coverage=f"{BENCH}/runs/coverage", purity=f"{BENCH}/runs/purity", out=f"{BENCH}/verdict", labels=None, models=None):
    """One model per held-out CASTLE pair, applied by `janus classify` to the pair's full, coverage and purity runs; HG008 and
    COLO829 by the model of all five pairs. Writes out/<pair>.vcf, out/coverage_matched, out/purity_matched, folds.json and
    log.txt, and the models (text) into `models` (default: out)."""
    t0 = time.time(); data = Data(runs, coverage, purity, labels)
    O = f"{out}/coverage_matched"; os.makedirs(O, exist_ok=True); os.makedirs(f"{out}/purity_matched", exist_ok=True); os.makedirs(models or out, exist_ok=True)
    log = open(f"{out}/log.txt", "w"); folds = {}; fits = {}
    covs = open(f"{O}/scores.tsv", "w"); covs.write("pair\tdepth\tcut\ttrain_pairs\tpass\ttruth\ttp_calls\ttp_truth\tprecision\trecall\tf1\n")
    purs = open(f"{out}/purity_matched/scores.tsv", "w"); purs.write("pair\tpurity\tcut\ttrain_pairs\tpass\ttruth\ttp_calls\ttp_truth\tprecision\trecall\tf1\n")
    for p in data.pairs:
        tr = tuple(q for q in CASTLE if q != p)
        if tr not in fits:
            m = train(data.training(tr)); path = f"{models or out}/janus{'' if len(tr) == len(CASTLE) else '_without_' + p}.model"
            export(m, path); fits[tr] = (m, path)
        m, path = fits[tr]
        for q, k in [k for k in data.D if k[0] == p]:
            L, X, K, y, keep = data.D[(q, k)]; mt = int(round(X[0, 56])); cut = m["cuts"].get(mt, m["cut"])
            v = f"{out}/{p}.vcf" if k == "full" else f"{O}/{p}_{k}.vcf" if k[-1] == "x" else f"{out}/purity_matched/{p}_{k}.vcf"
            calls = data.sets(p)[k]; pr, d, ti = classified(path, calls, f"{os.path.dirname(calls)}/t.cnt.hist", v)
            if k == "full":
                folds[p] = dict(train=list(tr), cut=round(float(cut), 2), train_s=round(m["train_s"], 3), train_records=m["n"], infer_s=round(ti, 4), records=len(L),
                                prob=[float(x) for x in pr], label=[int(x) for x in y], dupl=[int(x) for x in d], features=len(X[0]))
                s = f"{p} {'held-out' if p in CASTLE else 'not in training'} PASS {int(((pr >= cut) & ~d).sum())} cut {cut:.2f} minimum {mt} DUPL {int(d.sum())} records {len(L)}"
                print(s); log.write(s + "\n")
            else:
                n, tpc, tpt, pp, rr, ff = data.score(v, p)
                row = f"{cut:.2f}\t{','.join(tr)}\t{n}\t{len(data.TT[p])}\t{tpc}\t{tpt}\t{pp:.4f}\t{rr:.4f}\t{ff:.4f}\n"
                (covs if k[-1] == "x" else purs).write(f"{p}\t{k[:-1]}\t{row}")
    s = f"verdict ({len(data.pairs)} pairs, {SEEDS} seeds, each CASTLE pair held out) {time.time() - t0:.1f} s"
    print(s); log.write(s + "\n")
    json.dump(folds, open(f"{out}/folds.json", "w"))


# ---------------------------------------------------------------- the supplementary experiments on the classifier (python predictions)
def supp(job, out=f"{BENCH}/supp", shard="0/1"):
    from sklearn.metrics import roc_auc_score, average_precision_score
    data = Data(f"{BENCH}/runs/full", f"{BENCH}/runs/coverage", f"{BENCH}/runs/purity"); D = data.D; os.makedirs(out, exist_ok=True)
    tmp = f"{out}/.{job}.{shard.replace('/', '_')}.pass.vcf"

    def fit(tr, fitm=fit_models):
        m = train(data.training(tuple(tr)), fitm); return m["models"], m["cuts"], m["cut"]

    def outcome(ms, cuts, c0, p):
        r = {}
        for (q, k), (L, X, K, y, keep) in D.items():
            if q != p: continue
            pr = predict(ms, X); m = int(round(X[0, 56])); cut = cuts.get(m, c0); d = mark(K, pr)
            with open(tmp, "w") as o:
                for l, x, u in zip(L, pr, d):
                    if x >= cut and not u: g = l.split("\t"); g[6] = "PASS"; o.write("\t".join(g))
            n, tpc, tpt, P, R, F = data.score(tmp, p); ok = ~d
            auc = roc_auc_score(y[ok], pr[ok]) if 0 < y[ok].sum() < ok.sum() else float("nan")
            ap = average_precision_score(y[ok], pr[ok]) if y[ok].sum() > 0 else float("nan")
            r[k] = dict(calls=n, P=P, R=R, F1=F, roc_auc=auc, pr_auc=ap, cut=cut, minimum=m)
        return r

    def evaluate(name, splits, fitm=fit_models):
        rows = []
        for tr, te in splits:
            ms, cuts, c0 = fit(tr, fitm)
            for p in te:
                for k, x in outcome(ms, cuts, c0, p).items(): rows.append(dict(variant=name, train=",".join(tr), pair=p, setting=k, **x))
        return rows

    def write(rows, f):
        keys = list(rows[0])
        with open(f, "w") as o:
            o.write("\t".join(keys) + "\n")
            for r in rows: o.write("\t".join(str(round(r[k], 4)) if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")

    LOPO = [([q for q in CASTLE if q != p], [p]) for p in CASTLE]
    if job == "auc": write(evaluate("D", LOPO), f"{out}/heldout_auc.tsv")
    if job == "splits": write(evaluate("D", [(list(c), [p for p in CASTLE if p not in c]) for c in itertools.combinations(CASTLE, 3)]), f"{out}/splits_3_2.tsv")
    if job == "ablation":
        rows = []
        for key, vals in [("max_iter", [100, 600]), ("learning_rate", [0.02, 0.1]), ("max_leaf_nodes", [7, 31]), ("min_samples_leaf", [10, 160]), ("l2_regularization", [0.0, 5.0])]:
            for val in vals: rows += evaluate(f"{key}={val}", LOPO, lambda X, y, h=dict(HGB, **{key: val}): fit_models(X, y, h))
        rows += evaluate("seeds=1", LOPO, lambda X, y: [HistGradientBoostingClassifier(**HGB, random_state=0).fit(X, y)])
        write(rows, f"{out}/ablation.tsv")
    if job == "grid":
        # 3^5 parameter sets, every set run like the paper's classifier with each CASTLE pair held out; the threshold and the
        # inner F1 come from nested cross-validation inside the training pairs
        AX = [("learning_rate", [0.02, 0.05, 0.1]), ("max_iter", [150, 300, 600]), ("max_leaf_nodes", [7, 15, 31]), ("min_samples_leaf", [10, 40, 160]), ("l2_regularization", [0.0, 1.0, 5.0])]
        i0, n0 = map(int, shard.split("/")); f = f"{out}/grid.{i0}.tsv"; seen = {}
        for g in glob.glob(f"{out}/grid.*.tsv"):
            for l in open(g).read().split("\n")[1:]:
                if l: seen[l.split("\t")[0]] = seen.get(l.split("\t")[0], 0) + 1
        done = {v for v, n in seen.items() if n >= len(D)}
        o = open(f, "a")
        if not o.tell(): o.write("variant\tpair\tsetting\tcalls\tP\tR\tF1\troc_auc\tpr_auc\tcut\tminimum\tinner_f1\n")
        XY = lambda tr: (np.vstack([x for _, x, _ in data.training(tr)]), np.concatenate([y for _, _, y in data.training(tr)]))
        todo = [c for c in itertools.product(*[v for _, v in AX]) if ",".join(f"{a}={v}" for (a, _), v in zip(AX, c)) not in done]
        for c in todo[i0::n0]:
            h = dict(zip([a for a, _ in AX], c)); name = ",".join(f"{a}={v}" for a, v in h.items())
            M = {tr: fit_models(*XY(tr), dict(HGB, **h)) for r in (4, 3) for tr in itertools.combinations(CASTLE, r)}
            for p in CASTLE:
                tr = tuple(q for q in CASTLE if q != p); oof, ys, gs = [], [], []
                for q in tr:
                    mi = M[tuple(x for x in tr if x != q)]
                    for k in [k for k in D if k[0] == q and k[1] != "100p"]:
                        Xk, yk = D[k][1][D[k][4]], D[k][3][D[k][4]]; oof.append(predict(mi, Xk)); ys.append(yk); gs.append(np.round(Xk[:, 56]))
                oo, y, g = np.concatenate(oof), np.concatenate(ys), np.concatenate(gs)
                cuts = {int(m): float(max(CUTS, key=lambda t: f1(oo[g == m], y[g == m], t))) for m in set(g)}
                c0 = float(max(CUTS, key=lambda t: f1(oo, y, t)))
                inner = float(np.mean([f1(a, b, cuts.get(int(m[0]), c0)) for a, b, m in zip(oof, ys, gs)]))
                for k, x in outcome(M[tr], cuts, c0, p).items():
                    o.write("\t".join([name, p, k, str(x["calls"])] + [f"{x[v]:.4f}" for v in ("P", "R", "F1", "roc_auc", "pr_auc", "cut")] + [str(x["minimum"]), f"{inner:.4f}"]) + "\n")
            o.flush(); print(job, name, flush=True)
    if job == "firstpass":
        # the classifier trained with the records marked DUPL by janus call excluded (the paper) and kept
        orig = dict(D); rows = []
        for name in ("excluded", "kept"):
            for k, v in orig.items(): D[k] = v if name == "excluded" else v[:4] + (np.ones(len(v[3]), bool),)
            for p in CASTLE:
                ms, cuts, c0 = fit([q for q in CASTLE if q != p])
                for k, x in outcome(ms, cuts, c0, p).items(): rows.append(dict(variant=name, pair=p, setting=k, **x))
        D.update(orig); write(rows, f"{out}/firstpass.tsv")
    if job == "features":
        # permutation importance of the feature groups on every pair, with the model of the pair
        GROUPS = [("Tumor reads, single anchor pair", [1, 2, 3, 4, 5, 35]), ("Normal reads, single anchor pair", [6, 7, 8, 9, 10, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34]),
                  ("Tumor reads, anchor grid and paths", [39, 40, 49]), ("Normal reads, anchor grid and paths", [38, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52]),
                  ("Graph context", [12, 13, 36, 37, 53, 54, 55]), ("Alignment quality", [14, 15]), ("Reference context and depth", [16, 17, 23]),
                  ("Type and size", [18, 19, 20, 21, 22]), ("Minimum and count-stage context", [56, 57, 58, 59, 60])]
        F = json.load(open(f"{BENCH}/verdict/folds.json")); res = {}
        for p in CASTLE + ["hg008", "colo829"]:
            tr = [q for q in CASTLE if q != p]; ms = fit_models(np.vstack([x for _, x, _ in data.training(tr)]), np.concatenate([y for _, _, y in data.training(tr)]))
            L, X, K, y, keep = D[(p, "full")]; cut = F[p]["cut"]

            def f1of(Z):
                pr = predict(ms, Z); d = mark(K, pr)
                return float(f1(pr[~d], y[~d], cut)), pr
            base, pr = f1of(X)
            assert np.abs(pr - np.array(F[p]["prob"])).max() < 1e-3, p
            res[p] = dict(base=base, cut=cut, drop={})
            for name, ix in GROUPS:
                v = []
                for r in range(5):
                    Z = X.copy(); Z[:, ix] = X[np.random.RandomState(r).permutation(len(X))][:, ix]; v.append(base - f1of(Z)[0])
                res[p]["drop"][name] = v
            print(p, round(base, 4), {n: round(float(np.mean(v)), 4) for n, v in res[p]["drop"].items()}, flush=True)
        json.dump(dict(groups=GROUPS, pairs=res), open(f"{out}/features.json", "w"), indent=1)
    print("done", job)


if __name__ == "__main__":
    a = sys.argv[1:]; kw = {a[i][2:]: a[i + 1] for i in range(1, len(a) - 1) if a[i].startswith("--")}
    if a[0] == "verdict": verdict(**kw)
    elif a[0] == "supp": supp(a[1], **kw)
    else: sys.exit(__doc__)
