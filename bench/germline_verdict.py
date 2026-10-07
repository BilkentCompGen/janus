import bisect, collections, gzip, sys, numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier

CALLS, TRUTH, BED, TRF, OUT = sys.argv[1:6]
RD, PS = int(sys.argv[6]), float(sys.argv[7])
ODD = {str(i) for i in range(1, 23, 2)}
HGB = dict(max_iter=400, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=40, l2_regularization=1.0)
CUTS = np.arange(0.1, 0.9, 0.05)
SEEDS = int(__import__("os").environ.get("SEEDS", "3"))
GRID = [(k, w, r) for k in (1, 2) for w in (300, 500, 1000) for r in (0.5, 0.7)]


def intervals(path, gap=0):
    iv = collections.defaultdict(list)
    for l in open(path):
        c, a, b = l.split()[:3]
        iv[c].append((int(a), int(b)))
    for c in iv:
        iv[c].sort()
        m = []
        for a, b in iv[c]:
            if m and a <= m[-1][1] + gap: m[-1] = (m[-1][0], max(m[-1][1], b))
            else: m.append((a, b))
        iv[c] = m
    return iv, {c: [a for a, _ in iv[c]] for c in iv}


def inside(iv, st, c, p):
    i = bisect.bisect_right(st.get(c, []), p) - 1
    return i if i >= 0 and iv[c][i][1] >= p else -1


bed, bst = intervals(BED)
tr, tst = intervals(TRF, 10)
truth = collections.defaultdict(list)
for l in gzip.open(TRUTH, "rt"):
    if l[0] == "#": continue
    f = l.split("\t")
    d = len(f[4]) - len(f[3])
    if abs(d) >= 50: truth[f[0]].append((int(f[1]), "INS" if d > 0 else "DEL", abs(d), inside(bed, bst, f[0], int(f[1])) >= 0))
for c in truth: truth[c].sort()
tpos = {c: [x[0] for x in truth[c]] for c in truth}


def near(c, p, t, n):
    a, b = bisect.bisect_left(tpos.get(c, []), p - RD), bisect.bisect_right(tpos.get(c, []), p + RD)
    return [z for z in range(a, b) if truth[c][z][1] == t and n and PS <= truth[c][z][2] / n <= 1 / PS]


def repeat(c, p):
    i = inside(tr, tst, c, p)
    if i >= 0:
        a, b = tr[c][i]
        return [1, np.log2(b - a + 1), np.log2(min(p - a, b - p) + 1)]
    j = bisect.bisect_right(tst.get(c, []), p)
    gaps = ([p - tr[c][j - 1][1]] if j > 0 else []) + ([tr[c][j][0] - p] if c in tr and j < len(tr[c]) else [])
    return [0, 0, -np.log2(min(gaps, default=10 ** 6) + 1)]


lines, X, C, P, T, S = [], [], [], [], [], []
for l in open(CALLS):
    if l[0] == "#": continue
    f = l.split("\t")
    d = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
    t, n = d["SVTYPE"], int(d.get("SVLEN", 0))
    X.append([float(v) for v in d["F"].split(",")] + [f[6] == "DUPL", f[6] == "PASS", 1 - float(d["ERR"]), np.log2(n + 1)] + [t == k for k in ("INS", "DEL", "DUP", "BND", "INV")] + [np.log2(int(d.get("PIECES", 1)) + 1)] + repeat(f[0], int(f[1])))
    lines.append(f); C.append(f[0]); P.append(int(f[1])); T.append("INS" if t == "DUP" else t); S.append(n)
X = np.array(X, float)
by = collections.defaultdict(list)
for i, (c, p) in enumerate(zip(C, P)): by[c].append((p, i))
extra = np.zeros((len(lines), 7))
for c, v in by.items():
    v.sort()
    ps = [p for p, _ in v]
    for p, i in v:
        a, b = bisect.bisect_left(ps, p - 500), bisect.bisect_right(ps, p + 500)
        others = [v[k][1] for k in range(a, b) if v[k][1] != i]
        same = [j for j in others if T[j] == T[i]]
        sim = [j for j in same if S[i] and 0.7 <= S[j] / S[i] <= 1.43]
        nb = min(same, key=lambda j: abs(P[j] - p), default=None)
        pcs = [X[j][-4] for j in same]
        extra[i] = [len(others), len(same), len(sim), np.log2(abs(P[nb] - p) + 1) if nb is not None else 12, min(S[i], S[nb]) / max(S[i], S[nb], 1) if nb is not None else 0, (max(pcs) if pcs else 0) - X[i][-4], sum(1 for j in sim if X[j][-4] > X[i][-4])]
X = np.hstack([X, extra])
Y = np.array([bool(near(c, p, t, n)) for c, p, t, n in zip(C, P, T, S)], float)
INB = np.array([inside(bed, bst, c, p) >= 0 for c, p in zip(C, P)])
num = np.array([int(c[3:]) if c[3:].isdigit() else 23 for c in C])
odd = np.array([c[3:] in ODD for c in C])


def fit(m):
    return [HistGradientBoostingClassifier(**HGB, random_state=s).fit(X[m], Y[m]) for s in range(SEEDS)]


def predict(ms, m):
    return np.mean([x.predict_proba(X[m])[:, 1] for x in ms], 0)


def limit(idx, pr, keep, k, w, r):
    kept = collections.defaultdict(list)
    for i in idx[np.argsort(-pr[idx])]:
        if not keep[i]: continue
        if sum(1 for j in kept[(C[i], T[i])] if abs(P[j] - P[i]) <= w and S[i] and S[j] and min(S[i], S[j]) >= r * max(S[i], S[j])) >= k: keep[i] = False
        else: kept[(C[i], T[i])].append(i)


def score(idx, pr, keep, part):
    used, hit = set(), 0
    for i in idx[np.argsort(-pr[idx])]:
        if not keep[i] or not INB[i]: continue
        z = [(abs(truth[C[i]][z][0] - P[i]), z) for z in near(C[i], P[i], T[i], S[i]) if truth[C[i]][z][3] and (C[i], z) not in used]
        if z: used.add((C[i], min(z)[1])); hit += 1
    nt = sum(1 for c in truth if (c[3:] in ODD) == part for x in truth[c] if x[3])
    return 2 * hit / max((keep[idx] & INB[idx]).sum() + nt, 1)


pr, keep = np.zeros(len(Y)), np.zeros(len(Y), bool)
for test in (odd, ~odd):
    train = ~test
    oof = np.zeros(len(Y))
    q = train & (num % 4 < 2)
    for a, b in ((q, train & ~q), (train & ~q, q)): oof[b] = predict(fit(a), b)
    idx = np.where(train)[0]
    best, bestf = None, -1
    for k, w, r in GRID:
        for cut in CUTS:
            kk = oof >= cut
            limit(idx, oof, kk, k, w, r)
            f = score(idx, oof, kk, bool(odd[idx[0]]))
            if f > bestf: best, bestf = (cut, k, w, r), f
    pr[test] = predict(fit(train), test)
    cut, k, w, r = best
    keep[test] = pr[test] >= cut
    limit(np.where(test)[0], pr, keep, k, w, r)
    print(f"held-out {'odd' if test[0] else 'even'} chromosomes: cut {cut:.2f}, at most {k} per {w} bp at size ratio {r}", file=sys.stderr)
with open(OUT, "w") as o:
    for l in open(CALLS):
        if l[0] == "#": o.write(l)
    for f, k in zip(lines, keep):
        f[6] = "PASS" if k else "LOW"
        o.write("\t".join(f))
