"""Records of a janus VCF and the classifier features derived from them."""
import re
import numpy as np

ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
DUPT = 500  # two records of one sample are duplicates if both breakpoints lie within this distance


def hist(path):
    return {l.split()[0]: np.array(l.split()[1:], float) for l in open(path)}


def pct(h, f, lo):
    c = np.cumsum(h[lo:]) / max(h[lo:].sum(), 1)
    return float(lo + np.searchsorted(c, f))


def context(path, mt):
    """Depth and purity statistics of a count histogram: median shared-core counts of tumor and normal and the four context features."""
    h = hist(path)
    mn = int(np.nonzero(h["shared_normal"])[0][0])
    tm = pct(h["shared_tumor"], 0.5, mt); nm = pct(h["shared_normal"], 0.5, mn)
    return tm, nm, [np.log2(tm), np.log2(nm), pct(h["tumor_only"], 0.05, mt) / tm, pct(h["tumor_only"], 0.95, mt) / tm]


def minimum(X):
    """Minimum read support of the run, recovered from the smallest tumor read count in the records."""
    return int(round(2 ** X[:, 1].min() - 1))


def features(X, path, mt):
    """Input features of the classifier: the record features with the log2 counts normalised by the sample medians, the minimum read support, and the context features."""
    tm, nm, c4 = context(path, mt); Z = X.copy()
    for j in (1, 3, 39): Z[:, j] -= np.log2(tm + 1)
    for j in (26, 27, 33, 34, 38, 41, 47, 51): Z[:, j] -= np.log2(nm + 1)
    return np.hstack([Z, np.full((len(X), 1), mt, float), np.tile(c4, (len(X), 1))])


def records(path):
    """Parse a janus VCF: header lines, record lines, raw features, breakpoint keys, tumor index of each record, and the SAME and SPAN read counts per tumor."""
    H, L, X, K, own, SM, SP = [], [], [], [], [], [], []
    for l in open(path):
        if l[0] == "#": H.append(l); continue
        f = l.rstrip("\n").split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
        K.append((f[0], m.group(1) if m else f[0], int(f[1]), int(m.group(2)) if m else int(f[1]) + int(info["SVLEN"]), info["SVTYPE"], info["O"]))
        X.append([float(v) for v in info["F"].split(",")]); L.append(l)
        gt = [x.split(":")[0] for x in f[10:]]; own.append(gt.index("0/1") if "0/1" in gt else 0)
        SM.append([int(x) for x in info["SAME"].split(",")][1:]); SP.append([int(x) for x in info["SPAN"].split(",")][1:])
    return H, L, np.array(X), K, np.array(own), np.array(SM), np.array(SP)


def duplicates(K, pr):
    """Mark the record with the lower probability in each pair of duplicate records of one sample."""
    ix = sorted(range(len(K)), key=lambda i: (K[i][:4], i)); d = np.zeros(len(K), bool)
    for a in range(len(ix)):
        x = K[ix[a]]
        for b in ix[a + 1:]:
            y = K[b]
            if y[:2] != x[:2] or y[2] > x[2] + DUPT: break
            if abs(y[3] - x[3]) <= DUPT and (x[4] != y[4] or x[5] == y[5]): d[ix[a] if pr[ix[a]] < pr[b] else b] = True
    return d
