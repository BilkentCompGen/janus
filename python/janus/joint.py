"""Joint runs: acceptance of the records of several tumor samples, carriers from spanning reads, and the clade pass with the tree of the tumors."""
import argparse, collections, pickle
import numpy as np
from . import features as F
from . import classify as C

RAISE, FLOOR, MU = 0.2, 0.1, 0.5
NEAR = 500  # records of different samples describe the same SV if both breakpoints lie within this distance


def neighbors(K, du):
    """Records of all samples within NEAR of each non-duplicate record, and the SV groups (transitive closure)."""
    key = lambda q, da=0, db=0: (q[0], q[1], q[2] // NEAR + da, q[3] // NEAR + db)
    grid = collections.defaultdict(list)
    for i in np.where(~du)[0]: grid[key(K[i])].append(i)
    NB = {i: [j for da in (-1, 0, 1) for db in (-1, 0, 1) for j in grid.get(key(K[i], da, db), ()) if abs(K[j][2] - K[i][2]) <= NEAR and abs(K[j][3] - K[i][3]) <= NEAR] for i in np.where(~du)[0]}
    root = list(range(len(K)))

    def find(i):
        while root[i] != i: root[i] = root[root[i]]; i = root[i]
        return i
    for i in np.where(~du)[0]:
        for j in NB.get(i, ()): root[find(i)] = find(j)
    return NB, np.array([find(i) for i in range(len(K))])


def accept(pr, cut, du, own, NB, SM, SP, T, raise_=RAISE, floor=FLOOR):
    """First pass. A record is accepted if p >= cut + raise_, or if p >= floor and another tumor has a record of the SV at or above its cut.
    A tumor without a record carries an accepted SV if at least three tumors have records, one spanning read matches, and matching reads are
    at least a quarter of the spanning reads. Returns the accepted mask, the records counted for each tumor, and the carried records among them."""
    N = len(pr)
    conf = np.array([len({own[j] for j in NB.get(i, ()) if own[j] != own[i] and pr[j] >= cut[j]}) for i in range(N)])
    nown = np.array([len({own[j] for j in NB.get(i, ())} | {own[i]}) for i in range(N)])
    anyown = [set() for _ in range(T)]
    for i in np.where(~du)[0]: anyown[own[i]].update(NB[i])
    acc = ~du & ((pr >= cut + raise_) | ((pr >= floor) & (conf >= 1)))
    lists, adds = [], []
    for k in range(T):
        mine = list(np.where(acc & (own == k))[0]); have = set(mine)
        for i in mine: have.update(NB[i])
        add, seen = [], set()
        ok = acc & (own != k) & (SM[:, k] >= 1) & (0.25 * SP[:, k] <= SM[:, k]) & (nown >= 3)
        for i in np.where(ok)[0]:
            if i in anyown[k] or i in have or i in seen: continue
            add.append(i); seen.update(NB[i])
        lists.append(mine + add); adds.append(add)
    return acc, lists, adds


def nj_sides(cnt, T):
    """Neighbor joining on the Manhattan distances between the presence/absence vectors of the carrier sets; both sides of every edge and the full set."""
    ALL = frozenset(range(T)); sets = list(cnt.items()); V = np.zeros((T, len(sets))); w = np.array([n for _, n in sets], float)
    for c, (q, _) in enumerate(sets):
        for k in q: V[k, c] = 1
    D = {(i, j): float((np.abs(V[i] - V[j]) * w).sum()) for i in range(T) for j in range(T)}
    nodes = list(range(T)); children = {}; nxt = T
    while len(nodes) > 2:
        m = len(nodes); r = {i: sum(D[(i, j)] for j in nodes if j != i) for i in nodes}; best = None
        for a in range(m):
            for b in range(a + 1, m):
                i, j = nodes[a], nodes[b]; q = (m - 2) * D[(i, j)] - r[i] - r[j]
                if best is None or q < best[0]: best = (q, i, j)
        _, i, j = best; u = nxt; nxt += 1; children[u] = (i, j)
        for k in nodes:
            if k not in (i, j): D[(u, k)] = D[(k, u)] = 0.5 * (D[(i, k)] + D[(j, k)] - D[(i, j)])
        D[(u, u)] = 0.0; nodes = [k for k in nodes if k not in (i, j)] + [u]
    children[nxt] = tuple(nodes)
    leaves = lambda u: frozenset([u]) if u < T else frozenset().union(*(leaves(c) for c in children[u]))
    sides = {ALL}
    for u in children:
        for c in children[u]:
            q = leaves(c); sides.add(q); sides.add(ALL - q)
    return sides


def clade_pass(acc, lists, pr, cut, du, own, grp, T, mu=MU, cut0=None):
    """Second pass. The carriers of each SV group are re-decided with the tree built from the first-pass carrier sets: every side of every edge
    is a candidate carrier set, scored by the evidence of its members (logit p - logit cut for a tumor with a record, -mu otherwise) plus a
    prior from the number of SVs whose carriers are exactly that set. The members of the best candidate without an accepted record become
    carriers of the group's best accepted record; carriers outside it are dropped."""
    members = collections.defaultdict(list)
    for i in np.where(~du)[0]: members[grp[i]].append(i)
    GL = sorted(members); gi = {g: n for n, g in enumerate(GL)}; NG = len(GL); gidx = np.array([gi.get(g, -1) for g in grp])
    R = np.full((NG, T), np.nan)
    for g in GL:
        for j in members[g]:
            if not (R[gi[g], own[j]] >= pr[j]): R[gi[g], own[j]] = pr[j]
    carr = collections.defaultdict(set)
    for k in range(T):
        for i in lists[k]: carr[grp[i]].add(k)
    cnt = collections.Counter(frozenset(v) for v in carr.values())
    order = sorted(nj_sides(cnt, T), key=lambda c: (len(c), sorted(c)))
    K = np.array([[1.0 if k in c else 0.0 for k in range(T)] for c in order]); nK = np.array([cnt.get(c, 0) for c in order], float)
    prior = np.log((nK + 1) / (sum(cnt.values()) + len(order)) * len(order))
    tau = np.array([cut[own == k][0] if (own == k).any() else (np.median(cut) if cut0 is None else cut0) for k in range(T)])
    logit = lambda x: np.log(x / (1 - x))
    E = np.where(~np.isnan(R), logit(np.clip(np.nan_to_num(R, nan=0.5), 1e-6, 1 - 1e-6)) - logit(tau)[None, :], -mu)
    S = E @ K.T + prior[None, :]; bi = S.argmax(1); Kstar = K[bi] * (S[np.arange(NG), bi] > 0)[:, None]
    accg = collections.defaultdict(set); best = {}
    for i in np.where(acc)[0]:
        g = gidx[i]; accg[g].add(int(own[i]))
        if g not in best or pr[i] > pr[best[g]]: best[g] = i
    lists2, adds2 = [], []
    for k in range(T):
        add = [best[g] for g in np.where(Kstar[:, k] > 0)[0] if g in accg and k not in accg[g]]
        lists2.append(list(np.where(acc & (own == k))[0]) + add); adds2.append(add)
    return lists2, adds2


def run(model, X, K, own, SM, SP, hists, minimums, mu=MU, clade=True):
    """Classify the records of each tumor with its own counts, then the two passes. Returns probabilities, duplicate marks, thresholds,
    the accepted mask, and the carried records per tumor."""
    T = len(hists); N = len(X); pr = np.zeros(N); cut = np.zeros(N); du = np.zeros(N, bool)
    for k in range(T):
        ix = np.where(own == k)[0]
        if not len(ix): continue
        pr[ix], du[ix], c = C.classify(model, X[ix], [K[i] for i in ix], hists[k], minimums[k]); cut[ix] = c
    NB, grp = neighbors(K, du)
    acc, lists, adds = accept(pr, cut, du, own, NB, SM, SP, T)
    if clade: lists, adds = clade_pass(acc, lists, pr, cut, du, own, grp, T, mu, model["cut"])
    return pr, du, cut, acc, lists, adds


def write(out, H, L, pr, du, cut, acc, adds, names):
    carriers = collections.defaultdict(list)
    for k, add in enumerate(adds):
        for i in add: carriers[i].append(names[k])
    with open(out, "w") as o:
        o.write("".join(H))
        for i, (l, r, u) in enumerate(zip(L, pr, du)):
            g = l.split("\t"); g[6] = "DUPL" if u else "PASS" if acc[i] else "LOW"; g[7] += ";VP=%.4f" % r
            if i in carriers: g[7] += ";CARRIERS=" + ",".join(carriers[i])
            o.write("\t".join(g))


def main():
    a = argparse.ArgumentParser(description="Accept the records of a joint janus VCF and report the carriers of each SV among the tumor samples.")
    a.add_argument("-i", required=True, help="calls.vcf of a joint janus call"); a.add_argument("-o", required=True, help="output VCF")
    a.add_argument("-m", required=True, help="model (models/janus.pkl)")
    a.add_argument("-H", required=True, help="t.cnt.hist files of the tumor samples, in the order of the VCF columns, comma-separated")
    a.add_argument("-s", required=True, help="minimum read support of each tumor sample, comma-separated")
    a.add_argument("--mu", type=float, default=MU, help=f"penalty of a tumor without a record in the clade pass (default {MU})")
    a.add_argument("--no-clade", action="store_true", help="stop after the first pass")
    a = a.parse_args()
    H, L, X, K, own, SM, SP = F.records(a.i); model = pickle.load(open(a.m, "rb")); hists = a.H.split(","); mins = [int(x) for x in a.s.split(",")]
    names = H[-1].rstrip("\n").split("\t")[10:]
    pr, du, cut, acc, lists, adds = run(model, X, K, own, SM, SP, hists, mins, a.mu, not a.no_clade)
    write(a.o, H, L, pr, du, cut, acc, adds, names)
    print(f"{a.o}: {int(acc.sum())} PASS of {len(L)} records, {sum(len(x) for x in adds)} carriers")


if __name__ == "__main__":
    main()
