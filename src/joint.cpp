#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include "cmd.h"
#include "io.h"
#include "model.h"

// Joint runs: the records of several tumor samples of one patient. A record is accepted when its probability clears the
// threshold with a margin, or when it is low but another tumor has a record of the same SV; a tumor without a record
// carries an accepted SV when its spanning reads match; then the carriers of every SV are re-decided with the tree of
// the tumors built from the first-pass carrier sets (the clade pass).

static const int64_t NEAR = 500;

typedef uint64_t Set;

static std::vector<int> members(Set s) {
    std::vector<int> v;
    for (int k = 0; s; k++, s >>= 1)
        if (s & 1) v.push_back(k);
    return v;
}

static bool by_set(Set a, Set b) {
    int ca = __builtin_popcountll(a), cb = __builtin_popcountll(b);
    return ca != cb ? ca < cb : members(a) < members(b);
}

class Joint {
public:
    Joint(const Vcf& v, int T, double raise_, double floor_, double mu) : v(v), r(v.r), N(v.r.size()), T(T), RAISE(raise_), FLOOR(floor_), MU(mu), nb(N), root(N), pr(N, 0), cut(N, 0), du(N, 0), acc(N, 0), lists(T), adds(T) {}
    void neighbors();
    void accept();
    void clade(double cut0);

private:
    const Vcf& v;
    const std::vector<Sv>& r;
    size_t N;
    int T;
    double RAISE, FLOOR, MU;
    std::vector<std::vector<size_t>> nb;
    std::vector<size_t> root;

public:
    std::vector<double> pr, cut;
    std::vector<uint8_t> du, acc;
    std::vector<std::vector<size_t>> lists, adds;

private:
    size_t find(size_t i);
    std::set<Set> sides(const std::map<Set, int>& cnt) const;
};

size_t Joint::find(size_t i) {
    while (root[i] != i) { root[i] = root[root[i]]; i = root[i]; }
    return i;
}

void Joint::neighbors() {
    std::map<std::tuple<std::string, std::string, int64_t, int64_t>, std::vector<size_t>> grid;
    for (size_t i = 0; i < N; i++)
        if (!du[i]) grid[{r[i].c1, r[i].c2, r[i].p1 / NEAR, r[i].p2 / NEAR}].push_back(i);
    for (size_t i = 0; i < N; i++) {
        if (du[i]) continue;
        for (int da = -1; da <= 1; da++)
            for (int db = -1; db <= 1; db++) {
                std::map<std::tuple<std::string, std::string, int64_t, int64_t>, std::vector<size_t>>::const_iterator it = grid.find({r[i].c1, r[i].c2, r[i].p1 / NEAR + da, r[i].p2 / NEAR + db});
                if (it == grid.end()) continue;
                for (size_t j : it->second)
                    if (std::llabs(r[j].p1 - r[i].p1) <= NEAR && std::llabs(r[j].p2 - r[i].p2) <= NEAR) nb[i].push_back(j);
            }
    }
    for (size_t i = 0; i < N; i++) root[i] = i;
    for (size_t i = 0; i < N; i++)
        if (!du[i]) for (size_t j : nb[i]) root[find(i)] = find(j);
}

void Joint::accept() {
    std::vector<int> conf(N, 0), nown(N, 0);
    std::vector<std::set<size_t>> anyown(T);
    for (size_t i = 0; i < N; i++) {
        std::set<int> c, o = {r[i].own};
        for (size_t j : nb[i]) {
            o.insert(r[j].own);
            if (r[j].own != r[i].own && pr[j] >= cut[j]) c.insert(r[j].own);
        }
        conf[i] = c.size();
        nown[i] = o.size();
        if (!du[i]) anyown[r[i].own].insert(nb[i].begin(), nb[i].end());
    }
    for (size_t i = 0; i < N; i++) acc[i] = !du[i] && (pr[i] >= cut[i] + RAISE || (pr[i] >= FLOOR && conf[i] >= 1));
    for (int k = 0; k < T; k++) {
        std::vector<size_t> mine;
        std::set<size_t> have, seen;
        for (size_t i = 0; i < N; i++)
            if (acc[i] && r[i].own == k) { mine.push_back(i); have.insert(i); have.insert(nb[i].begin(), nb[i].end()); }
        std::vector<size_t> add;
        for (size_t i = 0; i < N; i++) {
            if (!(acc[i] && r[i].own != k && r[i].same[k] >= 1 && 0.25 * r[i].span[k] <= r[i].same[k] && nown[i] >= 3)) continue;
            if (anyown[k].count(i) || have.count(i) || seen.count(i)) continue;
            add.push_back(i);
            seen.insert(nb[i].begin(), nb[i].end());
        }
        lists[k] = mine;
        lists[k].insert(lists[k].end(), add.begin(), add.end());
        adds[k] = add;
    }
}

std::set<Set> Joint::sides(const std::map<Set, int>& cnt) const {
    // neighbor joining on the Manhattan distances between the presence/absence vectors of the carrier sets
    const Set ALL = T == 64 ? ~0ull : (1ull << T) - 1;
    std::map<std::pair<int, int>, double> D;
    for (int i = 0; i < T; i++)
        for (int j = 0; j < T; j++) {
            double d = 0;
            for (const std::pair<const Set, int>& s : cnt) d += (double)(((s.first >> i) ^ (s.first >> j)) & 1) * s.second;
            D[{i, j}] = d;
        }
    std::vector<int> nodes(T);
    for (int i = 0; i < T; i++) nodes[i] = i;
    std::map<int, std::vector<int>> children;
    int nxt = T;
    while (nodes.size() > 2) {
        int m = nodes.size();
        std::map<int, double> rs;
        for (int i : nodes) {
            double s = 0;
            for (int j : nodes) if (j != i) s += D[{i, j}];
            rs[i] = s;
        }
        double best = 0;
        int bi = -1, bj = -1;
        for (int a = 0; a < m; a++)
            for (int b = a + 1; b < m; b++) {
                int i = nodes[a], j = nodes[b];
                double q = (m - 2) * D[{i, j}] - rs[i] - rs[j];
                if (bi < 0 || q < best) { best = q; bi = i; bj = j; }
            }
        int u = nxt++;
        children[u] = {bi, bj};
        for (int k : nodes)
            if (k != bi && k != bj) D[{u, k}] = D[{k, u}] = 0.5 * (D[{bi, k}] + D[{bj, k}] - D[{bi, bj}]);
        D[{u, u}] = 0;
        std::vector<int> rest;
        for (int k : nodes) if (k != bi && k != bj) rest.push_back(k);
        rest.push_back(u);
        nodes = rest;
    }
    children[nxt] = nodes;
    std::map<int, Set> leaves;
    for (int u = 0; u < T; u++) leaves[u] = 1ull << u;
    for (int u = T; u <= nxt; u++) {
        Set s = 0;
        for (int c : children[u]) s |= leaves[c];
        leaves[u] = s;
    }
    std::set<Set> out = {ALL};
    for (const std::pair<const int, std::vector<int>>& u : children)
        for (int c : u.second) { out.insert(leaves[c]); out.insert(ALL & ~leaves[c]); }
    return out;
}

void Joint::clade(double cut0) {
    std::map<size_t, std::vector<size_t>> groups;
    for (size_t i = 0; i < N; i++)
        if (!du[i]) groups[find(i)].push_back(i);
    std::vector<size_t> gl;
    std::map<size_t, size_t> gi;
    for (const std::pair<const size_t, std::vector<size_t>>& g : groups) { gi[g.first] = gl.size(); gl.push_back(g.first); }
    size_t NG = gl.size();
    std::vector<std::vector<double>> R(NG, std::vector<double>(T, NAN));
    for (size_t g = 0; g < NG; g++)
        for (size_t j : groups[gl[g]])
            if (!(R[g][r[j].own] >= pr[j])) R[g][r[j].own] = pr[j];
    std::map<size_t, Set> carr;
    for (int k = 0; k < T; k++)
        for (size_t i : lists[k]) carr[find(i)] |= 1ull << k;
    std::map<Set, int> cnt;
    int total = 0;
    for (const std::pair<const size_t, Set>& c : carr) { cnt[c.second]++; total++; }
    std::set<Set> sd = sides(cnt);
    std::vector<Set> order(sd.begin(), sd.end());
    std::sort(order.begin(), order.end(), by_set);
    std::vector<double> prior(order.size()), tau(T);
    for (size_t c = 0; c < order.size(); c++) {
        std::map<Set, int>::const_iterator it = cnt.find(order[c]);
        prior[c] = std::log((double)((it == cnt.end() ? 0 : it->second) + 1) / (total + (double)order.size()) * order.size());
    }
    for (int k = 0; k < T; k++) {
        tau[k] = cut0;
        for (size_t i = 0; i < N; i++)
            if (r[i].own == k) { tau[k] = cut[i]; break; }
    }
    auto logit = [](double x) { return std::log(x / (1 - x)); };
    std::vector<Set> kstar(NG, 0);
    for (size_t g = 0; g < NG; g++) {
        std::vector<double> e(T);
        for (int k = 0; k < T; k++) e[k] = std::isnan(R[g][k]) ? -MU : logit(std::min(std::max(R[g][k], 1e-6), 1 - 1e-6)) - logit(tau[k]);
        double best = 0;
        size_t bi = 0;
        for (size_t c = 0; c < order.size(); c++) {
            double s = 0;
            for (int k = 0; k < T; k++) s += e[k] * (double)(order[c] >> k & 1);
            s += prior[c];
            if (!c || s > best) { best = s; bi = c; }
        }
        if (best > 0) kstar[g] = order[bi];
    }
    std::vector<Set> accg(NG, 0);
    std::vector<int64_t> bestrec(NG, -1);
    for (size_t i = 0; i < N; i++) {
        if (!acc[i]) continue;
        size_t g = gi[find(i)];
        accg[g] |= 1ull << r[i].own;
        if (bestrec[g] < 0 || pr[i] > pr[bestrec[g]]) bestrec[g] = i;
    }
    for (int k = 0; k < T; k++) {
        std::vector<size_t> add;
        for (size_t g = 0; g < NG; g++)
            if ((kstar[g] >> k & 1) && accg[g] && !(accg[g] >> k & 1)) add.push_back(bestrec[g]);
        lists[k].clear();
        for (size_t i = 0; i < N; i++)
            if (acc[i] && r[i].own == k) lists[k].push_back(i);
        lists[k].insert(lists[k].end(), add.begin(), add.end());
        adds[k] = add;
    }
}

int cmd_joint(const Args& a) {
    Vcf v(a.need(a.calls, "--calls"));
    Model m(a.need(a.model, "--model"));
    std::vector<std::string> hists = split(a.need(a.hist, "--hist")), mins = split(a.need(a.minimum, "--minimum")), names = v.tumors();
    int T = hists.size();
    if ((int)mins.size() != T || (int)names.size() != T) { std::fprintf(stderr, "janus joint: %zu hist files, %zu minimums and %zu tumor columns\n", hists.size(), mins.size(), names.size()); return 2; }
    Joint j(v, T, a.raise_, a.floor_, a.mu);
    for (int k = 0; k < T; k++) {
        std::vector<size_t> ix;
        for (size_t i = 0; i < v.r.size(); i++)
            if (v.r[i].own == k) ix.push_back(i);
        if (ix.empty()) continue;
        int mt = std::atoi(mins[k].c_str());
        classify(m, v.r, ix, hists[k].c_str(), mt, j.pr);
        std::vector<Sv> sub;
        std::vector<double> sp;
        for (size_t i : ix) { sub.push_back(v.r[i]); sp.push_back(j.pr[i]); }
        std::vector<uint8_t> d = duplicates(sub, sp);
        for (size_t q = 0; q < ix.size(); q++) { j.du[ix[q]] = d[q]; j.cut[ix[q]] = m.threshold(mt); }
    }
    j.neighbors();
    j.accept();
    if (!a.noclade) j.clade(m.cut);
    std::vector<std::string> extra(v.r.size());
    size_t carriers = 0;
    for (int k = 0; k < T; k++)
        for (size_t i : j.adds[k]) { extra[i] += (extra[i].empty() ? ";CARRIERS=" : ",") + names[k]; carriers++; }
    v.write(a.need(a.out, "--out"), j.pr, j.du, j.acc, extra);
    std::fprintf(stderr, "[joint] %zu PASS of %zu records, %zu carriers\n", (size_t)std::count(j.acc.begin(), j.acc.end(), 1), v.r.size(), carriers);
    return 0;
}
