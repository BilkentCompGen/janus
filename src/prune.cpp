#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <vector>
#include "cmd.h"
#include "files.h"
#include "io.h"
#include "keys.h"
#include "reads.h"

bool alike(const std::string& t, const std::string& n);

static const uint32_t X = 20;

static uint32_t edit(const char* x, uint32_t n, const char* y, uint32_t m) {
    if (n > m + X || m > n + X) return X;
    std::vector<uint32_t> p(m + 2, X), c(m + 2, X);
    for (uint32_t j = 0; j <= m && j < X; j++) p[j] = j;
    for (uint32_t i = 1; i <= n; i++) {
        uint32_t lo = i > X ? i - X : 1, hi = std::min(m, i + X), best = X;
        c[lo - 1] = lo == 1 ? std::min(i, X) : X;
        for (uint32_t j = lo; j <= hi; j++) best = std::min(best, c[j] = std::min(std::min(p[j] + 1, c[j - 1] + 1), std::min(p[j - 1] + (x[i - 1] != y[j - 1]), X)));
        if (hi < m) c[hi + 1] = X;
        if (best >= X) return X;
        p.swap(c);
    }
    return std::min(p[m], X);
}

bool alike(const std::string& t, const std::string& n) {
    uint32_t k = std::min(t.size(), n.size()), a = 0, b = 0;
    while (a < k && t[a] == n[a]) a++;
    while (b < k - a && t[t.size() - 1 - b] == n[n.size() - 1 - b]) b++;
    uint32_t dt = t.size() - a - b, dn = n.size() - a - b;
    return std::max(dt, dn) - std::min(dt, dn) < X && edit(t.data() + a, dt, n.data() + a, dn) < X;
}

struct Arm {
    std::vector<uint32_t> path;
    bool variant;
    uint64_t score;
    uint8_t lock;
    std::string bases;
};

class Bubbles : public Job {
public:
    std::vector<uint8_t> drop;
    std::vector<Arm> arms;
    Bubbles(const Graph& g, const Adjacency& a, int t);
    void find();
    void match(Chain& c, const uint8_t* rec, uint32_t m, uint32_t rc);
    void judge();

private:
    const Graph& g;
    const Adjacency& a;
    int t, phase;
    std::vector<uint64_t> mask;
    uint64_t full;
    std::vector<std::vector<uint32_t>> tips;
    std::vector<std::vector<Arm>> found;
    std::vector<uint32_t> first, index;
    std::vector<uint8_t> same;
    void task(int i) override;
    void tip(int i, bool back);
    uint32_t trace(uint32_t u, std::vector<uint32_t>& p, uint32_t lim) const;
    bool solid(uint32_t s, uint32_t m, uint32_t lim, std::vector<uint32_t>& q) const;
    bool variant(const std::vector<uint32_t>& p, const std::vector<uint32_t>& q) const;
    void bubbles(int i);
};

static const uint32_t ML = 4, WL = 40, TL = 12;

Bubbles::Bubbles(const Graph& g, const Adjacency& a, int t) : drop(g.K, 0), g(g), a(a), t(t), phase(0), mask(g.K, 0), full((1ull << g.S) - 1), tips(t), found(t) {
    for (uint64_t s = 0; s < g.S; s++)
        for (uint64_t v = 0; v < g.K; v++) mask[v] |= (uint64_t)(g.reads(s, v) > 0) << s;
}

void Bubbles::tip(int i, bool back) {
    const std::vector<uint32_t>& fo = back ? a.first_in : a.first_out, &fn = back ? a.source : a.target, &bo = back ? a.first_out : a.first_in;
    std::vector<uint32_t> ch;
    for (uint64_t v = share(g.K, i, t); v < share(g.K, i + 1, t); v++)
        for (uint32_t j = fo[v]; fo[v + 1] - fo[v] > 1 && j < fo[v + 1]; j++) {
            ch.clear();
            for (uint32_t u = fn[j]; ch.size() < TL && mask[u] != full && mask[u] == mask[fn[j]] && bo[u + 1] - bo[u] == 1; u = fn[fo[u]]) {
                ch.push_back(u);
                if (fo[u + 1] - fo[u] != 1) break;
            }
            if (!ch.empty() && fo[ch.back() + 1] == fo[ch.back()]) tips[i].insert(tips[i].end(), ch.begin(), ch.end());
        }
}

bool Bubbles::solid(uint32_t s, uint32_t m, uint32_t lim, std::vector<uint32_t>& q) const {
    std::unordered_map<uint32_t, uint32_t> from;
    std::vector<uint32_t> cur = {s}, nxt;
    for (uint32_t d = 0; d <= lim && !cur.empty() && from.size() < 256; d++) {
        nxt.clear();
        for (uint32_t u : cur)
            for (uint32_t j = a.first_out[u]; j < a.first_out[u + 1]; j++) {
                uint32_t v = a.target[j];
                if (v == m) {
                    for (uint32_t x = u; x != s; x = from[x]) q.push_back(x);
                    std::reverse(q.begin(), q.end());
                    return true;
                }
                if (mask[v] != full || from.count(v) || v == s) continue;
                from[v] = u;
                nxt.push_back(v);
            }
        cur.swap(nxt);
    }
    return false;
}

uint32_t Bubbles::trace(uint32_t u, std::vector<uint32_t>& p, uint32_t lim) const {
    for (uint32_t s = 0; s <= lim; s++) {
        p.push_back(u);
        if (a.outs(u) != 1) return NONE;
        uint32_t v = a.target[a.first_out[u]];
        if (a.ins(v) >= 2) return v;
        u = v;
    }
    return NONE;
}

bool Bubbles::variant(const std::vector<uint32_t>& p, const std::vector<uint32_t>& q) const {
    for (uint32_t v : p) if (mask[v] != mask[p[0]]) return false;
    for (uint32_t v : q) if (mask[v] == mask[p[0]]) return false;
    return mask[p[0]] != full;
}

static Arm arm(uint32_t s, const std::vector<uint32_t>& p, uint32_t m, bool v) {
    Arm r = {{s}, v, 0, 0, ""};
    r.path.insert(r.path.end(), p.begin(), p.end());
    r.path.push_back(m);
    return r;
}

void Bubbles::bubbles(int i) {
    for (uint64_t s = share(g.K, i, t); s < share(g.K, i + 1, t); s++) {
        uint32_t n = a.outs(s);
        if (n < 2) continue;
        std::vector<std::vector<uint32_t>> p(n);
        std::vector<uint32_t> m(n);
        std::vector<uint8_t> paired(n, 0);
        for (uint32_t k = 0; k < n; k++) {
            uint32_t u = a.target[a.first_out[s] + k];
            m[k] = trace(u, p[k], mask[u] != full ? WL : ML);
        }
        for (uint32_t k = 0; k < n; k++)
            for (uint32_t j = k + 1; j < n; j++) {
                if (m[k] == NONE || m[k] != m[j]) continue;
                bool vk = variant(p[k], p[j]), vj = variant(p[j], p[k]);
                if (!vk && !vj) continue;
                found[i].push_back(arm(s, p[k], m[k], vk));
                found[i].push_back(arm(s, p[j], m[k], vj));
                paired[k] = paired[j] = 1;
            }
        for (uint32_t k = 0; k < n; k++) {
            if (m[k] == NONE || paired[k] || mask[p[k][0]] == full) continue;
            std::vector<uint32_t> q;
            if (!solid(s, m[k], p[k].size() + 8, q) || !variant(p[k], q.empty() ? std::vector<uint32_t>(1, s) : q)) continue;
            found[i].push_back(arm(s, p[k], m[k], true));
            found[i].push_back(arm(s, q, m[k], false));
        }
    }
}

void Bubbles::task(int i) {
    if (phase < 2) tip(i, phase == 1);
    else if (phase == 2) bubbles(i);
    else
        for (uint64_t b = 2 * i; b < arms.size(); b += 2 * t) same[b / 2] = !arms[b].bases.empty() && !arms[b + 1].bases.empty() && alike(arms[b].bases, arms[b + 1].bases);
}

void Bubbles::find() {
    for (phase = 0; phase < 3; phase++) par(t);
    for (int i = 0; i < t; i++) {
        for (uint32_t v : tips[i]) drop[v] = 1;
        arms.insert(arms.end(), found[i].begin(), found[i].end());
    }
    tips.clear();
    found.clear();
    mask.clear();
    mask.shrink_to_fit();
    first.assign(g.K + 1, 0);
    index.resize(arms.size());
    for (const Arm& r : arms) first[r.path[0] + 1]++;
    for (uint32_t i = 0; i < g.K; i++) first[i + 1] += first[i];
    std::vector<uint32_t> p(first.begin(), first.end() - 1);
    for (uint32_t i = 0; i < arms.size(); i++) index[p[arms[i].path[0]]++] = i;
}

void Bubbles::match(Chain& c, const uint8_t* rec, uint32_t m, uint32_t rc) {
    for (uint32_t i = 0; i < m; i++)
        for (uint32_t q = first[c.id[i]]; q < first[c.id[i] + 1]; q++) {
            Arm* p = &arms[index[q]];
            uint32_t n = p->path.size(), j = 0;
            while (j < n && i + j < m && c.id[i + j] == p->path[j]) j++;
            if (j < n) continue;
            uint32_t st = c.start(i), en = c.end(i + n - 1);
            uint64_t sc = 1ull << 32 | (uint32_t)(c.r.head() ^ st ^ c.r.len);
            if (sc <= __atomic_load_n(&p->score, __ATOMIC_RELAXED)) continue;
            lock(p->lock);
            if (sc > p->score) { p->score = sc; p->bases = c.r.bases(rec, st, en, rc); }
            unlock(p->lock);
        }
}

void Bubbles::judge() {
    same.assign(arms.size() / 2, 0);
    phase = 3;
    par(t);
    for (uint32_t b = 0; b < arms.size(); b++)
        if (same[b / 2] && arms[b].variant)
            for (uint32_t j = 1; j + 1 < arms[b].path.size(); j++) drop[arms[b].path[j]] = 1;
}

class Prune : public Walker {
public:
    explicit Prune(const Args& a);
    int run();

private:
    const Args& a;
    Graph g;
    Adjacency adj;
    Bubbles bubbles;
    Combos combos;
    void visit(Chain& c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) override;
    bool odd(uint64_t i) const;
    void gather(std::vector<uint64_t>& rows) const;
    std::vector<uint8_t> seeds(const std::vector<uint8_t>& drop) const;
};

Prune::Prune(const Args& a) : Walker(a.threads), a(a), g(a.in), adj(g), bubbles(g, adj, a.threads), combos(g, adj) {}

std::vector<uint8_t> Prune::seeds(const std::vector<uint8_t>& drop) const {
    std::vector<uint8_t> sd(g.K, 0);
    for (uint64_t i = 0; i < g.E; i++) {
        uint32_t u = Graph::from(g.edges[i]), v = Graph::to(g.edges[i]);
        if (g.novel[i] && !drop[u] && !drop[v]) sd[u] = sd[v] = 1;
    }
    return sd;
}

void Prune::visit(Chain& c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) {
    bubbles.match(c, rec, m, rc);
    combos.tally(c, m, sample, rid);
}

bool Prune::odd(uint64_t i) const {
    bool on = false;
    for (uint64_t s = 1; s < combos.S; s++) on |= combos.n[i * (combos.S - 1) + s - 1] >= a.support;
    return on;
}

void Prune::gather(std::vector<uint64_t>& rows) const {
    const std::vector<uint8_t>& drop = bubbles.drop;
    const uint64_t S = combos.S;
    std::vector<std::array<uint64_t, 3>> w;
    for (uint64_t i = 0; i < combos.ZT; i++) {
        uint64_t k = combos.key[i];
        if (!k || !odd(i) || !combos.kept(k)) continue;
        uint64_t ei = g.edges[Combos::in(k)], eo = g.edges[Combos::out(k)];
        if (drop[Graph::from(ei)] || drop[Graph::to(ei)] || drop[Graph::from(eo)] || drop[Graph::to(eo)]) continue;
        for (uint64_t s = 1; s < S; s++)
            if (combos.n[i * (S - 1) + s - 1] >= a.support) w.push_back({Combos::central_read(s, combos.best[i * (S - 1) + s - 1]), k, i});
    }
    std::sort(w.begin(), w.end());
    for (const std::array<uint64_t, 3>& e : w) {
        rows.push_back(e[0]);
        rows.push_back(e[1]);
        rows.push_back(0);
        rows.insert(rows.end(), combos.n.begin() + e[2] * (S - 1), combos.n.begin() + (e[2] + 1) * (S - 1));
    }
}

int Prune::run() {
    bubbles.find();
    combos.table(combos.keep(seeds(std::vector<uint8_t>(g.K, 0))));
    index(g.hashes);
    std::vector<std::string> f = split(a.reads), tumors(f.begin() + 1, f.end());
    stream(std::vector<std::string>(1, f[0]), 0);
    for (;;) {
        stream(tumors, 1);
        if (!combos.full) break;
        std::fprintf(stderr, "[prune] the combination table filled up at %lu slots; counting the tumor reads again with %lu\n", (unsigned long)combos.ZT, (unsigned long)combos.ZT * 2);
        combos.grow();
    }
    bubbles.judge();
    std::vector<uint64_t> d, rows;
    for (uint64_t v = 0; v < g.K; v++)
        if (bubbles.drop[v]) d.push_back(g.hashes[v]);
    std::sort(d.begin(), d.end());
    std::string o(a.out);
    Output::save((o + ".drp").c_str(), Header('P', 1, d.size(), bubbles.arms.size() / 2), d.data(), d.size() * 8);
    combos.keep(seeds(bubbles.drop));
    gather(rows);
    Header h('S', 1, rows.size() / (2 + g.S), g.S);
    Output::save((o + ".sites").c_str(), h, rows.data(), rows.size() * 8);
    std::fprintf(stderr, "[prune] %zu dropped nodes, %lu sites\n", d.size(), (unsigned long)h.n);
    return 0;
}

int cmd_prune(const Args& a) {
    return Prune(a).run();
}
