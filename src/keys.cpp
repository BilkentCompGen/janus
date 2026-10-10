#include <algorithm>
#include "io.h"
#include "keys.h"

static const uint64_t ZF = 16;
static const uint64_t ZT0 = 8;

Combos::Combos(const Graph& g, const Adjacency& a) : S(g.S), Z(0), ZT(0), full(0), seg(g.K), g(g), a(a) {
    for (uint32_t u = 0; u + 1 < g.starts.size(); u++)
        for (uint32_t v = g.starts[u]; v < g.starts[u + 1]; v++) seg[v] = u;
}

uint64_t Combos::keep(const std::vector<uint8_t>& sd) {
    const std::vector<uint32_t>& st = g.starts;
    uint64_t np = 0, K = g.K;
    branch.assign(K, 0);
    for (uint64_t v = 0; v < K; v++) np += branch[v] = sd[v] && (a.ins(v) > 1 || a.outs(v) > 1);
    keeps.assign(st.size() - 1, 0);
    for (uint64_t v = 0; v < K; v++) keeps[seg[v]] |= sd[v] && st[seg[v] + 1] - st[seg[v]] > 1;
    for (uint8_t k : keeps) np += k;
    on.assign(K, 0);
    for (uint64_t v = 0; v < K; v++) on[v] = branch[v] || keeps[seg[v]];
    if (fork.empty()) {
        std::vector<uint32_t> has(K, 0), in(K, 0), out(K, 0);
        for (uint64_t v = 0; v < K; v++)
            for (uint64_t s = 1; s < S; s++) has[v] |= g.reads(0, v) ? ~0u : (uint32_t)(g.reads(s, v) > 0) << s;
        fork.assign(K, 0);
        for (uint64_t e : g.edges) {
            uint64_t u = Graph::from(e), v = Graph::to(e);
            uint32_t m = has[u] & has[v];
            fork[u] |= out[u] & m;
            out[u] |= m;
            fork[v] |= in[v] & m;
            in[v] |= m;
        }
    }
    return np;
}

void Combos::table(uint64_t sites) {
    for (Z = 1; Z < ZF * sites; Z <<= 1);
    seen.assign(Z, 0);
    ZT = std::max<uint64_t>(Z / ZT0, 1 << 16);
    key.assign(ZT, 0);
    n.assign(ZT * (S - 1), 0);
    best.assign(ZT * (S - 1), 0);
}

void Combos::grow() {
    ZT *= 2;
    full = 0;
    key.assign(ZT, 0);
    n.assign(ZT * (S - 1), 0);
    best.assign(ZT * (S - 1), 0);
}

bool Combos::kept(uint64_t k) const {
    uint32_t b = Graph::to(g.edges[in(k)]);
    return Graph::from(g.edges[out(k)]) == b ? branch[b] : keeps[seg[b]];
}

uint32_t Combos::edge(uint32_t u, uint32_t v) const {
    for (uint32_t i = a.first_out[u]; i < a.first_out[u + 1]; i++)
        if (a.target[i] == v) return i;
    return NONE;
}

void Combos::find(Chain& c, uint32_t m) const {
    std::vector<Hit>& out = c.hits;
    std::vector<uint32_t>& e = c.ed;
    const std::vector<uint32_t>& id = c.id;
    out.clear();
    e.resize(m);
    for (uint32_t i = 0; i + 1 < m; i++) __builtin_prefetch(&a.first_out[id[i]]);
    for (uint32_t i = 0; i + 1 < m; i++) __builtin_prefetch(&a.target[a.first_out[id[i]]]);
    for (uint32_t i = 0; i + 1 < m; i++) e[i] = edge(id[i], id[i + 1]);
    for (uint32_t i = 1; i + 1 < m; i++)
        if (branch[id[i]] && e[i - 1] != NONE && e[i] != NONE) out.push_back({(uint64_t)(e[i - 1] + 1) << 32 | e[i], i - 1, i + 1});
    for (uint32_t i = 1, j; i + 1 < m; i = j + 1) {
        uint32_t u = seg[id[i]];
        for (j = i; j + 1 < m && e[j] != NONE && seg[id[j + 1]] == u; j++);
        if (keeps[u] && g.starts[u] == id[i] && j - i + 1 == g.starts[u + 1] - g.starts[u] && j + 1 < m && e[i - 1] != NONE && e[j] != NONE)
            out.push_back({(uint64_t)(e[i - 1] + 1) << 32 | e[j], i - 1, j + 1});
    }
}

void Combos::remember(uint64_t k) {
    for (uint64_t i = home(k, Z);; i = (i + 1) & (Z - 1)) {
        uint64_t v = __atomic_load_n(&seen[i], __ATOMIC_RELAXED);
        if (v == k) return;
        if (!v && (__atomic_compare_exchange_n(&seen[i], &v, k, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED) || v == k)) return;
    }
}

bool Combos::known(uint64_t k) const {
    for (uint64_t i = home(k, Z);; i = (i + 1) & (Z - 1)) {
        uint64_t v = seen[i];
        if (v == k) return true;
        if (!v) return false;
    }
}

void Combos::add(uint64_t k, uint64_t s, uint64_t sc) {
    uint64_t i = home(k, ZT);
    for (uint32_t p = 0;; i = (i + 1) & (ZT - 1)) {
        uint64_t v = __atomic_load_n(&key[i], __ATOMIC_ACQUIRE);
        if (v == k) break;
        if (!v && __atomic_compare_exchange_n(&key[i], &v, k, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) break;
        if (v == k) break;
        if (++p == ZT / 2) { __atomic_store_n(&full, 1, __ATOMIC_RELAXED); return; }
    }
    uint64_t j = i * (S - 1) + s - 1, b = __atomic_load_n(&best[j], __ATOMIC_RELAXED);
    __atomic_fetch_add(&n[j], 1, __ATOMIC_RELAXED);
    while (sc > b && !__atomic_compare_exchange_n(&best[j], &b, sc, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
}

void Combos::tally(Chain& c, uint32_t m, uint64_t s, uint64_t rid) {
    if (s && __atomic_load_n(&full, __ATOMIC_RELAXED)) return;
    uint32_t i = 0;
    while (i < m && !on[c.id[i]]) i++;
    if (i == m) return;
    find(c, m);
    for (const Hit& h : c.hits) __builtin_prefetch(&seen[home(h.key, Z)]);
    for (const Hit& h : c.hits)
        if (!s) remember(h.key);
        else if (h.b - h.a == 2 && branch[c.id[h.a + 1]] && !keeps[seg[c.id[h.a + 1]]] && !(fork[c.id[h.a + 1]] >> s & 1)) continue;
        else if (!known(h.key)) add(h.key, s, centrality(c, h, rid));
}
