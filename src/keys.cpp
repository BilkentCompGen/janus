#include <algorithm>
#include "io.h"
#include "keys.h"

static const uint64_t ZF = 16; // table slots per possible combination, rounded up to a power of two

void make_keys(keys* x, const gph* g, const adj* a) {
    x->S = g->S;
    x->g = g;
    x->a = a;
    x->seg.resize(g->K);
    for (uint32_t u = 0; u + 1 < g->st.size(); u++) for (uint32_t v = g->st[u]; v < g->st[u + 1]; v++) x->seg[v] = u;
}

uint64_t keep_keys(keys* x, const std::vector<uint8_t>& sd) {
    const adj* a = x->a;
    const std::vector<uint32_t>& st = x->g->st;
    uint64_t np = 0, K = x->g->K;
    x->br.assign(K, 0);
    for (uint64_t v = 0; v < K; v++) np += x->br[v] = sd[v] && (a->io[v + 1] - a->io[v] > 1 || a->oo[v + 1] - a->oo[v] > 1);
    x->ks.assign(st.size() - 1, 0);
    for (uint64_t v = 0; v < K; v++) x->ks[x->seg[v]] |= sd[v] && st[x->seg[v] + 1] - st[x->seg[v]] > 1;
    for (uint8_t k : x->ks) np += k;
    x->on.assign(K, 0);
    for (uint64_t v = 0; v < K; v++) x->on[v] = x->br[v] || x->ks[x->seg[v]];
    return np;
}

void make_table(keys* x, uint64_t np) {
    for (x->Z = 1; x->Z < ZF * np; x->Z <<= 1);
    x->key.assign(x->Z, 0);
    x->n.assign(x->Z * x->S, 0);
    x->best.assign(x->Z * x->S, 0);
}

bool kept(const keys* x, uint64_t key) {
    uint64_t ei = (key >> 32) - 1, eo = key & 0xffffffff, b = x->g->e[ei] & 0xffffffff;
    return (x->g->e[eo] >> 32) == b ? x->br[b] : x->ks[x->seg[b]];
}

static uint32_t edge(const adj* a, uint32_t u, uint32_t v) {
    for (uint32_t i = a->oo[u]; i < a->oo[u + 1]; i++) if (a->on[i] == v) return i;
    return NONE;
}

void chain_keys(const keys* x, rth* c, uint32_t m, std::vector<hit>& out) {
    out.clear();
    std::vector<uint32_t>& e = c->ed;
    e.resize(m);
    for (uint32_t i = 0; i + 1 < m; i++) __builtin_prefetch(&x->a->oo[c->id[i]]);
    for (uint32_t i = 0; i + 1 < m; i++) __builtin_prefetch(&x->a->on[x->a->oo[c->id[i]]]);
    for (uint32_t i = 0; i + 1 < m; i++) e[i] = edge(x->a, c->id[i], c->id[i + 1]);
    for (uint32_t i = 1; i + 1 < m; i++)
        if (x->br[c->id[i]] && e[i - 1] != NONE && e[i] != NONE) out.push_back({(uint64_t)(e[i - 1] + 1) << 32 | e[i], i - 1, i + 1});
    for (uint32_t i = 1, j; i + 1 < m; i = j + 1) {
        uint32_t u = x->seg[c->id[i]];
        for (j = i; j + 1 < m && e[j] != NONE && x->seg[c->id[j + 1]] == u; j++);
        if (x->ks[u] && x->g->st[u] == c->id[i] && j - i + 1 == x->g->st[u + 1] - x->g->st[u] && j + 1 < m && e[i - 1] != NONE && e[j] != NONE)
            out.push_back({(uint64_t)(e[i - 1] + 1) << 32 | e[j], i - 1, j + 1});
    }
}

static void add(keys* x, uint64_t key, uint64_t s, uint64_t sc) {
    uint64_t i = (key * 0x9e3779b97f4a7c15ull) >> 20 & (x->Z - 1);
    for (;; i = (i + 1) & (x->Z - 1)) {
        uint64_t v = __atomic_load_n(&x->key[i], __ATOMIC_ACQUIRE);
        if (v == key) break;
        if (!v && __atomic_compare_exchange_n(&x->key[i], &v, key, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) break;
        if (v == key) break;
    }
    uint64_t j = i * x->S + s, b = __atomic_load_n(&x->best[j], __ATOMIC_RELAXED);
    __atomic_fetch_add(&x->n[j], 1, __ATOMIC_RELAXED);
    while (sc > b && !__atomic_compare_exchange_n(&x->best[j], &b, sc, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
}

void tally(keys* x, rth* c, uint32_t m, uint64_t s, uint64_t rid, std::vector<hit>& h) {
    uint32_t i = 0;
    while (i < m && !x->on[c->id[i]]) i++;
    if (i == m) return;
    chain_keys(x, c, m, h);
    for (const hit& o : h) {
        uint32_t st = c->r.a[c->ix[o.a]], en = c->r.b[c->ix[o.b]];
        add(x, o.key, s, (uint64_t)std::min(st, c->r.len - en) << 32 | (0xffffffffu - (uint32_t)rid));
    }
}
