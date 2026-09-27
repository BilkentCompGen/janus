#include <algorithm>
#include "keys.h"

static const uint32_t NONE = ~0u;

void make_keys(keys* x, const gph* g, const adj* a, int mode) {
    uint64_t nb = 0, nu = 0;
    x->S = g->S;
    x->g = g;
    x->a = a;
    x->mode = mode;
    x->unit.assign(g->K, NONE);
    x->off.assign(g->K, 0);
    x->br.assign(g->K, 0);
    for (uint64_t v = 0; v < g->K; v++) nb += x->br[v] = a->io[v + 1] - a->io[v] > 1 || a->oo[v + 1] - a->oo[v] > 1;
    for (int pass = 0; pass < 2; pass++)
        for (uint64_t v = 0; v < g->K; v++) {
            bool head = a->io[v + 1] - a->io[v] != 1 || a->oo[a->in[a->io[v]] + 1] - a->oo[a->in[a->io[v]]] != 1;
            if (x->unit[v] != NONE || (!pass && !head)) continue;
            uint32_t u = x->len.size(), w = v, i = 0;
            for (;;) {
                x->unit[w] = u;
                x->off[w] = i++;
                if (a->oo[w + 1] - a->oo[w] != 1) break;
                uint32_t nx = a->on[a->oo[w]];
                if (a->io[nx + 1] - a->io[nx] != 1 || x->unit[nx] != NONE) break;
                w = nx;
            }
            x->len.push_back(i);
            nu++;
        }
    for (x->Z = 1; x->Z < 4 * (nu + nb); x->Z <<= 1);
    x->k.assign(x->Z, 0);
    x->n.assign(x->Z * x->S, 0);
    x->best.assign(x->Z * x->S, 0);
}

static uint32_t edge(const adj* a, uint32_t u, uint32_t v) {
    for (uint32_t i = a->oo[u]; i < a->oo[u + 1]; i++) if (a->on[i] == v) return i;
    return NONE;
}

void chain_keys(const keys* x, rth* c, uint32_t m, std::vector<hit>& out) {
    out.clear();
    std::vector<uint32_t>& e = c->ed;
    e.resize(m);
    for (uint32_t i = 0; i + 1 < m; i++) e[i] = edge(x->a, c->id[i], c->id[i + 1]);
    for (uint32_t i = 1; i + 1 < m; i++)
        if (x->br[c->id[i]] && e[i - 1] != NONE && e[i] != NONE) out.push_back({(uint64_t)(e[i - 1] + 1) << 32 | e[i], i - 1, i + 1});
    for (uint32_t i = 1, j; i + 1 < m; i = j + 1) {
        uint32_t u = x->unit[c->id[i]];
        for (j = i; j + 1 < m && e[j] != NONE && x->unit[c->id[j + 1]] == u; j++);
        if (x->len[u] > 1 && x->off[c->id[i]] == 0 && j - i + 1 == x->len[u] && j + 1 < m && e[i - 1] != NONE && e[j] != NONE)
            out.push_back({(uint64_t)(e[i - 1] + 1) << 32 | e[j], i - 1, j + 1});
    }
}

uint64_t slot_of(const keys* x, uint64_t key) {
    uint64_t i = (key * 0x9e3779b97f4a7c15ull) >> 20 & (x->Z - 1);
    while (x->k[i] && x->k[i] != key) i = (i + 1) & (x->Z - 1);
    return i;
}

static void add(keys* x, uint64_t key, uint64_t s, uint64_t sc) {
    uint64_t i = (key * 0x9e3779b97f4a7c15ull) >> 20 & (x->Z - 1);
    for (;; i = (i + 1) & (x->Z - 1)) {
        uint64_t v = __atomic_load_n(&x->k[i], __ATOMIC_ACQUIRE);
        if (v == key) break;
        if (!v && __atomic_compare_exchange_n(&x->k[i], &v, key, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) break;
        if (v == key) break;
    }
    uint64_t j = i * x->S + s, b = __atomic_load_n(&x->best[j], __ATOMIC_RELAXED);
    __atomic_fetch_add(&x->n[j], 1, __ATOMIC_RELAXED);
    while (sc > b && !__atomic_compare_exchange_n(&x->best[j], &b, sc, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
}

struct kctx { keys* x; std::vector<std::vector<hit>> h; };

static void tally(rpass* p, rth* c, const uint8_t*, uint32_t m, uint32_t, uint64_t rid) {
    kctx* q = (kctx*)p->ctx;
    std::vector<hit>& h = q->h[c - p->t.data()];
    chain_keys(q->x, c, m, h);
    for (const hit& o : h) {
        uint32_t st = c->r.a[c->ix[o.a]], en = c->r.b[c->ix[o.b]], f = std::min(st, c->r.len - en), w = q->x->mode == 2 ? std::min(f, 4095u) << 19 | std::min(c->r.len >> 2, (1u << 19) - 1) : q->x->mode ? f : c->r.len;
        add(q->x, o.key, p->f, (uint64_t)w << 32 | (0xffffffffu - (uint32_t)rid));
    }
}

void count_keys(keys* x, const char* reads, int t) {
    kctx q = {x, std::vector<std::vector<hit>>(t)};
    rpass p;
    each_read(&p, x->g->g, reads, t, tally, &q);
}
