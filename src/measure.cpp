#include <algorithm>
#include <cstdlib>
#include "measure.h"

static const uint32_t CLT = 20, NONE = ~0u;

struct mctx { std::vector<site>* w; const keys* x; std::vector<uint32_t> ao, ai, tag; std::vector<std::vector<hit>> h; std::vector<std::vector<std::pair<uint32_t, uint32_t>>> ps; std::vector<std::vector<uint32_t>> st; std::vector<std::vector<uint64_t>> ks; };

static uint32_t only(const std::vector<std::pair<uint32_t, uint32_t>>& ps, uint32_t v) {
    auto a = std::lower_bound(ps.begin(), ps.end(), std::make_pair(v, 0u)), b = std::lower_bound(ps.begin(), ps.end(), std::make_pair(v, NONE));
    return b - a == 1 ? a->second : NONE;
}

static void take(site* w, uint32_t s, bool route, int32_t d, uint64_t sc, uint32_t lo, uint32_t hi, uint32_t a, uint32_t b) {
    while (__atomic_exchange_n(&w->lk, 1, __ATOMIC_ACQUIRE));
    if (!s) { if (w->nspan < DN) w->dn[w->nspan] = d; w->nspan++; }
    else if (s == w->s) {
        w->tspan++;
        if (route) {
            uint32_t k = 0, j;
            while (k < CL && w->cl[k].n && std::abs(w->cl[k].d - d) > (int32_t)CLT) k++;
            if (k == CL) for (k = 0, j = 1; j < CL; j++) if (std::abs(w->cl[j].d - d) < std::abs(w->cl[k].d - d)) k = j;
            clus* c = &w->cl[k];
            if (!c->n) c->d = d;
            c->n++;
            if (sc > c->best) { c->best = sc; c->lo = lo; c->hi = hi; c->a = a; c->b = b; }
        }
    }
    __atomic_store_n(&w->lk, 0, __ATOMIC_RELEASE);
}

static void look(rpass* p, rth* t, const uint8_t*, uint32_t m, uint32_t rc, uint64_t rid) {
    mctx* c = (mctx*)p->ctx;
    uint32_t th = t - p->t.data(), tag = ++c->tag[th];
    std::vector<std::pair<uint32_t, uint32_t>>& ps = c->ps[th];
    std::vector<uint64_t>& ks = c->ks[th];
    ps.clear();
    ks.clear();
    for (uint32_t i = 0; i < m; i++) ps.push_back({t->id[i], i});
    std::sort(ps.begin(), ps.end());
    chain_keys(c->x, t, m, c->h[th]);
    for (const hit& q : c->h[th]) ks.push_back(q.key);
    std::sort(ks.begin(), ks.end());
    const rdcore* r = &t->r;
    for (uint32_t i = 0; i < m; i++)
        for (uint32_t q = c->ao[t->id[i]]; q < c->ao[t->id[i] + 1]; q++) {
            uint32_t j = c->ai[q];
            if (c->st[th][j] == tag) continue;
            c->st[th][j] = tag;
            site* w = &(*c->w)[j];
            uint32_t li = 0, ri = 0, pi = NONE, pj = NONE;
            while (li < w->la.size() && (pi = only(ps, w->la[li])) == NONE) li++;
            while (ri < w->ra.size() && (pj = only(ps, w->ra[ri])) == NONE) ri++;
            if (pi == NONE || pj == NONE || pi >= pj) continue;
            uint32_t lo = pi, hi = pj;
            for (uint32_t k = li + 1; k < w->la.size(); k++) lo = std::min(lo, only(ps, w->la[k]));
            for (uint32_t k = ri + 1; k < w->ra.size(); k++) { uint32_t v = only(ps, w->ra[k]); if (v != NONE) hi = std::max(hi, v); }
            int32_t d = (int32_t)(r->a[t->ix[pj]] - r->a[t->ix[pi]]) - (int32_t)(w->rp[ri] - w->lp[li]);
            uint64_t sc = (uint64_t)std::min(r->a[t->ix[pi]], r->len - r->b[t->ix[pj]]) << 33 | (uint64_t)(rc != 0) << 32 | (0xffffffffu - (uint32_t)rid);
            take(w, p->f, std::binary_search(ks.begin(), ks.end(), c->x->k[w->key]), d, sc, r->a[t->ix[lo]], r->b[t->ix[hi]], r->a[t->ix[pi]], r->b[t->ix[pj]]);
        }
}

void measure(std::vector<site>& w, const keys* x, const char* reads, int t) {
    mctx c = {&w, x, std::vector<uint32_t>(x->g->K + 1, 0), {}, std::vector<uint32_t>(t, 0), std::vector<std::vector<hit>>(t), std::vector<std::vector<std::pair<uint32_t, uint32_t>>>(t), std::vector<std::vector<uint32_t>>(t, std::vector<uint32_t>(w.size(), 0)), std::vector<std::vector<uint64_t>>(t)};
    for (const site& s : w) { for (uint32_t v : s.la) c.ao[v + 1]++; for (uint32_t v : s.ra) c.ao[v + 1]++; }
    for (uint64_t v = 0; v < x->g->K; v++) c.ao[v + 1] += c.ao[v];
    c.ai.resize(c.ao[x->g->K]);
    std::vector<uint32_t> p(c.ao.begin(), c.ao.end() - 1);
    for (uint32_t i = 0; i < w.size(); i++) { for (uint32_t v : w[i].la) c.ai[p[v]++] = i; for (uint32_t v : w[i].ra) c.ai[p[v]++] = i; }
    rpass q;
    each_read(&q, x->g->g, reads, t, look, &c);
}
