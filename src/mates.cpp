#include <algorithm>
#include "mates.h"
#include "reads.h"

static const uint32_t NONE = ~0u;

struct mctx { std::vector<mate>* w; std::vector<uint32_t> lo, li, ro, ri; std::vector<std::vector<std::pair<uint32_t, uint32_t>>> pos; };

static void offer(uint64_t* b, uint64_t v) {
    uint64_t o = __atomic_load_n(b, __ATOMIC_RELAXED);
    while (v > o && !__atomic_compare_exchange_n(b, &o, v, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
}

static uint64_t reach(uint32_t in, uint32_t out, uint32_t site) {
    return (uint64_t)std::min(in, site + 1000) << 10 | std::min(out, 1023u);
}

static bool seen(const std::vector<std::pair<uint32_t, uint32_t>>& ps, uint32_t v, uint32_t lo, uint32_t hi) {
    auto it = std::lower_bound(ps.begin(), ps.end(), std::make_pair(v, lo));
    return v == NONE || (it != ps.end() && it->first == v && it->second < hi);
}

static void score(rpass* p, rth* t, const uint8_t*, uint32_t m, uint32_t rc, uint64_t rid) {
    mctx* c = (mctx*)p->ctx;
    std::vector<std::pair<uint32_t, uint32_t>>& ps = c->pos[t - p->t.data()];
    ps.clear();
    for (uint32_t i = 0; i < m; i++) ps.push_back({t->id[i], i});
    std::sort(ps.begin(), ps.end());
    const rdcore* r = &t->r;
    uint64_t tag = (uint64_t)(rc != 0) << 32 | (0xffffffffu - (uint32_t)rid);
    for (uint32_t i = 0; i < m; i++) {
        uint32_t v = t->id[i], st = r->a[t->ix[i]], en = r->b[t->ix[i]];
        for (uint32_t q = c->lo[v]; q < c->lo[v + 1]; q++) {
            mate* w = &(*c->w)[c->li[q]];
            if (w->f != (uint32_t)p->f) continue;
            if (seen(ps, w->L2, 0, i)) offer(&w->best[1], reach(r->len - en, st, w->site) << 33 | tag);
            auto it = std::lower_bound(ps.begin(), ps.end(), std::make_pair(w->R, i + 1));
            if (w->R != NONE && it != ps.end() && it->first == w->R) offer(&w->best[0], (uint64_t)std::min(st, r->len - r->b[t->ix[it->second]]) << 33 | tag);
        }
        for (uint32_t q = c->ro[v]; q < c->ro[v + 1]; q++) {
            mate* w = &(*c->w)[c->ri[q]];
            if (w->f == (uint32_t)p->f && seen(ps, w->R2, i + 1, m)) offer(&w->best[2], reach(st, r->len - en, w->site) << 33 | tag);
        }
    }
}

static void by_node(const std::vector<mate>& w, uint64_t K, bool right, std::vector<uint32_t>& o, std::vector<uint32_t>& ix) {
    o.assign(K + 1, 0);
    for (const mate& x : w) if ((right ? x.R : x.L) != NONE) o[(right ? x.R : x.L) + 1]++;
    for (uint64_t i = 0; i < K; i++) o[i + 1] += o[i];
    ix.resize(o[K]);
    std::vector<uint32_t> p(o.begin(), o.end() - 1);
    for (uint32_t i = 0; i < w.size(); i++) if ((right ? w[i].R : w[i].L) != NONE) ix[p[right ? w[i].R : w[i].L]++] = i;
}

void find_mates(std::vector<mate>& w, const std::vector<uint64_t>& g, const char* reads, int t) {
    mctx c;
    c.w = &w;
    c.pos.resize(t);
    by_node(w, g.size(), false, c.lo, c.li);
    by_node(w, g.size(), true, c.ro, c.ri);
    rpass p;
    each_read(&p, g, reads, t, score, &c);
}
