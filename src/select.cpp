#include <algorithm>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include "io.h"
#include "keys.h"
#include "measure.h"
#include "select.h"

static const uint32_t MM = 3, PAD = 1000, SIDE = 48, SAME = 20, NONE = ~0u;

struct job { uint64_t read; uint32_t w, k, rc, lo, hi, a, b, pc; };

struct rout { std::vector<std::pair<uint32_t, uint32_t>> iv; std::vector<std::string> s; };

struct sctx { keys* x; std::vector<site>* w; std::unordered_map<uint64_t, uint32_t> dix; std::vector<uint32_t> dep; std::vector<uint64_t> mask; const uint8_t* sh; std::vector<job> jb; std::vector<uint64_t> rd, all; std::vector<uint32_t> ro; std::vector<rout> out; std::vector<std::vector<hit>> h; std::vector<std::vector<uint32_t>> ids; };

static bool by_read(const site& p, const site& q) { return p.read < q.read; }

static bool by_job(const job& p, const job& q) { return p.read < q.read; }

static bool odd(const keys* x, uint64_t i) {
    bool on = false;
    for (uint64_t s = 1; s < x->S; s++) on |= x->n[i * x->S + s] >= MM;
    return on && !x->n[i * x->S];
}

static void gather(sctx* c, const std::vector<uint8_t>& drop) {
    keys* x = c->x;
    const gph* g = x->g;
    std::vector<std::pair<uint64_t, uint64_t>> at;
    for (uint64_t i = 0; i < x->Z; i++) {
        if (!x->k[i]) continue;
        uint64_t ei = (x->k[i] >> 32) - 1, eo = x->k[i] & 0xffffffff, b = g->e[ei] & 0xffffffff;
        at.push_back({(g->e[eo] >> 32) == b ? b : g->K + x->unit[b], i});
    }
    std::sort(at.begin(), at.end());
    for (size_t i = 0, j; i < at.size(); i = j) {
        std::vector<uint32_t> d(x->S, 0);
        bool any = false;
        for (j = i; j < at.size() && at[j].first == at[i].first; j++) {
            for (uint64_t s = 0; s < x->S; s++) d[s] += x->n[at[j].second * x->S + s];
            any |= odd(x, at[j].second);
        }
        if (!any) continue;
        for (size_t q = i; q < j; q++) c->dix[at[q].second] = c->dep.size();
        c->dep.insert(c->dep.end(), d.begin(), d.end());
    }
    at = {};
    for (uint64_t i = 0; i < x->Z; i++) {
        if (!x->k[i] || !odd(x, i)) continue;
        uint64_t ei = (x->k[i] >> 32) - 1, eo = x->k[i] & 0xffffffff;
        if (drop[g->e[ei] >> 32] || drop[g->e[ei] & 0xffffffff] || drop[g->e[eo] >> 32] || drop[g->e[eo] & 0xffffffff]) continue;
        for (uint64_t s = 1; s < x->S; s++)
            if (x->n[i * x->S + s] >= MM) { site w = {}; w.read = s << 40 | (0xffffffffu - (uint32_t)x->best[i * x->S + s]); w.key = i; w.s = s; c->w->push_back(w); }
    }
    std::sort(c->w->begin(), c->w->end(), by_read);
    for (size_t i = 0; i < c->w->size(); i++)
        if (!i || (*c->w)[i].read != (*c->w)[i - 1].read) { c->rd.push_back((*c->w)[i].read); c->ro.push_back(i); }
    c->ro.push_back(c->w->size());
}

static bool anc(const sctx* c, const std::vector<uint32_t>& ids, uint32_t v, uint64_t s) {
    auto r = std::equal_range(ids.begin(), ids.end(), v);
    return (c->mask[v] & (1ull | 1ull << s)) == (1ull | 1ull << s) && r.second - r.first == 1;
}

static void anchor(rpass* p, rth* t, const uint8_t*, uint32_t m, uint32_t rc, uint64_t rid) {
    sctx* c = (sctx*)p->ctx;
    uint64_t r = (uint64_t)p->f << 40 | rid, w = std::lower_bound(c->rd.begin(), c->rd.end(), r) - c->rd.begin();
    if (w == c->rd.size() || c->rd[w] != r) return;
    std::vector<hit>& h = c->h[t - p->t.data()];
    chain_keys(c->x, t, m, h);
    std::vector<uint32_t>& ids = c->ids[t - p->t.data()];
    ids.assign(t->id.begin(), t->id.begin() + m);
    std::sort(ids.begin(), ids.end());
    const rdcore* d = &t->r;
    for (uint32_t i = c->ro[w]; i < c->ro[w + 1]; i++)
        for (const hit& q : h) {
            site* e = &(*c->w)[i];
            if (q.key != c->x->k[e->key]) continue;
            uint32_t a = q.a, b = q.b;
            e->la.clear();
            e->ra.clear();
            e->lp.clear();
            e->rp.clear();
            for (int64_t j = q.a; j >= 0 && e->la.size() < SIDE; j--)
                if (anc(c, ids, t->id[j], p->f)) { if (e->la.empty()) a = j; e->la.push_back(t->id[j]); e->lp.push_back(d->a[t->ix[j]]); }
            for (uint32_t j = q.b; j < m && e->ra.size() < SIDE; j++)
                if (anc(c, ids, t->id[j], p->f)) { if (e->ra.empty()) b = j; e->ra.push_back(t->id[j]); e->rp.push_back(d->a[t->ix[j]]); }
            e->rc = rc != 0;
            e->fa = d->a[t->ix[q.a]];
            e->fb = d->b[t->ix[q.b]];
            e->a0 = d->a[t->ix[a]];
            e->b0 = d->b[t->ix[b]];
            e->span = e->b0 - e->a0;
            e->tonly = e->shadow = e->nbr = 0;
            for (uint32_t j = a; j <= b; j++) {
                uint32_t v = t->id[j];
                e->nbr += c->x->br[v];
                if (!(c->mask[v] & 1)) { e->tonly++; e->shadow = std::max<uint32_t>(e->shadow, c->sh[v]); }
            }
        }
}

static void stats(site* w) {
    uint32_t n = std::min(w->nspan, DN);
    if (!n) return;
    std::sort(w->dn, w->dn + n);
    w->lest = w->dn[n / 2];
    w->dmin = w->dn[0];
    w->dmax = w->dn[n - 1];
    for (uint32_t i = 0; i < n; i++) w->nsame += std::abs(w->dn[i]) <= std::max<int32_t>(SAME, std::abs(w->lest) / 2);
}

static void jobs(sctx* c) {
    std::vector<site>& w = *c->w;
    for (uint32_t i = 0; i < w.size(); i++) {
        uint32_t n = 0;
        for (uint32_t k = 0; k < CL; k++)
            if (w[i].cl[k].n) { const clus& q = w[i].cl[k]; c->jb.push_back({(uint64_t)w[i].s << 40 | (0xffffffffu - (uint32_t)q.best), i, k, (uint32_t)(q.best >> 32 & 1), q.lo, q.hi, q.a, q.b, 0}); n++; }
        if (!n) c->jb.push_back({w[i].read, i, NONE, w[i].rc, w[i].fa, w[i].fb, w[i].fa, w[i].fb, 0});
    }
    std::sort(c->jb.begin(), c->jb.end(), by_job);
    for (const job& j : c->jb) c->all.push_back(j.read);
    std::sort(c->all.begin(), c->all.end());
    c->all.erase(std::unique(c->all.begin(), c->all.end()), c->all.end());
}

static void cut(rpass* p, rth* t, const uint8_t* rec, uint32_t, uint32_t rc, uint64_t rid) {
    sctx* c = (sctx*)p->ctx;
    uint64_t r = (uint64_t)p->f << 40 | rid, w = std::lower_bound(c->all.begin(), c->all.end(), r) - c->all.begin();
    if (w == c->all.size() || c->all[w] != r) return;
    const rdcore* d = &t->r;
    rout* o = &c->out[w];
    for (auto it = std::lower_bound(c->jb.begin(), c->jb.end(), job{r, 0, 0, 0, 0, 0, 0, 0, 0}, by_job); it != c->jb.end() && it->read == r; it++) {
        if (it->rc != (rc != 0)) continue;
        uint32_t pd = it->k == NONE ? std::max(PAD, 3 * (it->hi - it->lo)) : PAD, lo = it->lo > pd ? it->lo - pd : 0, hi = std::min(d->len, it->hi + pd), a = it->a, b = it->b;
        if (rc) { uint32_t z = d->len - hi; hi = d->len - lo; lo = z; z = d->len - b; b = d->len - a; a = z; }
        it->pc = o->iv.size();
        it->a = a - lo;
        it->b = b - lo;
        o->iv.push_back({lo, hi});
        o->s.push_back(bases(d, rec, lo, hi, 0));
    }
}

static std::string list(const uint32_t* v, uint64_t n) {
    std::string s;
    for (uint64_t i = 0; i < n; i++) s += (i ? "," : "") + std::to_string(v[i]);
    return s;
}

int cmd_select(const char* in, const char* dr, const char* rd, const char* out, int t, int mode) {
    gph g;
    adj a;
    keys x;
    sctx c;
    std::vector<site> w;
    read_gph(in, &g);
    make_adj(&g, &a);
    masks(&g, c.mask);
    c.sh = g.sh.data();
    g.c = {};
    g.ec = {};
    make_keys(&x, &g, &a, mode);
    count_keys(&x, rd, t);
    std::vector<uint8_t> drop(g.K);
    for (uint64_t h : read_set(dr)) drop[std::lower_bound(g.g.begin(), g.g.end(), h) - g.g.begin()] = 1;
    c.x = &x;
    c.w = &w;
    gather(&c, drop);
    c.h.resize(t);
    c.ids.resize(t);
    rpass p;
    each_read(&p, g.g, rd, t, anchor, &c);
    measure(w, &x, rd, t);
    for (site& s : w) stats(&s);
    jobs(&c);
    c.out.resize(c.all.size());
    each_read(&p, g.g, rd, t, cut, &c);
    std::string o(out);
    FILE *f = std::fopen((o + ".fa").c_str(), "w"), *v = std::fopen((o + ".tsv").c_str(), "w");
    uint64_t nb = 0, np = 0;
    for (size_t i = 0; i < c.all.size(); i++)
        for (uint32_t k = 0; k < c.out[i].s.size(); k++) {
            std::fprintf(f, ">s%lu_%lu_%u\n%s\n", (uint64_t)(c.all[i] >> 40), (uint64_t)(c.all[i] & 0xffffffffff), k, c.out[i].s[k].c_str());
            nb += c.out[i].s[k].size();
            np++;
        }
    for (const job& j : c.jb) {
        const site& e = w[j.w];
        uint32_t ncl = 0;
        for (uint32_t k = 0; k < CL; k++) ncl += e.cl[k].n > 0;
        std::fprintf(v, "s%lu_%lu_%u\t%u\t%u\t%u\t%s\t%s\t%lu\t%u\t%u\t%u\t%u\t%u\t%u\t%d\t%d\t%d\t%u\t%u\t%u\t%u\n", (uint64_t)(j.read >> 40), (uint64_t)(j.read & 0xffffffffff), j.pc, j.a, j.b, j.w,
                     list(&x.n[e.key * x.S], x.S).c_str(), list(&c.dep[c.dix[e.key]], x.S).c_str(), (unsigned long)e.key, e.s, ncl, j.k == NONE ? 0 : e.cl[j.k].n, e.tspan, e.nspan, e.nsame, e.lest, e.dmin, e.dmax, e.shadow, e.tonly, e.nbr, e.span);
    }
    std::fclose(f);
    std::fclose(v);
    std::fprintf(stderr, "[select] %zu sites on %zu reads; %zu pieces on %zu reads, %lu bp\n", w.size(), c.rd.size(), np, c.all.size(), nb);
    return 0;
}
