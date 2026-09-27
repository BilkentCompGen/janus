#include <algorithm>
#include <string>
#include "io.h"
#include "keys.h"
#include "mates.h"
#include "select.h"

static const uint32_t MM = 3, PAD = 1000, SEG = 5000, NONE = ~0u;

struct unit { uint64_t read, key; uint32_t L, R, fa, fb, site; std::string seg; std::vector<uint32_t> mt; };

struct job { uint64_t read; uint32_t w, k, rc, fa, fb; };

struct rout { std::vector<std::pair<uint32_t, uint32_t>> iv; std::vector<std::string> s; std::vector<uint32_t> lo; };

struct sctx { keys* x; std::vector<mate>* w; std::vector<uint8_t> full; std::vector<unit> u; std::vector<job> jb; std::vector<uint64_t> rd, all; std::vector<uint32_t> ro; std::vector<rout> out; std::vector<std::vector<hit>> h; };

static bool by_read(const unit& p, const unit& q) { return p.read < q.read; }

static bool by_job(const job& p, const job& q) { return p.read < q.read; }

static bool odd(const keys* x, uint64_t i) {
    bool on = false, off = false;
    for (uint64_t s = 0; s < x->S; s++) { on |= x->n[i * x->S + s] >= MM; off |= x->n[i * x->S + s] < MM; }
    return on && off;
}

static void gather(sctx* c, const std::vector<uint8_t>& drop) {
    keys* x = c->x;
    const gph* g = x->g;
    for (uint64_t i = 0; i < x->Z; i++) {
        if (!x->k[i] || !odd(x, i)) continue;
        uint64_t ei = (x->k[i] >> 32) - 1, eo = x->k[i] & 0xffffffff;
        if (drop[g->e[ei] >> 32] || drop[g->e[ei] & 0xffffffff] || drop[g->e[eo] >> 32] || drop[g->e[eo] & 0xffffffff]) continue;
        for (uint64_t s = 0; s < x->S; s++)
            if (x->n[i * x->S + s]) c->u.push_back({s << 40 | (0xffffffffu - (uint32_t)x->best[i * x->S + s]), i, NONE, NONE, 0, 0, 0, "", {}});
    }
    std::sort(c->u.begin(), c->u.end(), by_read);
    for (size_t i = 0; i < c->u.size(); i++)
        if (!i || c->u[i].read != c->u[i - 1].read) { c->rd.push_back(c->u[i].read); c->ro.push_back(i); }
    c->ro.push_back(c->u.size());
}

static void anchor(rpass* p, rth* t, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) {
    sctx* c = (sctx*)p->ctx;
    uint64_t r = (uint64_t)p->f << 40 | rid, w = std::lower_bound(c->rd.begin(), c->rd.end(), r) - c->rd.begin();
    if (w == c->rd.size() || c->rd[w] != r) return;
    std::vector<hit>& h = c->h[t - p->t.data()];
    chain_keys(c->x, t, m, h);
    const rdcore* d = &t->r;
    for (uint32_t i = c->ro[w]; i < c->ro[w + 1]; i++)
        for (const hit& q : h) {
            if (q.key != c->x->k[c->u[i].key]) continue;
            unit* e = &c->u[i];
            uint32_t a = q.a, b = q.b;
            while (a > 0 && !c->full[t->id[a]]) a--;
            while (b + 1 < m && !c->full[t->id[b]]) b++;
            uint32_t st = d->a[t->ix[a]], en = d->b[t->ix[b]];
            e->L = c->full[t->id[a]] ? t->id[a] : NONE;
            e->R = c->full[t->id[b]] ? t->id[b] : NONE;
            e->site = en - st;
            e->fa = rc ? d->len - en : st;
            e->fb = rc ? d->len - st : en;
            if (e->L != NONE && e->R != NONE && en - st <= SEG) e->seg = bases(d, rec, st, en, rc);
        }
}

static void cut(rpass* p, rth* t, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) {
    sctx* c = (sctx*)p->ctx;
    uint64_t r = (uint64_t)p->f << 40 | rid, w = std::lower_bound(c->all.begin(), c->all.end(), r) - c->all.begin();
    if (w == c->all.size() || c->all[w] != r) return;
    const rdcore* d = &t->r;
    rout* o = &c->out[w];
    for (auto it = std::lower_bound(c->jb.begin(), c->jb.end(), job{r, 0, 0, 0, 0, 0}, by_job); it != c->jb.end() && it->read == r; it++) {
        if (it->rc != (rc != 0)) continue;
        mate* q = &(*c->w)[it->w];
        uint32_t k = it->k, i = 0, j;
        while (i < m && t->id[i] != (k == 2 ? q->R : q->L)) i++;
        for (j = k ? i : i + 1; k == 0 && j < m && t->id[j] != q->R; j++);
        uint32_t st = d->a[t->ix[i]], en = d->b[t->ix[j]], s = std::max(PAD, 3 * (k ? q->site : en - st));
        uint32_t lo = st > (k == 1 ? PAD : s) ? st - (k == 1 ? PAD : s) : 0, hi = std::min(d->len, en + (k == 2 ? PAD : s));
        o->iv.push_back(rc ? std::make_pair(d->len - hi, d->len - lo) : std::make_pair(lo, hi));
        it->fa = rc ? d->len - en : st;
        it->fb = rc ? d->len - st : en;
        if (k == 0 && en - st <= SEG) q->seg = bases(d, rec, st, en, rc);
    }
    if (!rc) return;
    uint64_t a = std::lower_bound(c->rd.begin(), c->rd.end(), r) - c->rd.begin();
    for (uint32_t i = a < c->rd.size() && c->rd[a] == r ? c->ro[a] : 0; a < c->rd.size() && c->rd[a] == r && i < c->ro[a + 1]; i++) {
        const unit& e = c->u[i];
        uint32_t pad = std::max(PAD, 3 * (e.fb - e.fa));
        o->iv.push_back({e.fa > pad ? e.fa - pad : 0, std::min(d->len, e.fb + pad)});
    }
    std::sort(o->iv.begin(), o->iv.end());
    std::vector<uint32_t> hi;
    for (auto& v : o->iv)
        if (!hi.empty() && v.first <= hi.back()) hi.back() = std::max(hi.back(), v.second);
        else { o->lo.push_back(v.first); hi.push_back(v.second); }
    for (uint32_t i = 0; i < o->lo.size(); i++) o->s.push_back(bases(d, rec, o->lo[i], hi[i], 0));
    o->iv = {};
}

static void pair_up(sctx* c) {
    std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> q;
    for (uint64_t i = 0; i < c->u.size(); i++)
        for (uint64_t s = 0; s < c->x->S; s++)
            if (s != c->u[i].read >> 40 && (c->u[i].L != NONE || c->u[i].R != NONE)) q.push_back({{(uint64_t)c->u[i].L << 32 | c->u[i].R, s}, i});
    std::sort(q.begin(), q.end());
    for (size_t i = 0; i < q.size(); i++) {
        const unit& e = c->u[q[i].second];
        if (!i || q[i].first != q[i - 1].first) c->w->push_back({e.L, e.R, (uint32_t)q[i].first.second, e.site, {0, 0, 0}, ""});
        c->u[q[i].second].mt.push_back(c->w->size() - 1);
    }
}

static void jobs(sctx* c) {
    const std::vector<mate>& w = *c->w;
    for (uint32_t i = 0; i < w.size(); i++)
        for (uint32_t k = 0; k < 3; k++)
            if (w[i].best[k] && (k == 0 || !w[i].best[0])) c->jb.push_back({(uint64_t)w[i].f << 40 | (0xffffffffu - (uint32_t)w[i].best[k]), i, k, (uint32_t)(w[i].best[k] >> 32 & 1), 0, 0});
    std::sort(c->jb.begin(), c->jb.end(), by_job);
    c->all = c->rd;
    for (const job& j : c->jb) c->all.push_back(j.read);
    std::sort(c->all.begin(), c->all.end());
    c->all.erase(std::unique(c->all.begin(), c->all.end()), c->all.end());
}

static std::string counts(const keys* x, uint64_t k) {
    std::string s;
    for (uint64_t i = 0; i < x->S; i++) s += (i ? "," : "") + std::to_string(x->n[k * x->S + i]);
    return s;
}

static std::string allele(const unit& e, const std::vector<mate>& w, uint64_t S) {
    std::string al;
    for (uint64_t s = 0, z = 0; s < S; s++) {
        if (s == e.read >> 40) { al += 'a'; continue; }
        const mate* q = z < e.mt.size() ? &w[e.mt[z++]] : nullptr;
        al += !q || q->seg.empty() || e.seg.empty() ? '-' : same(e.seg, q->seg) ? '1' : '0';
    }
    return al;
}

static std::string piece_of(const sctx* c, uint64_t read, uint32_t pos, uint32_t& off) {
    uint64_t i = std::lower_bound(c->all.begin(), c->all.end(), read) - c->all.begin();
    const std::vector<uint32_t>& lo = c->out[i].lo;
    uint32_t k = std::upper_bound(lo.begin(), lo.end(), pos) - lo.begin() - 1;
    off = lo[k];
    return "s" + std::to_string(read >> 40) + "_" + std::to_string(read & 0xffffffffffull) + "_" + std::to_string(k);
}

int cmd_select(const char* in, const char* dr, const char* rd, const char* out, int t, int mode) {
    gph g;
    adj a;
    keys x;
    sctx c;
    std::vector<mate> w;
    std::vector<uint64_t> m;
    read_gph(in, &g);
    make_adj(&g, &a);
    masks(&g, m);
    c.full.resize(g.K);
    for (uint64_t v = 0; v < g.K; v++) c.full[v] = m[v] == (1ull << g.S) - 1;
    m = {};
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
    rpass p;
    each_read(&p, g.g, rd, t, anchor, &c);
    pair_up(&c);
    find_mates(w, g.g, rd, t);
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
    std::vector<std::vector<uint32_t>> mj(w.size());
    for (uint32_t i = 0; i < c.jb.size(); i++) mj[c.jb[i].w].push_back(i);
    for (size_t j = 0; j < c.u.size(); j++) {
        const unit& e = c.u[j];
        uint32_t off;
        std::string nm = piece_of(&c, e.read, e.fa, off);
        std::fprintf(v, "%s\t%u\t%u\talt\t%zu\t%s\t%s\n", nm.c_str(), e.fa - off, e.fb - off, j, counts(&x, e.key).c_str(), allele(e, w, x.S).c_str());
        for (uint32_t i : e.mt)
            for (uint32_t q : mj[i]) {
                const job& b = c.jb[q];
                nm = piece_of(&c, b.read, b.fa, off);
                std::fprintf(v, "%s\t%u\t%u\tcp\t%zu\t%u\t%u\n", nm.c_str(), b.fa - off, b.fb - off, j, w[i].f, b.k);
            }
    }
    std::fclose(f);
    std::fclose(v);
    std::fprintf(stderr, "[select] %zu evidence units on %zu reads; %zu other-sample requests, %zu reads chosen for them; %zu reads, %lu pieces, %lu bp\n", c.u.size(), c.rd.size(), w.size(), c.jb.size(), c.all.size(), np, nb);
    return 0;
}
