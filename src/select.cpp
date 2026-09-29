#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <string>
#include "io.h"
#include "keys.h"
#include "cmd.h"

static const uint32_t MM = 3, PAD = 1000, NA = 3, MULT = 2;
static const double DL = 10000;

struct anc { double sc; uint32_t v, a, b; };

struct site { uint64_t read, key, slot; uint32_t s, lo, hi, a, b, pc; std::vector<anc> l, r; };

struct piece { uint64_t read; uint32_t lo, hi, k; std::string s; };

struct sctx { keys* x; std::vector<site> w; std::vector<piece> pc; std::vector<uint64_t> rd; std::vector<uint32_t> ro, med, gs; const uint16_t* cn; uint64_t K; std::vector<std::vector<hit>> h; std::vector<std::vector<uint32_t>> ids; std::vector<int> fs, fx; int t; };

static bool by_read(const site& p, const site& q) { return p.read != q.read ? p.read < q.read : p.key < q.key; }

static bool by_score(const anc& p, const anc& q) { return p.sc > q.sc; }

static bool by_near(const anc& p, const anc& q) { return p.sc < q.sc; }

static bool odd(const keys* x, uint64_t i) {
    bool on = false;
    for (uint64_t s = 1; s < x->S; s++) on |= x->n[i * x->S + s] >= MM;
    return on && !x->n[i * x->S];
}

static void gather(sctx* c, const std::vector<uint8_t>& drop) {
    keys* x = c->x;
    const gph* g = x->g;
    for (uint64_t i = 0; i < x->Z; i++) {
        if (!x->k[i] || !odd(x, i)) continue;
        uint64_t ei = (x->k[i] >> 32) - 1, eo = x->k[i] & 0xffffffff;
        if (drop[g->e[ei] >> 32] || drop[g->e[ei] & 0xffffffff] || drop[g->e[eo] >> 32] || drop[g->e[eo] & 0xffffffff]) continue;
        for (uint64_t s = 1; s < x->S; s++)
            if (x->n[i * x->S + s] >= MM) { site w = {}; w.read = s << 40 | (0xffffffffu - (uint32_t)x->best[i * x->S + s]); w.key = x->k[i]; w.slot = i; w.s = s; c->w.push_back(w); }
    }
    std::sort(c->w.begin(), c->w.end(), by_read);
    for (size_t i = 0; i < c->w.size(); i++)
        if (!i || c->w[i].read != c->w[i - 1].read) { c->rd.push_back(c->w[i].read); c->ro.push_back(i); }
    c->ro.push_back(c->w.size());
}

static bool full(const gph* g, uint64_t v) {
    for (uint64_t s = 0; s < g->S; s++) if (!g->c[s * g->K + v]) return false;
    return true;
}

static double fit(const sctx* c, const std::vector<uint32_t>& ids, uint32_t v, uint64_t s) {
    auto r = std::equal_range(ids.begin(), ids.end(), v);
    if (!c->cn[v] || !c->cn[s * c->K + v] || r.second - r.first != 1 || c->cn[v] > MULT * c->med[0] || c->cn[s * c->K + v] > MULT * c->med[s]) return 0;
    double n = (double)c->cn[v] / c->med[0], t = (double)c->cn[s * c->K + v] / c->med[s];
    return std::min(n, 1 / n) * std::min(t, 1 / t);
}

static void side(const sctx* c, const rth* t, const std::vector<uint32_t>& ids, uint64_t s, int64_t from, int64_t step, int64_t end, std::vector<anc>& o) {
    const rdcore* d = &t->r;
    uint32_t at = d->a[t->ix[from]];
    for (int64_t j = from; j != end; j += step) {
        uint32_t p = d->a[t->ix[j]];
        double w = std::exp(-(double)(p > at ? p - at : at - p) / DL), f;
        if (o.size() == NA && w < o.back().sc) break;
        if ((f = fit(c, ids, t->id[j], s)) == 0) continue;
        o.push_back({f * w, t->id[j], p, d->b[t->ix[j]]});
        std::sort(o.begin(), o.end(), by_score);
        if (o.size() > NA) o.pop_back();
    }
    for (anc& x : o) x.sc = x.a > at ? x.a - at : at - x.a;
    std::sort(o.begin(), o.end(), by_near);
}

static void flip(uint32_t len, uint32_t rc, uint32_t* a, uint32_t* b) {
    if (!rc) return;
    uint32_t z = len - *b;
    *b = len - *a;
    *a = z;
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
            site* e = &c->w[i];
            if (q.key != e->key) continue;
            e->l.clear();
            e->r.clear();
            side(c, t, ids, p->f, q.a, -1, -1, e->l);
            side(c, t, ids, p->f, q.b, 1, m, e->r);
            e->a = e->l.empty() ? d->a[t->ix[q.a]] : e->l[0].a;
            e->b = e->r.empty() ? d->b[t->ix[q.b]] : e->r[0].b;
            uint32_t pd = std::max(PAD, 3 * (e->b - e->a)), lo = e->l.empty() ? e->a : e->l.back().a, hi = e->r.empty() ? e->b : e->r.back().b;
            e->lo = lo > pd ? lo - pd : 0;
            e->hi = std::min(d->len, hi + pd);
            flip(d->len, rc, &e->lo, &e->hi);
            flip(d->len, rc, &e->a, &e->b);
            for (anc& x : e->l) flip(d->len, rc, &x.a, &x.b);
            for (anc& x : e->r) flip(d->len, rc, &x.a, &x.b);
            if (rc) e->l.swap(e->r);
        }
}

static bool before(const std::pair<uint64_t, uint32_t>& v, const piece& q) { return v.first < q.read || (v.first == q.read && v.second < q.lo); }

static void pieces(sctx* c) {
    for (uint32_t k = 0; k + 1 < c->ro.size(); k++) {
        std::vector<std::pair<uint32_t, uint32_t>> iv;
        for (uint32_t i = c->ro[k]; i < c->ro[k + 1]; i++) if (c->w[i].hi > c->w[i].lo) iv.push_back({c->w[i].lo, c->w[i].hi});
        std::sort(iv.begin(), iv.end());
        for (uint32_t j = 0; j < iv.size(); j++) {
            if (j && iv[j].first <= c->pc.back().hi) { c->pc.back().hi = std::max(c->pc.back().hi, iv[j].second); continue; }
            if (!j) c->gs.push_back(c->pc.size());
            c->pc.push_back({c->rd[k], iv[j].first, iv[j].second, j ? c->pc.back().k + 1 : 0, ""});
        }
    }
    c->gs.push_back(c->pc.size());
    for (site& e : c->w) e.pc = std::upper_bound(c->pc.begin(), c->pc.end(), std::make_pair(e.read, e.lo), before) - c->pc.begin() - 1;
}

static void cut(void* a, int i) {
    sctx* c = (sctx*)a;
    rdcore r;
    std::vector<uint8_t> rec;
    for (size_t k = i; k + 1 < c->gs.size(); k += c->t) {
        uint64_t rd = c->pc[c->gs[k]].read, off;
        uint32_t tag, s = rd >> 40;
        pread(c->fx[s], &off, 8, sizeof(bin_header) + 8 * (rd & 0xffffffffff));
        pread(c->fs[s], &tag, 4, off);
        rec.assign(4 + (tag >> 31 ? tag & 0x7fffffff : ((tag & 0x7fffffff) + 3) / 4) + 8, 0);
        pread(c->fs[s], rec.data(), rec.size() - 8, off);
        collapse(&r, rec.data());
        for (uint32_t j = c->gs[k]; j < c->gs[k + 1]; j++) c->pc[j].s = bases(&r, rec.data(), c->pc[j].lo, c->pc[j].hi, 0);
    }
}

static std::string list(const std::vector<anc>& v, uint32_t lo) {
    std::string s;
    for (const anc& x : v) s += (s.empty() ? "" : ",") + std::to_string(x.v) + ":" + std::to_string(x.a - lo);
    return s.empty() ? "-" : s;
}

static std::string name(const piece& q) {
    return "s" + std::to_string(q.read >> 40) + "_" + std::to_string(q.read & 0xffffffffff) + "_" + std::to_string(q.k);
}

int cmd_select(const char* in, const char* dr, const char* rd, const char* out, int t) {
    gph g;
    adj a;
    keys x;
    sctx c;
    read_gph(in, &g);
    make_adj(&g, &a);
    c.cn = g.c.data();
    c.K = g.K;
    c.t = t;
    for (uint64_t s = 0; s < g.S; s++) {
        std::vector<uint32_t> smp;
        for (uint64_t v = 0; v < g.K; v += 997) if (full(&g, v)) smp.push_back(g.c[s * g.K + v]);
        std::sort(smp.begin(), smp.end());
        c.med.push_back(smp[smp.size() / 2]);
    }
    std::vector<uint8_t> drop(g.K);
    for (uint64_t h : read_set(dr)) drop[std::lower_bound(g.g.begin(), g.g.end(), h) - g.g.begin()] = 1;
    make_keys(&x, &g, &a, drop);
    g.ec = {};
    count_keys(&x, rd, t);
    c.x = &x;
    gather(&c, drop);
    c.h.resize(t);
    c.ids.resize(t);
    rpass p;
    each_read(&p, g.g, rd, t, anchor, &c);
    pieces(&c);
    for (std::string& s : split(rd)) { c.fs.push_back(fileno(std::fopen(s.c_str(), "rb"))); c.fx.push_back(fileno(std::fopen((s + ".idx").c_str(), "rb"))); }
    par(t, cut, &c);
    std::string o(out);
    FILE *f = std::fopen((o + ".fa").c_str(), "w"), *v = std::fopen((o + ".tsv").c_str(), "w");
    uint64_t nb = 0;
    for (const piece& q : c.pc) { std::fprintf(f, ">%s\n%s\n", name(q).c_str(), q.s.c_str()); nb += q.s.size(); }
    for (uint32_t i = 0; i < c.w.size(); i++) {
        const site& e = c.w[i];
        if (e.hi <= e.lo) continue;
        const piece& q = c.pc[e.pc];
        std::string n;
        for (uint64_t s = 0; s < x.S; s++) n += (s ? "," : "") + std::to_string(x.n[e.slot * x.S + s]);
        std::fprintf(v, "%s\t%u\t%u\t%u\t%s\t%lu\t%u\t%s\t%s\n", name(q).c_str(), e.a - q.lo, e.b - q.lo, i, n.c_str(), (unsigned long)e.key, e.s, list(e.l, q.lo).c_str(), list(e.r, q.lo).c_str());
    }
    std::fclose(f);
    std::fclose(v);
    std::fprintf(stderr, "[select] %zu sites on %zu reads; %zu pieces, %lu bp\n", c.w.size(), c.rd.size(), c.pc.size(), nb);
    return 0;
}
