#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <string>
#include "cmd.h"
#include "extract.h"
#include "io.h"
#include "keys.h"

static const uint32_t PAD = 1000; // least flank cut on each side of a site's anchors
static const uint32_t FL = 3; // the flank is at least FL times the site's length (innermost anchor to innermost anchor)
static const uint32_t NA = 3; // anchors kept on each side of a site
static const uint32_t MULT = 2; // an anchor's count may be at most MULT times its sample's median node count
static const double DL = 15000; // an anchor's weight falls linearly to zero at DL bases from the combination

struct anc { double sc; uint32_t v, a, b; };

struct site { uint64_t read, key; uint32_t s, lo, hi, a, b, pc; std::string n; std::vector<anc> l, r; };

struct piece { uint64_t read; uint32_t lo, hi, k; std::string s; };

struct sctx { keys* x; rpass* p; std::vector<site> w; std::vector<piece> pc; std::vector<uint64_t> rd; std::vector<uint32_t> ro, med, gs; const uint16_t* cn; uint64_t K, lo, hi; std::vector<std::vector<hit>> h; std::vector<std::vector<uint32_t>> ids; std::vector<scr> sc; std::vector<int> fs, fx; int t; };

static bool by_cost(const anc& p, const anc& q) { return p.sc < q.sc; }

static void load(sctx* c, const char* path) {
    bin_header h;
    FILE* f = open_in(path, &h);
    std::vector<uint64_t> r(2 + h.total);
    for (uint64_t i = 0; i < h.nreads; i++) {
        std::fread(r.data(), 8, r.size(), f);
        site w = {};
        w.read = r[0];
        w.key = r[1];
        w.s = r[0] >> 40;
        for (uint64_t s = 0; s < h.total; s++) w.n += (s ? "," : "") + std::to_string(r[2 + s]);
        c->w.push_back(w);
    }
    std::fclose(f);
    for (size_t i = 0; i < c->w.size(); i++)
        if (!i || c->w[i].read != c->w[i - 1].read) { c->rd.push_back(c->w[i].read); c->ro.push_back(i); }
    c->ro.push_back(c->w.size());
}

static void medians(sctx* c, const gph* g) {
    std::vector<uint32_t> fv, smp;
    for (uint64_t v = 0; v < g->K; v++) {
        bool full = true;
        for (uint64_t s = 0; s < g->S; s++) full &= g->c[s * g->K + v] > 0;
        if (full) fv.push_back(v);
    }
    for (uint64_t s = 0; s < g->S; s++) {
        smp.clear();
        for (uint32_t v : fv) smp.push_back(g->c[s * g->K + v]);
        std::nth_element(smp.begin(), smp.begin() + smp.size() / 2, smp.end());
        c->med.push_back(smp[smp.size() / 2]);
    }
}

static double fit(const sctx* c, const std::vector<uint32_t>& ids, uint32_t v, uint64_t s) {
    auto r = std::equal_range(ids.begin(), ids.end(), v);
    if (!c->cn[v] || !c->cn[s * c->K + v] || r.second - r.first != 1 || c->cn[v] > MULT * c->med[0] || c->cn[s * c->K + v] > MULT * c->med[s]) return 0;
    double n = (double)c->cn[v] / c->med[0], t = (double)c->cn[s * c->K + v] / c->med[s];
    return std::sqrt(std::min(1.0, n) * std::min(1.0, t));
}

static void side(const sctx* c, const rth* t, const std::vector<uint32_t>& ids, uint64_t s, int64_t from, int64_t step, int64_t end, std::vector<anc>& o) {
    const rdcore* d = &t->r;
    uint32_t at = d->a[t->ix[from]];
    for (int64_t j = from; j != end; j += step) {
        uint32_t p = d->a[t->ix[j]];
        double w = std::max(0.0, 1 - (p > at ? p - at : at - p) / DL), f;
        if (w <= 0 || (o.size() == NA && w <= -o.back().sc)) break;
        if ((f = fit(c, ids, t->id[j], s)) == 0) continue;
        o.push_back({-(f * w), t->id[j], p, d->b[t->ix[j]]});
        std::sort(o.begin(), o.end(), by_cost);
        if (o.size() > NA) o.pop_back();
    }
    for (anc& x : o) x.sc = x.a > at ? x.a - at : at - x.a;
    std::sort(o.begin(), o.end(), by_cost);
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
            uint32_t pd = std::max(PAD, FL * (e->b - e->a)), lo = e->l.empty() ? e->a : e->l.back().a, hi = e->r.empty() ? e->b : e->r.back().b;
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

static void fetch(const sctx* c, uint64_t rd, std::vector<uint8_t>& rec) {
    uint64_t off;
    uint32_t tag, s = rd >> 40;
    pread(c->fx[s], &off, 8, sizeof(bin_header) + 8 * (rd & 0xffffffffff));
    pread(c->fs[s], &tag, 4, off);
    rec.assign(svb_bytes(tag) + 8, 0);
    pread(c->fs[s], rec.data(), rec.size() - 8, off);
}

static void grab(void* a, int i) {
    sctx* c = (sctx*)a;
    rth* t = &c->p->t[i];
    std::vector<uint8_t> rec;
    for (uint64_t k = c->lo + i; k < c->hi; k += c->t) {
        fetch(c, c->rd[k], rec);
        cores(&c->sc[i], &t->r, rec.data(), LEV);
        lookup(c->p, t);
        chains(c->p, t, rec.data(), c->rd[k] & 0xffffffffff);
    }
}

static void cut(void* a, int i) {
    sctx* c = (sctx*)a;
    rdcore r;
    std::vector<uint8_t> rec;
    for (size_t k = i; k + 1 < c->gs.size(); k += c->t) {
        fetch(c, c->pc[c->gs[k]].read, rec);
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
    LCP_INIT();
    gph g;
    adj a;
    keys x;
    sctx c;
    read_gph(in, &g);
    make_adj(&g, &a);
    c.cn = g.c.data();
    c.K = g.K;
    c.t = t;
    medians(&c, &g);
    make_keys(&x, &g, &a);
    keep_keys(&x, std::vector<uint8_t>(g.K, 1));
    c.x = &x;
    load(&c, dr);
    c.h.resize(t);
    c.ids.resize(t);
    for (std::string& f : split(rd)) { c.fs.push_back(fileno(std::fopen(f.c_str(), "rb"))); c.fx.push_back(fileno(std::fopen((f + ".idx").c_str(), "rb"))); }
    rpass p = {std::vector<rth>(t), &g.g, {}, 0, anchor, &c, 0, 0};
    make_map(&p);
    c.p = &p;
    c.sc.resize(t);
    for (uint64_t k = 0; k < c.rd.size(); p.f++) {
        c.lo = k;
        while (k < c.rd.size() && c.rd[k] >> 40 == (uint64_t)p.f) k++;
        c.hi = k;
        par(t, grab, &c);
    }
    pieces(&c);
    par(t, cut, &c);
    std::string o(out);
    FILE *f = std::fopen((o + ".fa").c_str(), "w"), *v = std::fopen((o + ".tsv").c_str(), "w");
    uint64_t nb = 0;
    for (const piece& q : c.pc) { std::fprintf(f, ">%s\n%s\n", name(q).c_str(), q.s.c_str()); nb += q.s.size(); }
    for (uint32_t i = 0; i < c.w.size(); i++) {
        const site& e = c.w[i];
        if (e.hi <= e.lo) continue;
        const piece& q = c.pc[e.pc];
        std::fprintf(v, "%s\t%u\t%u\t%u\t%s\t%lu\t%u\t%s\t%s\n", name(q).c_str(), e.a - q.lo, e.b - q.lo, i, e.n.c_str(), (unsigned long)e.key, e.s, list(e.l, q.lo).c_str(), list(e.r, q.lo).c_str());
    }
    std::fclose(f);
    std::fclose(v);
    std::fprintf(stderr, "[select] %zu sites on %zu reads; %zu pieces, %lu bp\n", c.w.size(), c.rd.size(), c.pc.size(), nb);
    return 0;
}
