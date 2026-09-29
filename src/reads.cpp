#include <algorithm>
#include "io.h"
#include "reads.h"

static const uint32_t X = 20;

struct fctx { std::vector<want>* w; std::vector<uint32_t> fo, fi; };

std::string bases(const rdcore* r, const uint8_t* rec, uint32_t lo, uint32_t hi, uint32_t rc) {
    const uint8_t* q = rec[3] >> 7 ? r->raw.data() : rec + 4;
    std::string s(hi - lo, 0);
    for (uint32_t i = lo; i < hi; i++) {
        uint32_t p = rc ? r->len - 1 - i : i;
        s[i - lo] = "ACGT"[(q[p >> 2] >> (6 - 2 * (p & 3)) & 3) ^ rc];
    }
    return s;
}

static const uint32_t NONE = ~0u;

static void make_map(rpass* x) {
    uint64_t z = 1, n = x->g->size();
    while (z < n + n / 2) z <<= 1;
    x->tab.assign(z, {0, NONE});
    x->mask = z - 1;
    for (uint32_t v = 0; v < n; v++) {
        uint64_t h = (*x->g)[v], i = h & x->mask;
        while (x->tab[i].h) i = (i + 1) & x->mask;
        x->tab[i] = {h, v};
    }
}

static void hashed(rpass* x, rth* c, const uint8_t* rec, const uint8_t* exc) {
    rdcore* r = &c->r;
    collapse(r, rec);
    load_cores(r, exc, x->wb);
    if (c->h.size() < r->n + r->k / 32 + 2) c->h.resize(2 * r->n + r->k / 16 + 4);
    if (c->nd.size() < r->n) c->nd.resize(2 * r->n);
    uint64_t* w = c->h.data() + r->n;
    for (uint32_t i = 0; i < r->n; i++) {
        c->h[i] = core_hash(w, core_cut(r, i, w));
        __builtin_prefetch(&x->tab[c->h[i] & x->mask]);
    }
    for (uint32_t i = 0; i < r->n; i++) {
        uint64_t j = c->h[i] & x->mask;
        while (x->tab[j].h && x->tab[j].h != c->h[i]) j = (j + 1) & x->mask;
        c->nd[i] = x->tab[j].v;
    }
}

static void rp_work(block* b) {
    rpass* x = (rpass*)b->ctx;
    rth* c = &x->t[b->id];
    rdcore* r = &c->r;
    uint64_t rid = b->rec;
    for (size_t o = 0, o2 = 0; o < b->inlen[0]; o += svb_size(b->in[0] + o, b->inlen[0] - o), o2 += exc_size(b->in[1] + o2, b->inlen[1] - o2), rid++) {
        const uint8_t* rec = (const uint8_t*)b->in[0] + o;
        hashed(x, c, rec, (const uint8_t*)b->in[1] + o2);
        if (c->id.size() < r->n) { c->id.resize(2 * r->n); c->ix.resize(2 * r->n); }
        for (uint32_t s = 0; s < 2; s++) {
            uint32_t m = 0;
            for (uint32_t i = s ? r->nf : 0; i < (s ? r->n : r->nf); i++)
                if (c->nd[i] != NONE) { c->id[m] = c->nd[i]; c->ix[m++] = i; }
            x->fn(x, c, rec, m, s ? 3 : 0, rid);
        }
    }
}

void each_read(rpass* x, const std::vector<uint64_t>& g, const char* reads, int t, read_fn fn, void* ctx) {
    *x = {std::vector<rth>(t), &g, {}, 0, fn, ctx, 0, 0};
    make_map(x);
    FILE* fo = std::fopen("/dev/null", "w");
    for (std::string& s : split(reads)) {
        bin_header h, e;
        FILE *fi = open_in(s.c_str(), &h), *fe = open_in((s + ".exc").c_str(), &e);
        x->wb = e.resv[1];
        stream_run(fi, fe, fo, t, svb_size, exc_size, rp_work, x);
        x->f++;
    }
    std::fclose(fo);
}

static void match(rpass* x, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t) {
    fctx* f = (fctx*)x->ctx;
    const rdcore* r = &c->r;
    for (uint32_t i = 0; i < m; i++)
        for (uint32_t q = f->fo[c->id[i]]; q < f->fo[c->id[i] + 1]; q++) {
            want* a = &(*f->w)[f->fi[q]];
            uint32_t n = a->p.size(), j = 0;
            while (j < n && i + j < m && c->id[i + j] == a->p[j]) j++;
            if (j < n || !(a->mask >> x->f & 1)) continue;
            uint32_t st = r->a[c->ix[i]], en = r->b[c->ix[i + n - 1]], tg = std::max(a->pad, a->mul * (en - st)), lp = std::min(tg, st), rp = std::min(tg, r->len - en);
            uint64_t sc = (uint64_t)(lp + rp + 1) << 32 | (uint32_t)(r->hs[0] ^ st ^ r->len);
            if (sc <= __atomic_load_n(&a->score, __ATOMIC_RELAXED)) continue;
            while (__atomic_exchange_n(&a->lk, 1, __ATOMIC_ACQUIRE));
            if (sc > a->score) { a->score = sc; a->lp = lp; a->rp = rp; a->s = bases(r, rec, st - lp, en + rp, rc); }
            __atomic_store_n(&a->lk, 0, __ATOMIC_RELEASE);
        }
}

void fetch(std::vector<want>& w, const std::vector<uint64_t>& g, const char* reads, int t) {
    fctx f = {&w, std::vector<uint32_t>(g.size() + 1), std::vector<uint32_t>(w.size())};
    for (want& a : w) f.fo[a.p[0] + 1]++;
    for (uint32_t i = 0; i < g.size(); i++) f.fo[i + 1] += f.fo[i];
    std::vector<uint32_t> p(f.fo.begin(), f.fo.end() - 1);
    for (uint32_t i = 0; i < w.size(); i++) f.fi[p[w[i].p[0]]++] = i;
    rpass x;
    each_read(&x, g, reads, t, match, &f);
}

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

bool same(const std::string& t, const std::string& n) {
    uint32_t k = std::min(t.size(), n.size()), a = 0, b = 0;
    while (a < k && t[a] == n[a]) a++;
    while (b < k - a && t[t.size() - 1 - b] == n[n.size() - 1 - b]) b++;
    uint32_t dt = t.size() - a - b, dn = n.size() - a - b;
    return std::max(dt, dn) - std::min(dt, dn) < X && edit(t.data() + a, dt, n.data() + a, dn) < X;
}
