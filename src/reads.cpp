#include <algorithm>
#include "io.h"
#include "reads.h"

static const uint32_t X = 20; // edits at which two sequences stop being the same allele

std::string bases(const rdcore* r, const uint8_t* rec, uint32_t lo, uint32_t hi, uint32_t rc) {
    const uint8_t* q = rec[3] >> 7 ? r->raw.data() : rec + 4;
    std::string s(hi - lo, 0);
    for (uint32_t i = lo; i < hi; i++) {
        uint32_t p = rc ? r->len - 1 - i : i;
        s[i - lo] = "ACGT"[(q[p >> 2] >> (6 - 2 * (p & 3)) & 3) ^ rc];
    }
    return s;
}

void make_map(rpass* x) {
    uint64_t z = 1, n = x->g->size();
    while (z < n + n / 2) z <<= 1;
    x->tab.assign(z, {0, NONE});
    x->mask = z - 1;
    for (uint32_t v = 0; v < n; v++) {
        if (v + PF < n) __builtin_prefetch(&x->tab[(*x->g)[v + PF] & x->mask]);
        uint64_t h = (*x->g)[v], i = h & x->mask;
        while (x->tab[i].h) i = (i + 1) & x->mask;
        x->tab[i] = {h, v};
    }
}

void lookup(const rpass* x, rth* c) {
    const rdcore* r = &c->r;
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

void chains(rpass* x, rth* c, const uint8_t* rec, uint64_t rid) {
    const rdcore* r = &c->r;
    if (c->id.size() < r->n) { c->id.resize(2 * r->n); c->ix.resize(2 * r->n); }
    for (uint32_t s = 0; s < 2; s++) {
        uint32_t m = 0;
        for (uint32_t i = s ? r->nf : 0; i < (s ? r->n : r->nf); i++)
            if (c->nd[i] != NONE) { c->id[m] = c->nd[i]; c->ix[m++] = i; }
        x->fn(x, c, rec, m, s ? 3 : 0, rid);
    }
}

static void rp_work(block* b) {
    rpass* x = (rpass*)b->ctx;
    rth* c = &x->t[b->id];
    uint64_t rid = b->rec;
    for (size_t o = 0, e = 0; o < b->inlen[0]; rid++) {
        const uint8_t* rec = next_read(&c->r, b, &o, &e, x->wb);
        lookup(x, c);
        chains(x, c, rec, rid);
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
        stream_run(fi, fe, fo, t, rp_work, x);
        x->f++;
    }
    std::fclose(fo);
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
