#include <algorithm>
#include <array>
#include <string>
#include <vector>
#include "cmd.h"
#include "io.h"
#include "keys.h"

static const uint32_t ML = 4; // longest tip removed, in nodes; bubble arms are traced up to ML + 1 nodes

struct want { std::vector<uint32_t> p; bool v; uint64_t score; uint8_t lk; std::string s; };

struct pctx { gph o; adj a; std::vector<uint64_t> m; uint64_t full; int t, back; std::vector<std::vector<uint32_t>> tp; std::vector<std::vector<want>> bw; std::vector<want> w; std::vector<uint8_t> jd; std::vector<uint32_t> fo, fi; keys k; std::vector<std::vector<hit>> h; };

static void masks(pctx* x) {
    const gph* g = &x->o;
    x->m.assign(g->K, 0);
    for (uint64_t s = 0; s < g->S; s++) for (uint64_t v = 0; v < g->K; v++) x->m[v] |= (uint64_t)(g->c[s * g->K + v] > 0) << s;
    x->full = (1ull << g->S) - 1;
}

static void tips(void* a, int t) {
    pctx* x = (pctx*)a;
    const std::vector<uint32_t>& fo = x->back ? x->a.io : x->a.oo, &fn = x->back ? x->a.in : x->a.on, &bo = x->back ? x->a.oo : x->a.io;
    std::vector<uint32_t> ch;
    for (uint64_t v = x->o.K * t / x->t; v < x->o.K * (t + 1) / x->t; v++)
        for (uint32_t i = fo[v]; fo[v + 1] - fo[v] > 1 && i < fo[v + 1]; i++) {
            ch.clear();
            for (uint32_t u = fn[i]; ch.size() < ML && x->m[u] != x->full && x->m[u] == x->m[fn[i]] && bo[u + 1] - bo[u] == 1; u = fn[fo[u]]) {
                ch.push_back(u);
                if (fo[u + 1] - fo[u] != 1) break;
            }
            if (!ch.empty() && fo[ch.back() + 1] == fo[ch.back()]) x->tp[t].insert(x->tp[t].end(), ch.begin(), ch.end());
        }
}

static uint32_t trace(const pctx* x, uint32_t u, std::vector<uint32_t>& a) {
    for (uint32_t s = 0; s <= ML; s++) {
        a.push_back(u);
        if (x->a.oo[u + 1] - x->a.oo[u] != 1) return NONE;
        uint32_t v = x->a.on[x->a.oo[u]], d = x->a.io[v + 1] - x->a.io[v];
        if (d >= 2) return v;
        u = v;
    }
    return NONE;
}

static bool variant(const pctx* x, const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    for (uint32_t v : a) if (x->m[v] != x->m[a[0]]) return false;
    for (uint32_t v : b) if (x->m[v] == x->m[a[0]]) return false;
    return x->m[a[0]] != x->full;
}

static void add_arm(std::vector<want>& w, uint32_t s, const std::vector<uint32_t>& a, uint32_t m, bool v) {
    want r = {{s}, v, 0, 0, ""};
    r.p.insert(r.p.end(), a.begin(), a.end());
    r.p.push_back(m);
    w.push_back(r);
}

static void bubbles(void* p, int t) {
    pctx* x = (pctx*)p;
    for (uint64_t s = x->o.K * t / x->t; s < x->o.K * (t + 1) / x->t; s++) {
        uint32_t n = x->a.oo[s + 1] - x->a.oo[s];
        if (n < 2) continue;
        std::vector<std::vector<uint32_t>> a(n);
        std::vector<uint32_t> m(n);
        for (uint32_t i = 0; i < n; i++) m[i] = trace(x, x->a.on[x->a.oo[s] + i], a[i]);
        for (uint32_t i = 0; i < n; i++)
            for (uint32_t j = i + 1; j < n; j++) {
                if (m[i] == NONE || m[i] != m[j]) continue;
                bool vi = variant(x, a[i], a[j]), vj = variant(x, a[j], a[i]);
                if (!vi && !vj) continue;
                add_arm(x->bw[t], s, a[i], m[i], vi);
                add_arm(x->bw[t], s, a[j], m[i], vj);
            }
    }
}

static void index_arms(pctx* x) {
    const std::vector<want>& w = x->w;
    x->fo.assign(x->o.K + 1, 0);
    x->fi.resize(w.size());
    for (const want& a : w) x->fo[a.p[0] + 1]++;
    for (uint32_t i = 0; i < x->o.K; i++) x->fo[i + 1] += x->fo[i];
    std::vector<uint32_t> p(x->fo.begin(), x->fo.end() - 1);
    for (uint32_t i = 0; i < w.size(); i++) x->fi[p[w[i].p[0]]++] = i;
}

static void match(pctx* x, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc) {
    const rdcore* r = &c->r;
    for (uint32_t i = 0; i < m; i++)
        for (uint32_t q = x->fo[c->id[i]]; q < x->fo[c->id[i] + 1]; q++) {
            want* a = &x->w[x->fi[q]];
            uint32_t n = a->p.size(), j = 0;
            while (j < n && i + j < m && c->id[i + j] == a->p[j]) j++;
            if (j < n) continue;
            uint32_t st = r->a[c->ix[i]], en = r->b[c->ix[i + n - 1]];
            uint64_t sc = 1ull << 32 | (uint32_t)(r->hs[0] ^ st ^ r->len);
            if (sc <= __atomic_load_n(&a->score, __ATOMIC_RELAXED)) continue;
            while (__atomic_exchange_n(&a->lk, 1, __ATOMIC_ACQUIRE));
            if (sc > a->score) { a->score = sc; a->s = bases(r, rec, st, en, rc); }
            __atomic_store_n(&a->lk, 0, __ATOMIC_RELEASE);
        }
}

static void judge(void* p, int t) {
    pctx* x = (pctx*)p;
    const std::vector<want>& w = x->w;
    for (uint64_t b = 2 * t; b < w.size(); b += 2 * x->t) x->jd[b / 2] = !w[b].s.empty() && !w[b + 1].s.empty() && same(w[b].s, w[b + 1].s);
}

static std::vector<uint8_t> seeds(const gph* g, const std::vector<uint8_t>& drop) {
    std::vector<uint8_t> sd(g->K, 0);
    for (uint64_t i = 0; i < g->E; i++) {
        uint64_t u = g->e[i] >> 32, v = g->e[i] & 0xffffffff;
        if (g->nv[i] && !drop[u] && !drop[v]) sd[u] = sd[v] = 1;
    }
    return sd;
}

static void both(rpass* p, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) {
    pctx* x = (pctx*)p->ctx;
    match(x, c, rec, m, rc);
    tally(&x->k, c, m, p->f, rid, x->h[c - p->t.data()]);
}

static bool odd(const keys* x, uint64_t i) {
    bool on = false;
    for (uint64_t s = 1; s < x->S; s++) on |= x->n[i * x->S + s] >= MM;
    return on && !x->n[i * x->S];
}

static void gather(const keys* x, const std::vector<uint8_t>& drop, std::vector<uint64_t>& rows) {
    const gph* g = x->g;
    std::vector<std::array<uint64_t, 3>> w;
    for (uint64_t i = 0; i < x->Z; i++) {
        if (!x->key[i] || !odd(x, i) || !kept(x, x->key[i])) continue;
        uint64_t ei = (x->key[i] >> 32) - 1, eo = x->key[i] & 0xffffffff;
        if (drop[g->e[ei] >> 32] || drop[g->e[ei] & 0xffffffff] || drop[g->e[eo] >> 32] || drop[g->e[eo] & 0xffffffff]) continue;
        for (uint64_t s = 1; s < x->S; s++)
            if (x->n[i * x->S + s] >= MM) w.push_back({s << 40 | (0xffffffffu - (uint32_t)x->best[i * x->S + s]), x->key[i], i});
    }
    std::sort(w.begin(), w.end());
    for (const std::array<uint64_t, 3>& e : w) {
        rows.push_back(e[0]);
        rows.push_back(e[1]);
        rows.insert(rows.end(), x->n.begin() + e[2] * x->S, x->n.begin() + (e[2] + 1) * x->S);
    }
}

int cmd_prune(const char* in, const char* rd, const char* out, int t) {
    pctx x;
    std::vector<uint64_t> d, rows;
    std::vector<uint8_t> drop;
    read_gph(in, &x.o);
    make_adj(&x.o, &x.a);
    masks(&x);
    x.t = t;
    x.tp.resize(t);
    x.bw.resize(t);
    drop.assign(x.o.K, 0);
    for (x.back = 0; x.back < 2; x.back++) par(t, tips, &x);
    par(t, bubbles, &x);
    for (int i = 0; i < t; i++) {
        for (uint32_t v : x.tp[i]) drop[v] = 1;
        x.w.insert(x.w.end(), x.bw[i].begin(), x.bw[i].end());
    }
    x.tp = {};
    x.bw = {};
    x.m = {};
    make_keys(&x.k, &x.o, &x.a);
    make_table(&x.k, keep_keys(&x.k, seeds(&x.o, std::vector<uint8_t>(x.o.K, 0))));
    x.h.resize(t);
    index_arms(&x);
    rpass p;
    each_read(&p, x.o.g, rd, t, both, &x);
    x.jd.assign(x.w.size() / 2, 0);
    par(t, judge, &x);
    for (uint32_t b = 0; b < x.w.size(); b++)
        if (x.jd[b / 2] && x.w[b].v) for (uint32_t j = 1; j + 1 < x.w[b].p.size(); j++) drop[x.w[b].p[j]] = 1;
    for (uint64_t v = 0; v < x.o.K; v++) if (drop[v]) d.push_back(x.o.g[v]);
    std::sort(d.begin(), d.end());
    std::string o(out);
    bin_header hd = header('P', 1, d.size(), x.w.size() / 2), hs = header('S', 1, 0, x.o.S);
    save((o + ".drp").c_str(), &hd, d.data(), d.size() * 8);
    keep_keys(&x.k, seeds(&x.o, drop));
    gather(&x.k, drop, rows);
    hs.nreads = rows.size() / (2 + x.o.S);
    save((o + ".sites").c_str(), &hs, rows.data(), rows.size() * 8);
    std::fprintf(stderr, "[prune] %zu dropped nodes, %lu sites\n", d.size(), (unsigned long)hs.nreads);
    return 0;
}
