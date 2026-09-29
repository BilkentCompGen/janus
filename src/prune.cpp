#include <algorithm>
#include <string>
#include <vector>
#include "cmd.h"
#include "io.h"
#include "keys.h"

static const uint32_t ML = 4, MM = 3, ZF = 16, NONE = ~0u;

struct graph { gph o; std::vector<uint64_t> m; adj a; uint64_t full; };

struct pctx { fctx f; keys k; std::vector<std::vector<hit>> h; };

struct found { uint64_t read, key; std::vector<uint64_t> n; };

static bool by_read(const found& p, const found& q) { return p.read != q.read ? p.read < q.read : p.key < q.key; }

static void tips(const graph* x, const std::vector<uint32_t>& fo, const std::vector<uint32_t>& fn, const std::vector<uint32_t>& bo, std::vector<uint64_t>& d) {
    std::vector<uint32_t> ch;
    for (uint64_t v = 0; v < x->o.K; v++)
        for (uint32_t i = fo[v]; fo[v + 1] - fo[v] > 1 && i < fo[v + 1]; i++) {
            ch.clear();
            for (uint32_t u = fn[i]; ch.size() < ML && x->m[u] != x->full && x->m[u] == x->m[fn[i]] && bo[u + 1] - bo[u] == 1; u = fn[fo[u]]) {
                ch.push_back(u);
                if (fo[u + 1] - fo[u] != 1) break;
            }
            if (!ch.empty() && fo[ch.back() + 1] == fo[ch.back()]) for (uint32_t w : ch) d.push_back(x->o.g[w]);
        }
}

static uint32_t trace(const graph* x, uint32_t u, std::vector<uint32_t>& a) {
    for (uint32_t s = 0; s <= ML; s++) {
        a.push_back(u);
        if (x->a.oo[u + 1] - x->a.oo[u] != 1) return NONE;
        uint32_t v = x->a.on[x->a.oo[u]], d = x->a.io[v + 1] - x->a.io[v];
        if (d >= 2) return v;
        u = v;
    }
    return NONE;
}

static bool variant(const graph* x, const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    for (uint32_t v : a) if (x->m[v] != x->m[a[0]]) return false;
    for (uint32_t v : b) if (x->m[v] == x->m[a[0]]) return false;
    return x->m[a[0]] != x->full;
}

static void add_arm(std::vector<want>& w, uint32_t s, const std::vector<uint32_t>& a, uint32_t m) {
    want r = {{s}, ~0ull, 0, 0, 0, 0, 0, 0, ""};
    r.p.insert(r.p.end(), a.begin(), a.end());
    r.p.push_back(m);
    w.push_back(r);
}

static std::vector<uint8_t> seeds(const gph* g, const std::vector<uint8_t>& drop) {
    std::vector<uint8_t> sd(g->K, 0);
    for (uint64_t i = 0; i < g->E; i++) {
        uint64_t u = g->e[i] >> 32, v = g->e[i] & 0xffffffff;
        bool t = false;
        for (uint64_t s = 1; s < g->S; s++) t |= g->ec[s * g->E + i] >= MM;
        if (t && !g->ec[i] && !drop[u] && !drop[v]) sd[u] = sd[v] = 1;
    }
    return sd;
}

static void both(rpass* p, rth* c, const uint8_t* rec, uint32_t m, uint32_t rc, uint64_t rid) {
    pctx* q = (pctx*)p->ctx;
    match(&q->f, c, rec, m, rc, p->f);
    tally(&q->k, c, m, p->f, rid, q->h[c - p->t.data()]);
}

static bool odd(const keys* x, uint64_t i) {
    bool on = false;
    for (uint64_t s = 1; s < x->S; s++) on |= x->n[i * x->S + s] >= MM;
    return on && !x->n[i * x->S];
}

static void gather(const keys* x, const std::vector<uint8_t>& drop, std::vector<found>& w) {
    const gph* g = x->g;
    for (uint64_t i = 0; i < x->Z; i++) {
        if (!x->k[i] || !odd(x, i) || !kept(x, x->k[i])) continue;
        uint64_t ei = (x->k[i] >> 32) - 1, eo = x->k[i] & 0xffffffff;
        if (drop[g->e[ei] >> 32] || drop[g->e[ei] & 0xffffffff] || drop[g->e[eo] >> 32] || drop[g->e[eo] & 0xffffffff]) continue;
        for (uint64_t s = 1; s < x->S; s++)
            if (x->n[i * x->S + s] >= MM) w.push_back({s << 40 | (0xffffffffu - (uint32_t)x->best[i * x->S + s]), x->k[i], std::vector<uint64_t>(x->n.begin() + i * x->S, x->n.begin() + (i + 1) * x->S)});
    }
    std::sort(w.begin(), w.end(), by_read);
}

int cmd_prune(const char* in, const char* rd, const char* out, int t) {
    graph x;
    pctx q;
    std::vector<want> w;
    std::vector<uint64_t> d;
    std::vector<uint8_t> vf, drop;
    read_gph(in, &x.o);
    make_adj(&x.o, &x.a);
    masks(&x.o, x.m);
    x.full = (1ull << x.o.S) - 1;
    tips(&x, x.a.oo, x.a.on, x.a.io, d);
    tips(&x, x.a.io, x.a.in, x.a.oo, d);
    for (uint32_t s = 0; s < x.o.K; s++) {
        uint32_t n = x.a.oo[s + 1] - x.a.oo[s];
        if (n < 2) continue;
        std::vector<std::vector<uint32_t>> a(n);
        std::vector<uint32_t> m(n);
        for (uint32_t i = 0; i < n; i++) m[i] = trace(&x, x.a.on[x.a.oo[s] + i], a[i]);
        for (uint32_t i = 0; i < n; i++)
            for (uint32_t j = i + 1; j < n; j++) {
                if (m[i] == NONE || m[i] != m[j]) continue;
                uint8_t f = variant(&x, a[i], a[j]) | variant(&x, a[j], a[i]) << 1;
                if (!f) continue;
                vf.push_back(f);
                add_arm(w, s, a[i], m[i]);
                add_arm(w, s, a[j], m[i]);
            }
    }
    x.m = {};
    make_keys(&q.k, &x.o, &x.a);
    uint64_t np = keep_keys(&q.k, seeds(&x.o, std::vector<uint8_t>(x.o.K, 0)));
    for (q.k.Z = 1; q.k.Z < ZF * np; q.k.Z <<= 1);
    q.k.k.assign(q.k.Z, 0);
    q.k.n.assign(q.k.Z * q.k.S, 0);
    q.k.best.assign(q.k.Z * q.k.S, 0);
    q.h.resize(t);
    wanted(&q.f, w, x.o.K);
    rpass p;
    each_read(&p, x.o.g, rd, nullptr, t, both, &q);
    for (uint32_t b = 0; b < w.size(); b += 2)
        if (!w[b].s.empty() && !w[b + 1].s.empty() && same(w[b].s, w[b + 1].s))
            for (uint32_t u = 0; u < 2; u++) if (vf[b / 2] >> u & 1) for (uint32_t j = 1; j + 1 < w[b + u].p.size(); j++) d.push_back(x.o.g[w[b + u].p[j]]);
    std::sort(d.begin(), d.end());
    d.erase(std::unique(d.begin(), d.end()), d.end());
    std::string o(out);
    write_set((o + ".drp").c_str(), "SVBP", d, w.size() / 2);
    drop.assign(x.o.K, 0);
    for (uint64_t h : d) drop[std::lower_bound(x.o.g.begin(), x.o.g.end(), h) - x.o.g.begin()] = 1;
    keep_keys(&q.k, seeds(&x.o, drop));
    std::vector<found> st;
    gather(&q.k, drop, st);
    FILE* fo = open_tmp((o + ".sites").c_str());
    std::fseek(fo, sizeof(bin_header), SEEK_SET);
    for (const found& e : st) { std::fwrite(&e.read, 8, 1, fo); std::fwrite(&e.key, 8, 1, fo); std::fwrite(e.n.data(), 8, x.o.S, fo); }
    bin_header ho = {{'S', 'V', 'B', 'S'}, 1, st.size(), x.o.S, {}};
    finish(fo, (o + ".sites").c_str(), &ho);
    std::fprintf(stderr, "[prune] %zu dropped nodes, %zu sites\n", d.size(), st.size());
    return 0;
}
