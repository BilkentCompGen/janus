#include <algorithm>
#include <string>
#include <vector>
#include "core.h"
#include "io.h"
#include "prune.h"

static const uint32_t ML = 4, X = 20, NONE = ~0u;

struct graph { std::vector<uint64_t> g; std::vector<uint32_t> oo, on, io, in; std::vector<uint8_t> som; };

struct arm { std::vector<uint64_t> p; std::string s; uint8_t got; };

struct pctx { std::vector<rdcore> r; std::vector<std::vector<uint64_t>> h; std::vector<arm> a; std::vector<std::pair<uint64_t, uint32_t>> first; int w; };

static void load(const char* path, graph* x) {
    bin_header h;
    FILE* f = open_in(path, &h);
    uint64_t K = h.nreads, E = h.total;
    std::vector<uint16_t> c(K);
    std::vector<uint64_t> e(E);
    x->g.resize(K);
    std::fread(x->g.data(), 8, K, f);
    std::fread(c.data(), 2, K, f);
    std::fseek(f, sizeof(bin_header) + K * (8 + 2 * h.resv[0]), SEEK_SET);
    std::fread(e.data(), 8, E, f);
    std::fclose(f);
    x->som.resize(K);
    x->oo.assign(K + 1, 0);
    x->io.assign(K + 1, 0);
    for (uint64_t i = 0; i < K; i++) x->som[i] = !c[i];
    for (uint64_t v : e) { x->oo[(v >> 32) + 1]++; x->io[(v & 0xffffffff) + 1]++; }
    for (uint64_t i = 0; i < K; i++) { x->oo[i + 1] += x->oo[i]; x->io[i + 1] += x->io[i]; }
    x->on.resize(E);
    x->in.resize(E);
    std::vector<uint32_t> p(x->io.begin(), x->io.end() - 1);
    for (uint64_t i = 0; i < E; i++) { x->on[i] = e[i] & 0xffffffff; x->in[p[e[i] & 0xffffffff]++] = e[i] >> 32; }
}

static void tips(const graph* x, const std::vector<uint32_t>& fo, const std::vector<uint32_t>& fn, const std::vector<uint32_t>& bo, std::vector<uint64_t>& d) {
    std::vector<uint32_t> ch;
    for (uint64_t v = 0; v < x->g.size(); v++)
        for (uint32_t i = fo[v]; fo[v + 1] - fo[v] > 1 && i < fo[v + 1]; i++) {
            ch.clear();
            for (uint32_t u = fn[i]; ch.size() < ML && x->som[u] && bo[u + 1] - bo[u] == 1; u = fn[fo[u]]) {
                ch.push_back(u);
                if (fo[u + 1] - fo[u] != 1) break;
            }
            if (!ch.empty() && fo[ch.back() + 1] == fo[ch.back()]) for (uint32_t w : ch) d.push_back(x->g[w]);
        }
}

static uint32_t trace(const graph* x, uint32_t u, std::vector<uint32_t>& a) {
    for (uint32_t s = 0; s <= ML; s++) {
        a.push_back(u);
        if (x->oo[u + 1] - x->oo[u] != 1 || x->som[u] != x->som[a[0]]) return NONE;
        uint32_t v = x->on[x->oo[u]], d = x->io[v + 1] - x->io[v];
        if (d >= 2) return v;
        u = v;
    }
    return NONE;
}

static void add_arm(pctx* p, const graph* x, uint32_t s, const std::vector<uint32_t>& a, uint32_t m) {
    arm r = {{x->g[s]}, "", 0};
    for (uint32_t v : a) r.p.push_back(x->g[v]);
    r.p.push_back(x->g[m]);
    p->first.push_back({r.p[1], (uint32_t)p->a.size()});
    p->a.push_back(r);
}

static std::string bases(const rdcore* r, const uint8_t* rec, uint32_t lo, uint32_t hi, uint32_t rc) {
    const uint8_t* q = rec[3] >> 7 ? r->raw.data() : rec + 4;
    std::string s(hi - lo, 0);
    for (uint32_t i = lo; i < hi; i++) {
        uint32_t p = rc ? r->len - 1 - i : i;
        s[i - lo] = "ACGT"[(q[p >> 2] >> (6 - 2 * (p & 3)) & 3) ^ rc];
    }
    return s;
}

static void match(pctx* x, const rdcore* r, const uint8_t* rec, const uint64_t* h, uint32_t lo, uint32_t hi) {
    for (uint32_t i = lo + 1; i < hi; i++)
        for (auto it = std::lower_bound(x->first.begin(), x->first.end(), std::make_pair(h[i], 0u)); it != x->first.end() && it->first == h[i]; it++) {
            arm* a = &x->a[it->second];
            uint32_t n = a->p.size(), j = 0;
            while (j < n && i + j - 1 < hi && h[i + j - 1] == a->p[j]) j++;
            if (j < n || __atomic_exchange_n(&a->got, 1, __ATOMIC_RELAXED)) continue;
            a->s = bases(r, rec, r->a[i - 1], r->b[i + n - 2], i >= r->nf ? 3 : 0);
        }
}

static void pr_work(block* b) {
    pctx* x = (pctx*)b->ctx;
    rdcore* r = &x->r[b->id];
    std::vector<uint64_t>& h = x->h[b->id];
    for (size_t o = 0, o2 = 0; o < b->inlen[0]; o += svb_size(b->in[0] + o, b->inlen[0] - o), o2 += exc_size(b->in[1] + o2, b->inlen[1] - o2)) {
        const uint8_t* rec = (const uint8_t*)b->in[0] + o;
        collapse(r, rec);
        load_cores(r, (const uint8_t*)b->in[1] + o2, x->w);
        if (h.size() < r->n + r->k / 32 + 2) h.resize(2 * r->n + r->k / 16 + 4);
        for (uint32_t i = 0; i < r->n; i++) h[i] = core_hash(h.data() + r->n, core_cut(r, i, h.data() + r->n));
        match(x, r, rec, h.data(), 0, r->nf);
        match(x, r, rec, h.data(), r->nf, r->n);
    }
}

static uint32_t edit(const char* x, uint32_t n, const char* y, uint32_t m) {
    std::vector<uint32_t> p(m + 1), c(m + 1);
    for (uint32_t j = 0; j <= m; j++) p[j] = j;
    for (uint32_t i = 1; i <= n; i++) {
        c[0] = i;
        for (uint32_t j = 1; j <= m; j++) c[j] = std::min(std::min(p[j] + 1, c[j - 1] + 1), p[j - 1] + (x[i - 1] != y[j - 1]));
        p.swap(c);
    }
    return p[m];
}

static bool small(const std::string& t, const std::string& n) {
    uint32_t k = std::min(t.size(), n.size()), a = 0, b = 0;
    while (a < k && t[a] == n[a]) a++;
    while (b < k - a && t[t.size() - 1 - b] == n[n.size() - 1 - b]) b++;
    uint32_t dt = t.size() - a - b, dn = n.size() - a - b;
    return std::max(dt, dn) - std::min(dt, dn) < X && edit(t.data() + a, dt, n.data() + a, dn) < X;
}

int cmd_prune(const char* in, const char* rd, const char* out, int t) {
    graph x;
    pctx p;
    std::vector<uint64_t> d;
    std::vector<uint32_t> bt;
    load(in, &x);
    tips(&x, x.oo, x.on, x.io, d);
    tips(&x, x.io, x.in, x.oo, d);
    for (uint32_t s = 0; s < x.g.size(); s++) {
        uint32_t n = x.oo[s + 1] - x.oo[s];
        if (n < 2) continue;
        std::vector<std::vector<uint32_t>> a(n);
        std::vector<uint32_t> m(n);
        for (uint32_t i = 0; i < n; i++) m[i] = trace(&x, x.on[x.oo[s] + i], a[i]);
        for (uint32_t i = 0; i < n; i++)
            for (uint32_t j = i + 1; j < n; j++) {
                if (m[i] == NONE || m[i] != m[j] || x.som[a[i][0]] == x.som[a[j][0]]) continue;
                uint32_t k = x.som[a[i][0]] ? i : j;
                bt.push_back(p.a.size());
                add_arm(&p, &x, s, a[k], m[k]);
                add_arm(&p, &x, s, a[i + j - k], m[k]);
            }
    }
    x = graph();
    std::sort(p.first.begin(), p.first.end());
    p.r.resize(t);
    p.h.resize(t);
    FILE* fo = open_tmp(out);
    for (std::string& s : split(rd)) {
        bin_header h, e;
        FILE *fi = open_in(s.c_str(), &h), *fe = open_in((s + ".exc").c_str(), &e);
        p.w = e.resv[1];
        stream_run(fi, fe, fo, t, svb_size, exc_size, pr_work, &p);
    }
    for (uint32_t b : bt)
        if (!p.a[b].s.empty() && !p.a[b + 1].s.empty() && small(p.a[b].s, p.a[b + 1].s)) d.insert(d.end(), p.a[b].p.begin() + 1, p.a[b].p.end() - 1);
    std::sort(d.begin(), d.end());
    d.erase(std::unique(d.begin(), d.end()), d.end());
    std::fwrite(d.data(), 8, d.size(), fo);
    bin_header o = {{'S', 'V', 'B', 'P'}, 1, d.size(), bt.size(), {}};
    finish(fo, out, &o);
    return 0;
}
