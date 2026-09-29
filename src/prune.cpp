#include <algorithm>
#include <string>
#include <vector>
#include "files.h"
#include "cmd.h"
#include "reads.h"

static const uint32_t ML = 4, NONE = ~0u;

struct graph { std::vector<uint64_t> g, m; adj a; uint64_t full; };

static void load(const char* path, graph* x) {
    gph o;
    read_gph(path, &o);
    make_adj(&o, &x->a);
    masks(&o, x->m);
    x->full = (1ull << o.S) - 1;
    x->g.swap(o.g);
}

static void tips(const graph* x, const std::vector<uint32_t>& fo, const std::vector<uint32_t>& fn, const std::vector<uint32_t>& bo, std::vector<uint64_t>& d) {
    std::vector<uint32_t> ch;
    for (uint64_t v = 0; v < x->g.size(); v++)
        for (uint32_t i = fo[v]; fo[v + 1] - fo[v] > 1 && i < fo[v + 1]; i++) {
            ch.clear();
            for (uint32_t u = fn[i]; ch.size() < ML && x->m[u] != x->full && x->m[u] == x->m[fn[i]] && bo[u + 1] - bo[u] == 1; u = fn[fo[u]]) {
                ch.push_back(u);
                if (fo[u + 1] - fo[u] != 1) break;
            }
            if (!ch.empty() && fo[ch.back() + 1] == fo[ch.back()]) for (uint32_t w : ch) d.push_back(x->g[w]);
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

int cmd_prune(const char* in, const char* rd, const char* out, int t) {
    graph x;
    std::vector<want> w;
    std::vector<uint64_t> d;
    std::vector<uint8_t> vf;
    load(in, &x);
    tips(&x, x.a.oo, x.a.on, x.a.io, d);
    tips(&x, x.a.io, x.a.in, x.a.oo, d);
    for (uint32_t s = 0; s < x.g.size(); s++) {
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
    x.a = adj();
    x.m = {};
    fetch(w, x.g, rd, t);
    for (uint32_t b = 0; b < w.size(); b += 2)
        if (!w[b].s.empty() && !w[b + 1].s.empty() && same(w[b].s, w[b + 1].s))
            for (uint32_t u = 0; u < 2; u++) if (vf[b / 2] >> u & 1) for (uint32_t j = 1; j + 1 < w[b + u].p.size(); j++) d.push_back(x.g[w[b + u].p[j]]);
    std::sort(d.begin(), d.end());
    d.erase(std::unique(d.begin(), d.end()), d.end());
    write_set(out, "SVBP", d, w.size() / 2);
    return 0;
}
