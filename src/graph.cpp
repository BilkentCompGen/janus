#include <algorithm>
#include <string>
#include <vector>
#include "files.h"
#include "cmd.h"
#include "io.h"

static const uint32_t MM = 3, NONE = ~0u, PF = 16, WK = 16;

static void unite(std::vector<uint64_t>& a, const std::vector<uint64_t>& b, std::vector<uint64_t>& t) {
    t.resize(a.size() + b.size());
    t.resize(std::set_union(a.begin(), a.end(), b.begin(), b.end(), t.begin()) - t.begin());
    a.swap(t);
}

static void place(const std::vector<uint64_t>& g, cnt* s, uint16_t* c) {
    std::vector<uint32_t> id(s->hs.size());
    for (uint64_t j = 0, k = 0; j < s->hs.size(); j++) {
        while (g[k] < s->hs[j]) k++;
        id[j] = k;
        c[k] = s->n[j];
    }
    for (uint64_t& v : s->e) v = (uint64_t)id[v >> 32] << 32 | id[v & 0xffffffff];
}

// Segments are the maximal non-branching paths. An edge u -> v joins them when u has one successor and v one
// predecessor; a segment starts at a node without a joining in-edge (a later pass takes the pure cycles). WK walks
// run side by side so their cache misses overlap; a segment gets consecutive new ids when its walk ends. The other
// passes scan the hash-sorted edges in order and prefetch the one random access each step makes.
static void segments(const gph* o, std::vector<uint32_t>& id, std::vector<uint32_t>& ord, std::vector<uint32_t>& st) {
    const uint64_t K = o->K, E = o->E, M = 0xffffffff;
    const std::vector<uint64_t>& e = o->e;
    std::vector<uint32_t> oo(K + 1, 0), in(K, 0), nx(K, NONE);
    for (uint64_t i = 0; i < E; i++) oo[(e[i] >> 32) + 1]++;
    for (uint64_t u = 0; u < K; u++) oo[u + 1] += oo[u];
    for (uint64_t i = 0; i < E; i++) {
        if (i + PF < E) __builtin_prefetch(&in[e[i + PF] & M]);
        in[e[i] & M]++;
    }
    for (uint64_t u = 0; u < K; u++) {
        if (u + PF < K && oo[u + PF] < E) __builtin_prefetch(&in[e[oo[u + PF]] & M]);
        uint64_t v = e[oo[u]] & M;
        if (oo[u + 1] - oo[u] == 1 && v != u && in[v] == 1) { nx[u] = v; in[v] = NONE; }
    }
    id.assign(K, NONE);
    std::vector<uint32_t> at(WK, NONE);
    std::vector<std::vector<uint32_t>> b(WK);
    for (uint64_t h = 0, busy = 1; busy;) {
        busy = 0;
        for (uint32_t j = 0; j < WK; j++) {
            if (at[j] == NONE) {
                while (h < K && in[h] == NONE) h++;
                if (h == K) continue;
                at[j] = h++;
            }
            busy = 1;
            b[j].push_back(at[j]);
            at[j] = nx[at[j]];
            if (at[j] != NONE) { __builtin_prefetch(&nx[at[j]]); continue; }
            st.push_back(ord.size());
            for (uint32_t w : b[j]) { id[w] = ord.size(); ord.push_back(w); }
            b[j].clear();
        }
    }
    for (uint64_t v = 0; v < K; v++) {
        if (id[v] != NONE) continue;
        st.push_back(ord.size());
        for (uint32_t w = v; id[w] == NONE; w = nx[w]) { id[w] = ord.size(); ord.push_back(w); }
    }
    st.push_back(K);
}

// Renumber nodes and edges to the segment order; edges stay sorted by (source, target) in the new ids.
static void renumber(gph* o, const std::vector<uint32_t>& id, const std::vector<uint32_t>& ord) {
    const uint64_t K = o->K, E = o->E, M = 0xffffffff;
    std::vector<uint64_t> g(K), e(E);
    std::vector<uint16_t> c(o->S * K);
    std::vector<uint8_t> nv(E);
    std::vector<uint32_t> at(K + 1, 0);
    for (uint64_t i = 0; i < E; i++) at[id[o->e[i] >> 32] + 1]++;
    for (uint64_t v = 0; v < K; v++) {
        if (v + PF < K) __builtin_prefetch(&o->g[ord[v + PF]]);
        g[v] = o->g[ord[v]];
        at[v + 1] += at[v];
    }
    for (uint64_t s = 0; s < o->S; s++)
        for (uint64_t v = 0; v < K; v++) {
            if (v + PF < K) __builtin_prefetch(&o->c[s * K + ord[v + PF]]);
            c[s * K + v] = o->c[s * K + ord[v]];
        }
    for (uint64_t i = 0; i < E; i++) {
        if (i + PF < E) { __builtin_prefetch(&id[o->e[i + PF] & M]); __builtin_prefetch(&at[id[o->e[i + PF] >> 32]]); }
        uint64_t u = id[o->e[i] >> 32], p = at[u]++;
        e[p] = u << 32 | id[o->e[i] & M];
        nv[p] = o->nv[i];
    }
    for (uint64_t p = 0, q; p < E; p = q) {
        for (q = p + 1; q < E && e[q] >> 32 == e[p] >> 32; q++);
        for (uint64_t i = p + 1; i < q; i++)
            for (uint64_t j = i; j > p && e[j] < e[j - 1]; j--) { std::swap(e[j], e[j - 1]); std::swap(nv[j], nv[j - 1]); }
    }
    o->g.swap(g);
    o->c.swap(c);
    o->e.swap(e);
    o->nv.swap(nv);
}

int cmd_graph(const char* in, const char* out) {
    std::vector<std::string> ps = split(in);
    gph o = {0, 0, ps.size(), {}, {}, {}, {}, {}};
    cnt s;
    std::vector<uint64_t> t;
    for (auto& p : ps) {
        read_cnt(p.c_str(), &s);
        unite(o.g, s.hs, t);
    }
    o.K = o.g.size();
    o.c.resize(o.S * o.K);
    for (uint64_t i = 0; i < o.S; i++) {
        read_cnt(ps[i].c_str(), &s);
        place(o.g, &s, o.c.data() + i * o.K);
        unite(o.e, s.e, t);
    }
    t = {};
    o.E = o.e.size();
    std::vector<uint8_t> tu(o.E, 0);
    o.nv.assign(o.E, 1);
    for (uint64_t i = 0; i < o.S; i++) {
        read_cnt(ps[i].c_str(), &s);
        place(o.g, &s, o.c.data() + i * o.K);
        for (uint64_t j = 0, k = 0; j < s.e.size(); j++) {
            while (o.e[k] < s.e[j]) k++;
            if (!i) o.nv[k] = !s.ec[j];
            else tu[k] |= s.ec[j] >= MM;
        }
    }
    s = {};
    for (uint64_t i = 0; i < o.E; i++) o.nv[i] &= tu[i];
    tu = {};
    std::vector<uint32_t> id, ord;
    segments(&o, id, ord, o.st);
    renumber(&o, id, ord);
    write_gph(out, &o);
    return 0;
}
