#include <algorithm>
#include <string>
#include <vector>
#include "files.h"
#include "graph.h"
#include "io.h"

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
    o.sh.assign(o.S * o.K, 0);
    for (uint64_t i = 0; i < o.S; i++) {
        std::vector<uint64_t> lo = read_set((ps[i] + ".low").c_str());
        for (uint64_t v = 0, k = 0; v < o.K; v++) {
            while (k < lo.size() && lo[k] < o.g[v]) k++;
            o.sh[i * o.K + v] = k < lo.size() && lo[k] == o.g[v];
        }
    }
    o.E = o.e.size();
    o.ec.resize(o.S * o.E);
    for (uint64_t i = 0; i < o.S; i++) {
        read_cnt(ps[i].c_str(), &s);
        place(o.g, &s, o.c.data() + i * o.K);
        for (uint64_t j = 0, k = 0; j < s.e.size(); j++) {
            while (o.e[k] < s.e[j]) k++;
            o.ec[i * o.E + k] = s.ec[j];
        }
    }
    write_gph(out, &o);
    return 0;
}
