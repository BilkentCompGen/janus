#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>
#include "cmd.h"
#include "files.h"
#include "io.h"

class Counts {
public:
    std::vector<uint64_t> hashes, edges;
    std::vector<uint16_t> reads, edge_reads;
    void read(const char* path, bool all);
    void to_union(const std::vector<uint64_t>& nodes, uint16_t* counts);
};

void Counts::read(const char* path, bool all) {
    Input in(path);
    uint64_t K = in.h.n, E = in.h.total;
    hashes.resize(K);
    std::fread(hashes.data(), 8, K, in.f);
    if (!all) { std::fclose(in.f); return; }
    reads.resize(K);
    edge_reads.resize(E);
    edges.clear();
    edges.reserve(E);
    std::fread(reads.data(), 2, K, in.f);
    std::fread(edge_reads.data(), 2, E, in.f);
    long at = std::ftell(in.f);
    std::fseek(in.f, 0, SEEK_END);
    std::vector<uint8_t> bb(std::ftell(in.f) - at + 8);
    std::fseek(in.f, at, SEEK_SET);
    std::fread(bb.data(), 1, bb.size() - 8, in.f);
    std::fclose(in.f);
    BitReader degree(bb.data(), 0), target(bb.data(), K * in.h.resv[1]);
    for (uint64_t i = 0; i < K; i++)
        for (uint32_t d = degree.get(in.h.resv[1]); d; d--) edges.push_back(i << 32 | target.get(in.h.resv[0]));
}

void Counts::to_union(const std::vector<uint64_t>& nodes, uint16_t* counts) {
    std::vector<uint32_t> id(hashes.size());
    for (uint64_t j = 0, k = 0; j < hashes.size(); j++) {
        while (nodes[k] < hashes[j]) k++;
        id[j] = k;
        counts[k] = reads[j];
    }
    for (uint64_t& v : edges) v = (uint64_t)id[v >> 32] << 32 | id[v & 0xffffffff];
}

class Segments {
public:
    std::vector<uint32_t> id, ord, starts;
    explicit Segments(const Graph& o);
};

static const uint32_t WK = 16;

Segments::Segments(const Graph& o) {
    const uint64_t K = o.K, E = o.E;
    const std::vector<uint64_t>& e = o.edges;
    std::vector<uint32_t> oo(K + 1, 0), in(K, 0), nx(K, NONE);
    std::vector<uint32_t>& st = starts;
    for (uint64_t i = 0; i < E; i++) oo[Graph::from(e[i]) + 1]++;
    for (uint64_t u = 0; u < K; u++) oo[u + 1] += oo[u];
    for (uint64_t i = 0; i < E; i++) {
        if (i + PF < E) __builtin_prefetch(&in[Graph::to(e[i + PF])]);
        in[Graph::to(e[i])]++;
    }
    for (uint64_t u = 0; u < K; u++) {
        if (u + PF < K && oo[u + PF] < E) __builtin_prefetch(&in[Graph::to(e[oo[u + PF]])]);
        if (oo[u + 1] - oo[u] != 1) continue;
        uint64_t v = Graph::to(e[oo[u]]);
        if (v != u && in[v] == 1) { nx[u] = v; in[v] = NONE; }
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

class MakeGraph {
public:
    explicit MakeGraph(const Args& a) : a(a) {}
    int run();

private:
    const Args& a;
    Graph o;
    void renumber(const std::vector<uint32_t>& id, const std::vector<uint32_t>& ord);
};

static void unite(std::vector<uint64_t>& a, const std::vector<uint64_t>& b, std::vector<uint64_t>& t) {
    t.resize(a.size() + b.size());
    t.resize(std::set_union(a.begin(), a.end(), b.begin(), b.end(), t.begin()) - t.begin());
    a.swap(t);
}

int MakeGraph::run() {
    std::vector<std::string> ps = split(a.in);
    std::vector<uint64_t> t;
    Counts s;
    o.S = ps.size();
    for (std::string& p : ps) {
        s.read(p.c_str(), false);
        unite(o.hashes, s.hashes, t);
    }
    o.K = o.hashes.size();
    o.counts.resize(o.S * o.K);
    for (uint64_t i = 0; i < o.S; i++) {
        s.read(ps[i].c_str(), true);
        s.to_union(o.hashes, o.counts.data() + i * o.K);
        unite(o.edges, s.edges, t);
    }
    t.clear();
    t.shrink_to_fit();
    o.E = o.edges.size();
    std::vector<uint8_t> tu(o.E, 0);
    o.novel.assign(o.E, 1);
    for (uint64_t i = 0; i < o.S; i++) {
        s.read(ps[i].c_str(), true);
        s.to_union(o.hashes, o.counts.data() + i * o.K);
        for (uint64_t j = 0, k = 0; j < s.edges.size(); j++) {
            while (o.edges[k] < s.edges[j]) k++;
            if (!i) o.novel[k] = !s.edge_reads[j];
            else tu[k] |= s.edge_reads[j] >= a.support;
        }
    }
    s = Counts();
    for (uint64_t i = 0; i < o.E; i++) o.novel[i] &= tu[i];
    tu.clear();
    tu.shrink_to_fit();
    Segments seg(o);
    o.starts.swap(seg.starts);
    renumber(seg.id, seg.ord);
    o.write(a.out);
    return 0;
}

void MakeGraph::renumber(const std::vector<uint32_t>& id, const std::vector<uint32_t>& ord) {
    const uint64_t K = o.K, E = o.E;
    std::vector<uint64_t> g(K), e(E);
    std::vector<uint16_t> c(o.S * K);
    std::vector<uint8_t> nv(E);
    std::vector<uint32_t> at(K + 1, 0);
    for (uint64_t i = 0; i < E; i++) at[id[Graph::from(o.edges[i])] + 1]++;
    for (uint64_t v = 0; v < K; v++) {
        if (v + PF < K) __builtin_prefetch(&o.hashes[ord[v + PF]]);
        g[v] = o.hashes[ord[v]];
        at[v + 1] += at[v];
    }
    for (uint64_t s = 0; s < o.S; s++)
        for (uint64_t v = 0; v < K; v++) {
            if (v + PF < K) __builtin_prefetch(&o.counts[s * K + ord[v + PF]]);
            c[s * K + v] = o.counts[s * K + ord[v]];
        }
    for (uint64_t i = 0; i < E; i++) {
        if (i + PF < E) { __builtin_prefetch(&id[Graph::to(o.edges[i + PF])]); __builtin_prefetch(&at[id[Graph::from(o.edges[i + PF])]]); }
        uint64_t u = id[Graph::from(o.edges[i])], p = at[u]++;
        e[p] = u << 32 | id[Graph::to(o.edges[i])];
        nv[p] = o.novel[i];
    }
    for (uint64_t p = 0, q; p < E; p = q) {
        for (q = p + 1; q < E && e[q] >> 32 == e[p] >> 32; q++);
        for (uint64_t i = p + 1; i < q; i++)
            for (uint64_t j = i; j > p && e[j] < e[j - 1]; j--) { std::swap(e[j], e[j - 1]); std::swap(nv[j], nv[j - 1]); }
    }
    o.hashes.swap(g);
    o.counts.swap(c);
    o.edges.swap(e);
    o.novel.swap(nv);
}

int cmd_graph(const Args& a) {
    return MakeGraph(a).run();
}
