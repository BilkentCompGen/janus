#include <algorithm>
#include <string>
#include <vector>
#include "graph.h"
#include "io.h"

static void unite(std::vector<uint64_t>& a, const std::vector<uint64_t>& b, std::vector<uint64_t>& t) {
    t.resize(a.size() + b.size());
    t.resize(std::set_union(a.begin(), a.end(), b.begin(), b.end(), t.begin()) - t.begin());
    a.swap(t);
}

static void hashes(const char* path, std::vector<uint64_t>& hs) {
    bin_header h;
    FILE* f = open_in(path, &h);
    hs.resize(h.nreads);
    std::fread(hs.data(), 8, h.nreads, f);
    std::fclose(f);
}

static void sample(const char* path, const std::vector<uint64_t>& g, uint16_t* c, std::vector<uint64_t>& e, std::vector<uint16_t>& ec) {
    bin_header h;
    FILE* f = open_in(path, &h);
    uint64_t K = h.nreads, E = h.total, at = sizeof(bin_header) + K * 10 + E * 2, pos = 0, q = K * h.resv[1];
    std::vector<uint64_t> hs(K);
    std::vector<uint16_t> n(K);
    std::vector<uint32_t> id(K);
    ec.resize(E);
    std::fread(hs.data(), 8, K, f);
    std::fread(n.data(), 2, K, f);
    std::fread(ec.data(), 2, E, f);
    std::fseek(f, 0, SEEK_END);
    std::vector<uint8_t> bb(std::ftell(f) - at + 8);
    std::fseek(f, at, SEEK_SET);
    std::fread(bb.data(), 1, bb.size() - 8, f);
    std::fclose(f);
    for (uint64_t j = 0, k = 0; j < K; j++) {
        while (g[k] < hs[j]) k++;
        id[j] = k;
        c[k] = n[j];
    }
    e.clear();
    for (uint64_t i = 0; i < K; i++)
        for (uint32_t d = rbits(bb.data(), &pos, h.resv[1]); d; d--)
            e.push_back((uint64_t)id[i] << 32 | id[rbits(bb.data(), &q, h.resv[0])]);
}

int cmd_graph(const char* in, const char* out) {
    std::vector<std::string> ps = split(in);
    std::vector<uint64_t> g, e, s, t;
    std::vector<uint16_t> sc;
    for (auto& x : ps) {
        hashes(x.c_str(), s);
        unite(g, s, t);
    }
    uint64_t K = g.size(), S = ps.size();
    std::vector<uint16_t> c(S * K);
    for (uint64_t i = 0; i < S; i++) {
        sample(ps[i].c_str(), g, c.data() + i * K, s, sc);
        unite(e, s, t);
    }
    t = {};
    uint64_t E = e.size();
    std::vector<uint16_t> ec(S * E);
    for (uint64_t i = 0; i < S; i++) {
        sample(ps[i].c_str(), g, c.data() + i * K, s, sc);
        for (uint64_t j = 0, k = 0; j < s.size(); j++) {
            while (e[k] < s[j]) k++;
            ec[i * E + k] = sc[j];
        }
    }
    FILE* fo = open_tmp(out);
    std::fseek(fo, sizeof(bin_header), SEEK_SET);
    std::fwrite(g.data(), 8, K, fo);
    std::fwrite(c.data(), 2, S * K, fo);
    std::fwrite(e.data(), 8, E, fo);
    std::fwrite(ec.data(), 2, S * E, fo);
    bin_header o = {{'S', 'V', 'B', 'G'}, 1, K, E, {(uint8_t)S}};
    finish(fo, out, &o);
    return 0;
}
