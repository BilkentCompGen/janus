#include <cstring>
#include "files.h"
#include "io.h"

void Graph::read(const char* path) {
    Input in(path);
    uint64_t U, X = in.h.total, j = 0;
    std::fread(&U, 8, 1, in.f);
    K = in.h.n;
    S = in.h.resv[0];
    hashes.assign(K, 0);
    counts.assign(S * K, 0);
    starts.assign(U + 1, 0);
    std::vector<uint64_t> ex(X);
    std::vector<uint8_t> nx(X), ni(K);
    std::fread(hashes.data(), 8, K, in.f);
    std::fread(counts.data(), 2, S * K, in.f);
    std::fread(starts.data(), 4, U + 1, in.f);
    std::fread(ex.data(), 8, X, in.f);
    std::fread(nx.data(), 1, X, in.f);
    std::fread(ni.data(), 1, K, in.f);
    std::fclose(in.f);
    edges.clear();
    novel.clear();
    edges.reserve(X + K - U);
    novel.reserve(X + K - U);
    for (uint64_t u = 0; u < U; u++)
        for (uint64_t v = starts[u]; v < starts[u + 1]; v++)
            if (v + 1 < starts[u + 1]) { edges.push_back(v << 32 | (v + 1)); novel.push_back(ni[v]); }
            else for (; j < X && from(ex[j]) == v; j++) { edges.push_back(ex[j]); novel.push_back(nx[j]); }
    E = edges.size();
}

void Graph::write(const char* path) const {
    uint64_t U = starts.size() - 1;
    std::vector<uint64_t> ex;
    std::vector<uint8_t> nx, ni(K, 0), inner(K, 0);
    for (uint64_t u = 0; u < U; u++)
        for (uint64_t v = starts[u]; v + 1 < starts[u + 1]; v++) inner[v] = 1;
    for (uint64_t i = 0; i < E; i++)
        if (inner[from(edges[i])]) ni[from(edges[i])] = novel[i];
        else { ex.push_back(edges[i]); nx.push_back(novel[i]); }
    Output out(path);
    std::fwrite(&U, 8, 1, out.f);
    std::fwrite(hashes.data(), 8, K, out.f);
    std::fwrite(counts.data(), 2, S * K, out.f);
    std::fwrite(starts.data(), 4, U + 1, out.f);
    std::fwrite(ex.data(), 8, ex.size(), out.f);
    std::fwrite(nx.data(), 1, nx.size(), out.f);
    std::fwrite(ni.data(), 1, K, out.f);
    Header h('G', 2, K, ex.size());
    h.resv[0] = S;
    out.finish(h);
}

Adjacency::Adjacency(const Graph& g) : first_out(g.K + 1, 0), target(g.E), first_in(g.K + 1, 0), source(g.E) {
    const std::vector<uint64_t>& e = g.edges;
    for (uint64_t i = 0; i < g.E; i++) {
        if (i + PF < g.E) __builtin_prefetch(&first_in[Graph::to(e[i + PF]) + 1]);
        first_out[Graph::from(e[i]) + 1]++;
        first_in[Graph::to(e[i]) + 1]++;
    }
    for (uint64_t i = 0; i < g.K; i++) { first_out[i + 1] += first_out[i]; first_in[i + 1] += first_in[i]; }
    std::vector<uint32_t> p(first_in.begin(), first_in.end() - 1);
    for (uint64_t i = 0; i < g.E; i++) {
        if (i + PF < g.E) { __builtin_prefetch(&p[Graph::to(e[i + PF])]); __builtin_prefetch(&source[p[Graph::to(e[i + PF / 2])]]); }
        target[i] = Graph::to(e[i]);
        source[p[Graph::to(e[i])]++] = Graph::from(e[i]);
    }
}

Fasta::Fasta(const char* path) {
    Lines in(path);
    while (in.next()) {
        ssize_t n = in.n;
        char* ln = in.s;
        while (n && (ln[n - 1] == '\n' || ln[n - 1] == '\r')) n--;
        if (ln[0] == '>') { name.push_back(std::string(ln + 1, std::strcspn(ln + 1, " \t\n"))); seq.push_back(""); }
        else {
            for (ssize_t i = 0; i < n; i++) ln[i] &= 0xdf;
            seq.back().append(ln, n);
        }
    }
}

Fai::Fai(const char* fasta) {
    char nm[256];
    unsigned long f[4];
    FILE* in = std::fopen((std::string(fasta) + ".fai").c_str(), "r");
    while (std::fscanf(in, "%255s %lu %lu %lu %lu", nm, &f[0], &f[1], &f[2], &f[3]) == 5) m[nm] = {f[0], f[1], f[2], f[3]};
    std::fclose(in);
}

std::string Fai::range(FILE* fa, const std::string& c, uint32_t a, uint32_t b) {
    const FaiEntry& v = m[c];
    uint64_t s = v.off + a / v.lb * v.lw + a % v.lb, e = v.off + b / v.lb * v.lw + b % v.lb;
    std::string r(e - s, 0), o;
    std::fseek(fa, s, SEEK_SET);
    r.resize(std::fread(&r[0], 1, e - s, fa));
    for (char ch : r)
        if (ch != '\n') o += ch & 0xdf;
    return o;
}
