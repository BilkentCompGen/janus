#include <cstdlib>
#include <cstring>
#include "files.h"
#include "io.h"

void read_cnt(const char* path, cnt* x) {
    bin_header h;
    FILE* f = open_in(path, &h);
    uint64_t K = h.nreads, E = h.total, pos = 0, q = K * h.resv[1];
    x->hs.resize(K);
    x->n.resize(K);
    x->ec.resize(E);
    x->e.clear();
    std::fread(x->hs.data(), 8, K, f);
    std::fread(x->n.data(), 2, K, f);
    std::fread(x->ec.data(), 2, E, f);
    long at = std::ftell(f);
    std::fseek(f, 0, SEEK_END);
    std::vector<uint8_t> bb(std::ftell(f) - at + 8);
    std::fseek(f, at, SEEK_SET);
    std::fread(bb.data(), 1, bb.size() - 8, f);
    std::fclose(f);
    for (uint64_t i = 0; i < K; i++)
        for (uint32_t d = rbits(bb.data(), &pos, h.resv[1]); d; d--) x->e.push_back(i << 32 | rbits(bb.data(), &q, h.resv[0]));
}

void read_gph(const char* path, gph* x) {
    bin_header h;
    FILE* f = open_in(path, &h);
    uint64_t U, X = h.total, j = 0;
    std::fread(&U, 8, 1, f);
    *x = {h.nreads, 0, h.resv[0], std::vector<uint64_t>(h.nreads), {}, std::vector<uint16_t>(h.resv[0] * h.nreads), {}, std::vector<uint32_t>(U + 1)};
    std::vector<uint64_t> ex(X);
    std::vector<uint8_t> nx(X), ni(x->K);
    std::fread(x->g.data(), 8, x->K, f);
    std::fread(x->c.data(), 2, x->S * x->K, f);
    std::fread(x->st.data(), 4, U + 1, f);
    std::fread(ex.data(), 8, X, f);
    std::fread(nx.data(), 1, X, f);
    std::fread(ni.data(), 1, x->K, f);
    std::fclose(f);
    x->e.reserve(X + x->K - U);
    x->nv.reserve(X + x->K - U);
    for (uint64_t u = 0; u < U; u++)
        for (uint64_t v = x->st[u]; v < x->st[u + 1]; v++)
            if (v + 1 < x->st[u + 1]) { x->e.push_back(v << 32 | (v + 1)); x->nv.push_back(ni[v]); }
            else for (; j < X && ex[j] >> 32 == v; j++) { x->e.push_back(ex[j]); x->nv.push_back(nx[j]); }
    x->E = x->e.size();
}

void write_gph(const char* path, const gph* x) {
    uint64_t U = x->st.size() - 1;
    std::vector<uint64_t> ex;
    std::vector<uint8_t> nx, ni(x->K, 0), in(x->K, 0);
    for (uint64_t u = 0; u < U; u++) for (uint64_t v = x->st[u]; v + 1 < x->st[u + 1]; v++) in[v] = 1;
    for (uint64_t i = 0; i < x->E; i++)
        if (in[x->e[i] >> 32]) ni[x->e[i] >> 32] = x->nv[i];
        else { ex.push_back(x->e[i]); nx.push_back(x->nv[i]); }
    FILE* fo = open_tmp(path);
    std::fwrite(&U, 8, 1, fo);
    std::fwrite(x->g.data(), 8, x->K, fo);
    std::fwrite(x->c.data(), 2, x->S * x->K, fo);
    std::fwrite(x->st.data(), 4, U + 1, fo);
    std::fwrite(ex.data(), 8, ex.size(), fo);
    std::fwrite(nx.data(), 1, nx.size(), fo);
    std::fwrite(ni.data(), 1, x->K, fo);
    bin_header o = header('G', 2, x->K, ex.size());
    o.resv[0] = x->S;
    finish(fo, path, &o);
}

void make_adj(const gph* x, adj* a) {
    const uint64_t M = 0xffffffff;
    a->oo.assign(x->K + 1, 0);
    a->io.assign(x->K + 1, 0);
    a->on.resize(x->E);
    a->in.resize(x->E);
    for (uint64_t i = 0; i < x->E; i++) {
        if (i + PF < x->E) __builtin_prefetch(&a->io[(x->e[i + PF] & M) + 1]);
        a->oo[(x->e[i] >> 32) + 1]++;
        a->io[(x->e[i] & M) + 1]++;
    }
    for (uint64_t i = 0; i < x->K; i++) { a->oo[i + 1] += a->oo[i]; a->io[i + 1] += a->io[i]; }
    std::vector<uint32_t> p(a->io.begin(), a->io.end() - 1);
    for (uint64_t i = 0; i < x->E; i++) {
        if (i + PF < x->E) { __builtin_prefetch(&p[x->e[i + PF] & M]); __builtin_prefetch(&a->in[p[x->e[i + PF / 2] & M]]); }
        a->on[i] = x->e[i] & M;
        a->in[p[x->e[i] & M]++] = x->e[i] >> 32;
    }
}

void read_fasta(const char* path, seqs* f) {
    FILE* in = std::fopen(path, "r");
    char* ln = nullptr;
    size_t cap = 0;
    for (ssize_t n; (n = getline(&ln, &cap, in)) > 0;) {
        while (n && (ln[n - 1] == '\n' || ln[n - 1] == '\r')) n--;
        if (ln[0] == '>') { f->name.push_back(std::string(ln + 1, std::strcspn(ln + 1, " \t\n"))); f->seq.push_back(""); }
        else {
            for (ssize_t i = 0; i < n; i++) ln[i] &= 0xdf;
            f->seq.back().append(ln, n);
        }
    }
    std::free(ln);
    std::fclose(in);
}
