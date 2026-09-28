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
    *x = {h.nreads, h.total, h.resv[0], std::vector<uint64_t>(h.nreads), std::vector<uint64_t>(h.total),
          std::vector<uint16_t>(h.resv[0] * h.nreads), std::vector<uint16_t>(h.resv[0] * h.total), std::vector<uint8_t>(h.resv[0] * h.nreads)};
    std::fread(x->g.data(), 8, x->K, f);
    std::fread(x->c.data(), 2, x->S * x->K, f);
    std::fread(x->e.data(), 8, x->E, f);
    std::fread(x->ec.data(), 2, x->S * x->E, f);
    std::fread(x->sh.data(), 1, x->S * x->K, f);
    std::fclose(f);
}

void write_gph(const char* path, const gph* x) {
    FILE* fo = open_tmp(path);
    std::fseek(fo, sizeof(bin_header), SEEK_SET);
    std::fwrite(x->g.data(), 8, x->K, fo);
    std::fwrite(x->c.data(), 2, x->S * x->K, fo);
    std::fwrite(x->e.data(), 8, x->E, fo);
    std::fwrite(x->ec.data(), 2, x->S * x->E, fo);
    std::fwrite(x->sh.data(), 1, x->S * x->K, fo);
    bin_header o = {{'S', 'V', 'B', 'G'}, 1, x->K, x->E, {(uint8_t)x->S}};
    finish(fo, path, &o);
}

void make_adj(const gph* x, adj* a) {
    a->oo.assign(x->K + 1, 0);
    a->io.assign(x->K + 1, 0);
    a->on.resize(x->E);
    a->in.resize(x->E);
    for (uint64_t v : x->e) { a->oo[(v >> 32) + 1]++; a->io[(v & 0xffffffff) + 1]++; }
    for (uint64_t i = 0; i < x->K; i++) { a->oo[i + 1] += a->oo[i]; a->io[i + 1] += a->io[i]; }
    std::vector<uint32_t> p(a->io.begin(), a->io.end() - 1);
    for (uint64_t i = 0; i < x->E; i++) { a->on[i] = x->e[i] & 0xffffffff; a->in[p[x->e[i] & 0xffffffff]++] = x->e[i] >> 32; }
}

void masks(const gph* x, std::vector<uint64_t>& m) {
    m.assign(x->K, 0);
    for (uint64_t i = 0; i < x->S * x->K; i++) m[i % x->K] |= (uint64_t)(x->c[i] > 0) << (i / x->K);
}

void make_look(const std::vector<uint64_t>& g, std::vector<uint32_t>& off) {
    off.assign((1 << 24) + 1, 0);
    for (uint64_t h : g) off[(h >> 40) + 1]++;
    for (uint32_t i = 0; i < 1 << 24; i++) off[i + 1] += off[i];
}

uint32_t look(const std::vector<uint64_t>& g, const std::vector<uint32_t>& off, uint64_t h) {
    uint32_t i = off[h >> 40], e = off[(h >> 40) + 1];
    while (i < e && g[i] < h) i++;
    return i;
}

std::vector<uint64_t> read_set(const char* path) {
    bin_header h;
    FILE* f = open_in(path, &h);
    std::vector<uint64_t> v(h.nreads);
    std::fread(v.data(), 8, h.nreads, f);
    std::fclose(f);
    return v;
}

void write_set(const char* path, const char* magic, const std::vector<uint64_t>& v, uint64_t aux) {
    FILE* fo = open_tmp(path);
    std::fseek(fo, sizeof(bin_header), SEEK_SET);
    std::fwrite(v.data(), 8, v.size(), fo);
    bin_header o = {{magic[0], magic[1], magic[2], magic[3]}, 1, v.size(), aux, {}};
    finish(fo, path, &o);
}
