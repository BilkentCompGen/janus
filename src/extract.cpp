#include <algorithm>
#include <cstdlib>
#include "cmd.h"
#include "extract.h"
#include "io.h"

static const int W = LEV + 4; // bits of a core's start delta and length in .exc; all ones escapes to two 32-bit fields

struct ectx { std::vector<scr> s; std::vector<uint64_t> nc; };

static void sgrow(scr* s, rdcore* r, uint32_t k) {
    if (r->a.size() < 2 * k + 2) { r->a.resize(4 * k + 4); r->b.resize(4 * k + 4); }
    if (s->map.size() > k) return;
    k *= 2;
    uint64_t c = (uint64_t)(k / LCP_CONSTANT_FACTOR) + 8, *p = s->mem = (uint64_t*)std::realloc(s->mem, c * 72);
    for (int i = 0; i < 2; i++, p += 4 * c) s->o[i] = {p, (uint32_t*)(p + 3 * c), p + c, p + 2 * c, (uint32_t*)(p + 3 * c) + c};
    s->lab = (uint32_t*)p;
    s->lsz = s->lab + c;
    s->map.resize(k + 1);
    s->hs.resize(k);
}

static void ascii(const uint64_t* w, uint32_t off, uint32_t k, char* out) {
    for (uint32_t i = off; i < off + k; i++) out[i - off] = "ACGT"[w[i >> 5] >> (62 - 2 * (i & 31)) & 3];
}

static uint32_t parse(scr* s, rdcore* r, uint32_t k, int lev, uint32_t at) {
    int n = parse1_soa(s->hs.data(), s->hs.data() + k, &s->o[0], 0), c = 0;
    for (int v = 1; v < lev && n > 0; v++, c ^= 1) {
        lcp_dct_soa(&s->o[c], n, v, s->lab, s->lsz);
        n = parse3_soa(s->lab, s->lsz, &s->o[c], n, &s->o[c ^ 1]);
    }
    for (int i = 0; i < n; i++) { r->a[at + i] = s->map[s->o[c].st[i]]; r->b[at + i] = s->map[s->o[c].en[i]]; }
    return n > 0 ? n : 0;
}

void cores(scr* s, rdcore* r, const uint8_t* rec, int lev) {
    collapse(r, rec);
    sgrow(s, r, r->k);
    uint32_t q = 0;
    for (uint32_t j = 0; j < (r->len + 31) / 32; j++)
        for (uint64_t m = r->mk[j]; m; m ^= 1ull << (63 - __builtin_clzll(m))) s->map[q++] = 32 * j + __builtin_clzll(m) / 2;
    s->map[r->k] = r->len;
    r->n = r->nf = 0;
    if (r->k < 3) return;
    ascii(r->hs.data(), 0, r->k, s->hs.data());
    r->nf = parse(s, r, r->k, lev, 0);
    ascii(r->rs.data(), r->pad, r->k, s->hs.data());
    std::reverse(s->map.begin(), s->map.begin() + r->k + 1);
    for (uint32_t i = 0; i <= r->k; i++) s->map[i] = r->len - s->map[i];
    r->n = r->nf + parse(s, r, r->k, lev, r->nf);
}

static void ex_work(block* bl) {
    ectx* c = (ectx*)bl->ctx;
    scr* s = &c->s[bl->id];
    const rdcore* r = &s->r;
    bits bw = {&bl->out, 0, 0};
    uint32_t m = (1u << W) - 1;
    for (size_t off = 0; off < bl->inlen[0]; off += svb_size(bl->in[0] + off, bl->inlen[0] - off)) {
        cores(s, &s->r, (const uint8_t*)bl->in[0] + off, LEV);
        size_t t0 = bl->out.size();
        bl->out.resize(t0 + 6);
        for (uint32_t i = 0, prev = 0; i < r->n; prev = r->a[i++]) {
            if (i == r->nf) prev = 0;
            uint32_t d = r->a[i] - prev, len = r->b[i] - r->a[i];
            bool e = d >= m || len > m;
            if (e) wbits(&bw, m, W);
            wbits(&bw, d, e ? 32 : W);
            wbits(&bw, len, e ? 32 : W);
        }
        wflush(&bw);
        uint8_t* u = (uint8_t*)bl->out.data() + t0;
        size_t ub = bl->out.size() - t0;
        u[0] = r->n >> 8; u[1] = r->n; u[2] = r->nf >> 8; u[3] = r->nf; u[4] = ub >> 8; u[5] = ub;
        c->nc[bl->id] += r->n;
    }
}

int cmd_extract(const char* in, const char* out, int t) {
    LCP_INIT();
    bin_header h;
    ectx c = {std::vector<scr>(t), std::vector<uint64_t>(t, 0)};
    FILE* fo = open_tmp(out);
    bin_header o = header('X', 3, stream_run(open_in(in, &h), nullptr, fo, t, ex_work, &c), 0);
    o.resv[0] = LEV; o.resv[1] = W;
    for (int i = 0; i < t; i++) { o.total += c.nc[i]; std::free(c.s[i].mem); }
    finish(fo, out, &o);
    return 0;
}
