#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include "core.h"
#include "extract.h"
#include "io.h"
#include "lps.h"

static const char SYM[5] = "ACGT";

struct scr { rdcore r; lcp_soa o[2]; bits bw; uint64_t nc, *mem; uint32_t *lab, *lsz; std::vector<uint32_t> map; std::vector<char> hs; };

struct ectx { scr* s; int lev, w; };

static void sgrow(scr* s, uint32_t k) {
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
    for (uint32_t i = off; i < off + k; i++) out[i - off] = SYM[w[i >> 5] >> (62 - 2 * (i & 31)) & 3];
}

static int parse(scr* s, uint32_t k, int lev, int w) {
    int n = parse1_soa(s->hs.data(), s->hs.data() + k, &s->o[0], 0), c = 0;
    for (int v = 1; v < lev && n > 0; v++, c ^= 1) {
        lcp_dct_soa(&s->o[c], n, v, s->lab, s->lsz);
        n = parse3_soa(s->lab, s->lsz, &s->o[c], n, &s->o[c ^ 1]);
    }
    uint32_t prev = 0, m = (1u << w) - 1;
    for (int i = 0; i < n; i++) {
        uint32_t a = s->map[s->o[c].st[i]], len = s->map[s->o[c].en[i]] - a, d = a - prev;
        bool e = d >= m || len > m;
        if (e) wbits(&s->bw, m, w);
        wbits(&s->bw, d, e ? 32 : w);
        wbits(&s->bw, len, e ? 32 : w);
        prev = a;
    }
    return n > 0 ? n : 0;
}

static void ex_work(block* bl) {
    ectx* c = (ectx*)bl->ctx;
    scr* s = &c->s[bl->id];
    rdcore* r = &s->r;
    s->bw.out = &bl->out;
    for (size_t off = 0; off < bl->inlen[0]; off += svb_size(bl->in[0] + off, bl->inlen[0] - off)) {
        collapse(r, (const uint8_t*)bl->in[0] + off);
        sgrow(s, r->k);
        uint32_t q = 0, nf = 0, nr = 0;
        for (uint32_t j = 0; j < (r->len + 31) / 32; j++)
            for (uint64_t m = r->mk[j]; m; m ^= 1ull << (63 - __builtin_clzll(m))) s->map[q++] = 32 * j + __builtin_clzll(m) / 2;
        s->map[r->k] = r->len;
        size_t t0 = bl->out.size();
        bl->out.resize(t0 + 6);
        if (r->k >= 3) {
            ascii(r->hs.data(), 0, r->k, s->hs.data());
            nf = parse(s, r->k, c->lev, c->w);
            ascii(r->rs.data(), r->pad, r->k, s->hs.data());
            std::reverse(s->map.begin(), s->map.begin() + r->k + 1);
            for (uint32_t i = 0; i <= r->k; i++) s->map[i] = r->len - s->map[i];
            nr = parse(s, r->k, c->lev, c->w);
        }
        wflush(&s->bw);
        uint8_t* u = (uint8_t*)bl->out.data() + t0;
        size_t ub = bl->out.size() - t0;
        u[0] = (nf + nr) >> 8; u[1] = nf + nr; u[2] = nf >> 8; u[3] = nf; u[4] = ub >> 8; u[5] = ub;
        s->nc += nf + nr;
    }
}

int cmd_extract(const char* in, const char* out, int t, int l) {
    LCP_INIT();
    bin_header h;
    std::vector<scr> s(t);
    ectx c = {s.data(), l, l + 4};
    FILE* fo = open_tmp(out);
    bin_header o = {{'S', 'V', 'B', 'X'}, 3, stream_run(open_in(in, &h), nullptr, fo, t, svb_size, nullptr, ex_work, &c), 0, {(uint8_t)l, (uint8_t)(l + 4)}};
    for (scr& x : s) { o.total += x.nc; std::free(x.mem); }
    finish(fo, out, &o);
    return 0;
}
