#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "core.h"
#include "extract.h"
#include "io.h"
#include "lps.h"

static const char SYM[5] = "ACGT";

struct scr { rdcore r; lcp_soa o[2]; uint32_t *lab, *lsz, *map; char* hs; bits bw; uint32_t cap; int cur; };

struct ectx { scr* s; uint64_t* nc; int lev, w; };

static void sgrow(scr* s, uint32_t n) {
    if (n <= s->cap) return;
    s->cap = 2 * n;
    for (int i = 0; i < 2; i++) {
        lcp_soa* o = &s->o[i];
        o->rep = (uint64_t*)std::realloc(o->rep, s->cap * 8);
        o->st = (uint64_t*)std::realloc(o->st, s->cap * 8);
        o->en = (uint64_t*)std::realloc(o->en, s->cap * 8);
        o->rep32 = (uint32_t*)std::realloc(o->rep32, s->cap * 4);
        o->sz = (uint32_t*)std::realloc(o->sz, s->cap * 4);
    }
    s->lab = (uint32_t*)std::realloc(s->lab, s->cap * 4);
    s->lsz = (uint32_t*)std::realloc(s->lsz, s->cap * 4);
    s->map = (uint32_t*)std::realloc(s->map, s->cap * 4);
    s->hs = (char*)std::realloc(s->hs, s->cap);
}

static int parse_lev(scr* s, int k, int lev) {
    int n = parse1_soa(s->hs, s->hs + k, &s->o[0], 0), c = 0;
    for (int v = 1; v < lev && n > 0; v++, c ^= 1) {
        lcp_dct_soa(&s->o[c], n, v, s->lab, s->lsz);
        n = parse3_soa(s->lab, s->lsz, &s->o[c], n, &s->o[c ^ 1]);
    }
    s->cur = c;
    return n > 0 ? n : 0;
}

static void emit(block* bl, scr* s, int n, int w) {
    uint32_t prev = 0, m = (1u << w) - 1;
    for (int i = 0; i < n; i++) {
        uint32_t a = s->map[s->o[s->cur].st[i]], len = s->map[s->o[s->cur].en[i]] - a, dst = a - prev;
        if (dst >= m || len > m) { wbits(bl, &s->bw, m, w); wbits(bl, &s->bw, dst, 32); wbits(bl, &s->bw, len, 32); }
        else { wbits(bl, &s->bw, dst, w); wbits(bl, &s->bw, len, w); }
        prev = a;
    }
}

static void rc(scr* s, uint32_t k, uint32_t len) {
    for (uint32_t i = 0, j = k - 1; i <= j; i++, j--) {
        char a = s->hs[i], b = s->hs[j];
        s->hs[i] = SYM[CODE.t[(uint8_t)b] ^ 3];
        s->hs[j] = SYM[CODE.t[(uint8_t)a] ^ 3];
    }
    for (uint32_t i = 0, j = k; i <= j; i++, j--) {
        uint32_t a = s->map[i], b = s->map[j];
        s->map[i] = len - b;
        s->map[j] = len - a;
    }
}

static void ex_work(block* bl) {
    ectx* c = (ectx*)bl->ctx;
    scr* s = &c->s[bl->id];
    rdcore* r = &s->r;
    for (size_t off = 0; off < bl->inlen; off += svb_size(bl->in + off, bl->inlen - off)) {
        collapse(r, (const uint8_t*)bl->in + off);
        sgrow(s, r->k + 2);
        uint32_t q = 0;
        for (uint32_t i = 0; i < r->k; i++) s->hs[i] = SYM[r->hs[i >> 5] >> (62 - 2 * (i & 31)) & 3];
        for (uint32_t j = 0; j < (r->len + 31) / 32; j++)
            for (uint64_t m = r->mk[j]; m; m ^= 1ull << (63 - __builtin_clzll(m))) s->map[q++] = 32 * j + __builtin_clzll(m) / 2;
        s->map[r->k] = r->len;
        size_t t0 = bl->outlen;
        char tb[6] = {0};
        bput(bl, tb, 6);
        int nf = 0, nr = 0;
        if (r->k >= 3) {
            nf = parse_lev(s, r->k, c->lev);
            emit(bl, s, nf, c->w);
            rc(s, r->k, r->len);
            nr = parse_lev(s, r->k, c->lev);
            emit(bl, s, nr, c->w);
        }
        wflush(bl, &s->bw);
        uint8_t* u = (uint8_t*)bl->out + t0;
        size_t ub = bl->outlen - t0;
        u[0] = (nf + nr) >> 8; u[1] = nf + nr; u[2] = nf >> 8; u[3] = nf; u[4] = ub >> 8; u[5] = ub;
        c->nc[bl->id] += nf + nr;
    }
}

int cmd_extract(const char* in, const char* out, int t, int l) {
    std::string path;
    if (!out) { path = std::string(in) + ".exc"; out = path.c_str(); }
    LCP_INIT();
    FILE* fi = std::fopen(in, "rb");
    bin_header h;
    std::fread(&h, 1, 32, fi);
    ectx c = {(scr*)std::calloc(t, sizeof(scr)), (uint64_t*)std::calloc(t, 8), l, l + 4};
    FILE* fo = open_tmp(out);
    uint64_t n = stream_run(fi, nullptr, fo, 32, t, svb_size, nullptr, ex_work, &c), cores = 0;
    std::fclose(fi);
    for (int i = 0; i < t; i++) {
        scr* s = &c.s[i];
        cores += c.nc[i];
        free_core(&s->r);
        for (int j = 0; j < 2; j++) { std::free(s->o[j].rep); std::free(s->o[j].rep32); std::free(s->o[j].st); std::free(s->o[j].en); std::free(s->o[j].sz); }
        std::free(s->lab); std::free(s->lsz); std::free(s->map); std::free(s->hs);
    }
    std::free(c.s);
    std::free(c.nc);
    bin_header o = {{'S', 'V', 'B', 'X'}, 3, n, cores, {(uint8_t)l, (uint8_t)(l + 4)}};
    std::fseek(fo, 0, SEEK_SET);
    std::fwrite(&o, 1, 32, fo);
    commit_file(fo, out);
    return 0;
}
