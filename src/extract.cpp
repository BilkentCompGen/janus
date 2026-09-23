#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "extract.h"
#include "io.h"
#include "lps.h"

static const char SYM[5] = "ACGT";

struct scr {
    uint64_t* rep[2]; uint32_t* rep32[2]; uint64_t* st[2]; uint64_t* en[2]; uint32_t* sz[2];
    uint32_t *lab, *lsz, *map; char* hs;
    bits bw; int cap, cur;
};

struct ectx { scr* s; uint64_t* nc; int lev, w; };

static void sgrow(scr* s, int n) {
    if (n <= s->cap) return;
    s->cap = n;
    for (int i = 0; i < 2; i++) {
        s->rep[i] = (uint64_t*)std::realloc(s->rep[i], (size_t)n * 8);
        s->st[i] = (uint64_t*)std::realloc(s->st[i], (size_t)n * 8);
        s->en[i] = (uint64_t*)std::realloc(s->en[i], (size_t)n * 8);
        s->rep32[i] = (uint32_t*)std::realloc(s->rep32[i], (size_t)n * 4);
        s->sz[i] = (uint32_t*)std::realloc(s->sz[i], (size_t)n * 4);
    }
    s->lab = (uint32_t*)std::realloc(s->lab, (size_t)n * 4);
    s->lsz = (uint32_t*)std::realloc(s->lsz, (size_t)n * 4);
    s->map = (uint32_t*)std::realloc(s->map, (size_t)n * 4);
    s->hs = (char*)std::realloc(s->hs, (size_t)n);
}

static void wrec(block* bl, scr* s, uint32_t dst, uint32_t len, int w) {
    uint32_t m = (1u << w) - 1;
    if (dst >= m || len > m) { wbits(bl, &s->bw, m, w); wbits(bl, &s->bw, dst, 32); wbits(bl, &s->bw, len, 32); }
    else { wbits(bl, &s->bw, dst, w); wbits(bl, &s->bw, len, w); }
}

static int parse_lev(scr* s, int k, int lev) {
    lcp_soa o = {s->rep[0], s->rep32[0], s->st[0], s->en[0], s->sz[0]};
    int n = parse1_soa(s->hs, s->hs + k, &o, 0), c = 0;
    for (int v = 1; v < lev && n > 0; v++) {
        lcp_soa a = {s->rep[c], s->rep32[c], s->st[c], s->en[c], s->sz[c]};
        lcp_soa d = {s->rep[c ^ 1], s->rep32[c ^ 1], s->st[c ^ 1], s->en[c ^ 1], s->sz[c ^ 1]};
        lcp_dct_soa(&a, n, v, s->lab, s->lsz);
        n = parse3_soa(s->lab, s->lsz, &a, n, &d);
        c ^= 1;
    }
    s->cur = c;
    return n > 0 ? n : 0;
}

static void emit(block* bl, scr* s, int n, int w) {
    const uint64_t* st = s->st[s->cur];
    const uint64_t* en = s->en[s->cur];
    uint32_t prev = 0;
    for (int i = 0; i < n; i++) {
        uint32_t a = s->map[st[i]], b = s->map[en[i]];
        wrec(bl, s, a - prev, b - a, w);
        prev = a;
    }
}

static void rc(scr* s, int k, uint32_t len) {
    for (int i = 0, j = k - 1; i <= j; i++, j--) {
        char a = s->hs[i], b = s->hs[j];
        s->hs[i] = SYM[CODE.t[(uint8_t)b] ^ 3];
        s->hs[j] = SYM[CODE.t[(uint8_t)a] ^ 3];
    }
    for (int i = 0, j = k; i <= j; i++, j--) {
        uint32_t a = s->map[i], b = s->map[j];
        s->map[i] = len - b;
        s->map[j] = len - a;
    }
}

static void ex_work(block* bl) {
    ectx* c = (ectx*)bl->ctx;
    scr* s = &c->s[bl->id];
    size_t off = 0, off2 = 0;
    while (off < bl->inlen) {
        const uint8_t* p = (const uint8_t*)bl->in + off;
        const uint8_t* h = (const uint8_t*)bl->in2 + off2;
        uint32_t tag, k;
        std::memcpy(&tag, p, 4);
        std::memcpy(&k, h, 4);
        uint32_t len = tag & 0x7fffffff;
        sgrow(s, (int)k + 2);
        for (uint32_t i = 0; i < k; i++) s->hs[i] = SYM[(h[4 + (i >> 2)] >> (6 - 2 * (i & 3))) & 3];
        uint32_t q = 0;
        int prev = -1;
        for (uint32_t i = 0; i < len; i++) {
            int v = tag & 0x80000000 ? CODE.t[(uint8_t)(p[4 + i] & 0xdf)] : (p[4 + (i >> 2)] >> (6 - 2 * (i & 3))) & 3;
            if (v != prev) s->map[q++] = i;
            prev = v;
        }
        s->map[k] = len;
        size_t t0 = bl->outlen;
        char tb[6] = {0};
        bput(bl, tb, 6);
        int nf = 0, nr = 0;
        if (k >= 3) {
            nf = parse_lev(s, (int)k, c->lev);
            emit(bl, s, nf, c->w);
            rc(s, (int)k, len);
            nr = parse_lev(s, (int)k, c->lev);
            emit(bl, s, nr, c->w);
        }
        wflush(bl, &s->bw);
        uint8_t* u = (uint8_t*)bl->out + t0;
        size_t ub = bl->outlen - t0;
        u[0] = (nf + nr) >> 8; u[1] = nf + nr; u[2] = nf >> 8; u[3] = nf; u[4] = ub >> 8; u[5] = ub;
        c->nc[bl->id] += (uint64_t)(nf + nr);
        off += svb_size((const char*)p, bl->inlen - off);
        off2 += svb_size((const char*)h, bl->inlen2 - off2);
    }
}

int cmd_extract(const char* in, const char* out, int t, int l) {
    std::string path, hpath = std::string(in) + ".hpc";
    if (!out) { path = std::string(in) + ".exc"; out = path.c_str(); }
    LCP_INIT();
    FILE* fi = std::fopen(in, "rb");
    FILE* fh = std::fopen(hpath.c_str(), "rb");
    bin_header h;
    std::fread(&h, 1, 32, fi);
    std::fseek(fh, 32, SEEK_SET);
    ectx c = {(scr*)std::calloc(t, sizeof(scr)), (uint64_t*)std::calloc(t, 8), l, l + 4};
    FILE* fo = open_tmp(out);
    uint64_t n = stream_run(fi, fh, fo, 32, t, svb_size, svb_size, ex_work, &c);
    std::fclose(fi);
    std::fclose(fh);
    uint64_t cores = 0;
    for (int i = 0; i < t; i++) {
        scr* s = &c.s[i];
        cores += c.nc[i];
        for (int j = 0; j < 2; j++) { std::free(s->rep[j]); std::free(s->rep32[j]); std::free(s->st[j]); std::free(s->en[j]); std::free(s->sz[j]); }
        std::free(s->lab); std::free(s->lsz); std::free(s->map); std::free(s->hs);
    }
    std::free(c.s);
    std::free(c.nc);
    bin_header o = {{'S', 'V', 'B', 'X'}, 3, n, cores, {0}};
    o.resv[0] = (uint8_t)l;
    o.resv[1] = (uint8_t)(l + 4);
    std::fseek(fo, 0, SEEK_SET);
    std::fwrite(&o, 1, 32, fo);
    commit_file(fo, out);
    return 0;
}
