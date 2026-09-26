#include <pthread.h>
#include <unistd.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <string>
#include "io.h"

static const size_t IBUF = 1 << 26;

code_table::code_table() {
    for (int i = 0; i < 256; i++) t[i] = 3;
    t[(uint8_t)'A'] = 0; t[(uint8_t)'C'] = 1; t[(uint8_t)'G'] = 2; t[(uint8_t)'T'] = 3;
}

const code_table CODE;

FILE* open_in(const char* path, bin_header* h) {
    FILE* f = std::fopen(path, "rb");
    std::fread(h, 1, sizeof(bin_header), f);
    return f;
}

FILE* open_tmp(const char* path) {
    FILE* fo = std::fopen((std::string(path) + ".tmp").c_str(), "wb+");
    std::setvbuf(fo, nullptr, _IOFBF, 1 << 22);
    return fo;
}

void finish(FILE* fo, const char* path, const bin_header* h) {
    std::fseek(fo, 0, SEEK_SET);
    std::fwrite(h, 1, sizeof(bin_header), fo);
    std::fclose(fo);
    std::rename((std::string(path) + ".tmp").c_str(), path);
}

struct pjob { par_fn f; void* ctx; int i; };

static void* prun(void* a) {
    pjob* j = (pjob*)a;
    j->f(j->ctx, j->i);
    return nullptr;
}

void par(int n, par_fn f, void* ctx) {
    pthread_t th[65];
    pjob j[65];
    for (int i = 0; i < n; i++) { j[i] = {f, ctx, i}; pthread_create(&th[i], nullptr, prun, &j[i]); }
    for (int i = 0; i < n; i++) pthread_join(th[i], nullptr);
}

uint64_t be64(const uint8_t* p) {
    uint64_t x;
    std::memcpy(&x, p, 8);
    return __builtin_bswap64(x);
}

uint32_t rbits(const uint8_t* p, uint64_t* pos, int w) {
    uint32_t v = (uint32_t)((be64(p + (*pos >> 3)) << (*pos & 7)) >> (64 - w));
    *pos += w;
    return v;
}

void wbits(bits* s, uint64_t v, int w) {
    s->acc = s->acc << w | v;
    for (s->pos += w; s->pos >= 8; s->pos -= 8) s->out->push_back((char)(s->acc >> (s->pos - 8)));
}

void wflush(bits* s) {
    if (s->pos) wbits(s, 0, 8 - s->pos);
}

std::vector<std::string> split(const char* s) {
    std::vector<std::string> v;
    std::string p = s;
    for (size_t a = 0, b; a < p.size(); a = b + 1) {
        b = std::min(p.find(',', a), p.size());
        v.push_back(p.substr(a, b - a));
    }
    return v;
}

size_t nlines(const char* p, size_t n, int k) {
    size_t i = 0;
    while (k--) {
        const char* e = (const char*)std::memchr(p + i, '\n', n - i);
        if (!e) return 0;
        i = e - p + 1;
    }
    return i;
}

size_t svb_size(const char* p, size_t n) {
    if (n < 4) return 0;
    uint32_t tag;
    std::memcpy(&tag, p, 4);
    uint32_t len = tag & 0x7fffffff;
    size_t sz = 4 + (tag & 0x80000000 ? len : (len + 3) / 4);
    return sz <= n ? sz : 0;
}

struct part { int fd; char* p; size_t n, got[8]; long o; };

static void read_part(void* a, int i) {
    part* r = (part*)a;
    size_t at = r->n * i / 8, end = r->n * (i + 1) / 8;
    for (ssize_t g = 1; g > 0 && at < end; at += g) g = pread(r->fd, r->p + at, end - at, r->o + at);
    r->got[i] = at - r->n * i / 8;
}

static size_t pread8(char* p, size_t n, FILE* f) {
    part r = {fileno(f), p, n, {}, std::ftell(f)};
    par(8, read_part, &r);
    size_t got = 0;
    for (size_t g : r.got) got += g;
    std::fseek(f, r.o + (long)got, SEEK_SET);
    return got;
}

struct stream { FILE* f; size_fn sf; char* b[2]; size_t len, off, start, next; };

struct batch { stream* s; block* b; work_fn wf; int n, nb, cur; size_t cap; };

static void run(void* a, int i) {
    batch* r = (batch*)a;
    if (i < r->nb) { r->wf(&r->b[i]); return; }
    for (int j = 0; j < r->n; j++) {
        stream* s = &r->s[j];
        size_t k = s->len - s->off;
        std::memcpy(s->b[r->cur ^ 1], s->b[r->cur] + s->off, k);
        s->next = k + pread8(s->b[r->cur ^ 1] + k, r->cap - k, s->f);
    }
}

static void cut(block* b, stream* s, int n, int cur) {
    for (int i = 0; i < n; i++) {
        b->in[i] = s[i].b[cur] + s[i].start;
        b->inlen[i] = s[i].off - s[i].start;
        s[i].start = s[i].off;
    }
}

uint64_t stream_run(FILE* fi, FILE* f2, FILE* fo, int t, size_fn sf, size_fn s2, work_fn wf, void* ctx) {
    stream s[2] = {{fi, sf, {}, 0, 0, 0, 0}, {f2, s2, {}, 0, 0, 0, 0}};
    block b[64];
    batch r = {s, b, wf, f2 ? 2 : 1, 0, 0, IBUF};
    uint64_t nrec = 0;
    t = std::min(t, 64);
    std::fseek(fo, sizeof(bin_header), SEEK_SET);
    for (int i = 0; i < t; i++) { b[i].id = i; b[i].ctx = ctx; b[i].out.reserve(IBUF); }
    for (int i = 0; i < r.n; i++) {
        s[i].b[0] = (char*)std::malloc(IBUF + 8);
        s[i].b[1] = (char*)std::malloc(IBUF + 8);
        s[i].len = pread8(s[i].b[0], IBUF, s[i].f);
    }
    while (s[0].len) {
        int i;
        r.nb = 0;
        for (i = 0; i < r.n; i++) s[i].off = s[i].start = 0;
        for (;;) {
            size_t z[2];
            for (i = 0; i < r.n && (z[i] = s[i].sf(s[i].b[r.cur] + s[i].off, s[i].len - s[i].off)); i++);
            if (i < r.n) break;
            for (i = 0; i < r.n; i++) s[i].off += z[i];
            nrec++;
            if (s[0].off - s[0].start > s[0].len / t && r.nb < t - 1) cut(&b[r.nb++], s, r.n, r.cur);
        }
        if (s[0].off > s[0].start) cut(&b[r.nb++], s, r.n, r.cur);
        if (!r.nb) {
            for (i = 0; i < r.n && s[i].len < r.cap; i++);
            if (i == r.n) break;
            r.cap *= 2;
            for (i = 0; i < r.n; i++) {
                s[i].b[0] = (char*)std::realloc(s[i].b[0], r.cap + 8);
                s[i].b[1] = (char*)std::realloc(s[i].b[1], r.cap + 8);
                s[i].len += pread8(s[i].b[r.cur] + s[i].len, r.cap - s[i].len, s[i].f);
            }
            continue;
        }
        par(r.nb + 1, run, &r);
        for (i = 0; i < r.nb; i++) { std::fwrite(b[i].out.data(), 1, b[i].out.size(), fo); b[i].out.clear(); }
        r.cur ^= 1;
        for (i = 0; i < r.n; i++) s[i].len = s[i].next;
    }
    for (int i = 0; i < r.n; i++) { std::fclose(s[i].f); std::free(s[i].b[0]); std::free(s[i].b[1]); }
    return nrec;
}
