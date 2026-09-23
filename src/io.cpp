#include <fcntl.h>
#include <pthread.h>
#include <sys/mman.h>
#include <unistd.h>
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

const uint8_t* map_file(const char* path) {
    int fd = open(path, O_RDONLY);
    const uint8_t* p = (const uint8_t*)mmap(nullptr, lseek(fd, 0, SEEK_END), PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    return p;
}

FILE* open_tmp(const char* path) {
    std::string tmp = std::string(path) + ".tmp";
    FILE* fo = std::fopen(tmp.c_str(), "wb+");
    std::setvbuf(fo, nullptr, _IOFBF, 1 << 22);
    return fo;
}

void commit_file(FILE* fo, const char* path) {
    std::fclose(fo);
    std::rename((std::string(path) + ".tmp").c_str(), path);
}

void bput(block* b, const void* d, size_t n) {
    if (b->outlen + n > b->outcap) {
        b->outcap = (b->outlen + n) * 2;
        b->out = (char*)std::realloc(b->out, b->outcap);
    }
    std::memcpy(b->out + b->outlen, d, n);
    b->outlen += n;
}

void bgrow(block* b, size_t n) {
    if (n > b->bufcap) { b->bufcap = 2 * n + 64; b->buf = (char*)std::realloc(b->buf, b->bufcap); }
}

void wbits(block* b, bits* s, uint64_t v, int w) {
    s->acc = (s->acc << w) | v;
    s->pos += w;
    while (s->pos >= 8) {
        s->pos -= 8;
        char c = (char)(s->acc >> s->pos);
        bput(b, &c, 1);
    }
}

void wflush(block* b, bits* s) {
    if (!s->pos) return;
    char c = (char)(s->acc << (8 - s->pos));
    bput(b, &c, 1);
    s->pos = 0;
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

static work_fn g_work;

static void* run_block(void* a) {
    g_work((block*)a);
    return nullptr;
}

struct part { int fd; char* p; size_t n, got; long o; };

static void* read_part(void* a) {
    part* r = (part*)a;
    for (ssize_t g = 1; g > 0 && r->got < r->n; r->got += g) g = pread(r->fd, r->p + r->got, r->n - r->got, r->o + r->got);
    return nullptr;
}

static size_t pread8(char* p, size_t n, FILE* f) {
    long o = std::ftell(f);
    part r[8];
    pthread_t th[8];
    size_t got = 0;
    for (int i = 0; i < 8; i++) {
        r[i] = {fileno(f), p + n * i / 8, n * (i + 1) / 8 - n * i / 8, 0, o + (long)(n * i / 8)};
        pthread_create(&th[i], nullptr, read_part, &r[i]);
    }
    for (int i = 0; i < 8; i++) { pthread_join(th[i], nullptr); got += r[i].got; }
    std::fseek(f, o + (long)got, SEEK_SET);
    return got;
}

struct filler { FILE* f1; FILE* f2; char* b1; char* b2; size_t k1, k2, cap, l1, l2; };

static void* fill(void* a) {
    filler* f = (filler*)a;
    f->l1 = f->k1 + pread8(f->b1 + f->k1, f->cap - f->k1, f->f1);
    if (f->f2) f->l2 = f->k2 + pread8(f->b2 + f->k2, f->cap - f->k2, f->f2);
    return nullptr;
}

uint64_t stream_run(FILE* fi, FILE* f2, FILE* fo, size_t hdrsz, int t, size_fn sf, size_fn s2, work_fn wf, void* ctx) {
    if (t > 64) t = 64;
    g_work = wf;
    std::fseek(fo, (long)hdrsz, SEEK_SET);
    size_t icap = IBUF;
    char* ib[2] = {(char*)std::malloc(icap), (char*)std::malloc(icap)};
    char* jb[2] = {nullptr, nullptr};
    if (f2) { jb[0] = (char*)std::malloc(icap); jb[1] = (char*)std::malloc(icap); }
    block* b = (block*)std::calloc(t, sizeof(block));
    pthread_t* th = (pthread_t*)std::malloc(t * sizeof(pthread_t));
    pthread_t rd;
    for (int i = 0; i < t; i++) { b[i].id = i; b[i].ctx = ctx; }
    int cur = 0;
    uint64_t nrec = 0;
    size_t ilen = pread8(ib[0], icap, fi), jlen = 0;
    if (f2) jlen = pread8(jb[0], icap, f2);
    while (ilen) {
        size_t target = ilen / t + 1, off = 0, start = 0, off2 = 0, start2 = 0;
        int nb = 0;
        for (;;) {
            size_t r = sf(ib[cur] + off, ilen - off);
            if (!r) break;
            if (f2) {
                size_t r2 = s2(jb[cur] + off2, jlen - off2);
                if (!r2) break;
                off2 += r2;
            }
            off += r;
            nrec++;
            if (off - start >= target && nb < t - 1) {
                b[nb].in = ib[cur] + start;
                b[nb].inlen = off - start;
                b[nb].in2 = f2 ? jb[cur] + start2 : nullptr;
                b[nb].inlen2 = off2 - start2;
                nb++;
                start = off;
                start2 = off2;
            }
        }
        if (off > start) {
            b[nb].in = ib[cur] + start;
            b[nb].inlen = off - start;
            b[nb].in2 = f2 ? jb[cur] + start2 : nullptr;
            b[nb].inlen2 = off2 - start2;
            nb++;
        }
        if (!nb) {
            if (ilen < icap && (!f2 || jlen < icap)) break;
            icap *= 2;
            ib[0] = (char*)std::realloc(ib[0], icap);
            ib[1] = (char*)std::realloc(ib[1], icap);
            ilen += pread8(ib[cur] + ilen, icap - ilen, fi);
            if (f2) {
                jb[0] = (char*)std::realloc(jb[0], icap);
                jb[1] = (char*)std::realloc(jb[1], icap);
                jlen += pread8(jb[cur] + jlen, icap - jlen, f2);
            }
            continue;
        }
        filler f = {fi, f2, ib[cur ^ 1], jb[cur ^ 1], ilen - off, jlen - off2, icap, 0, 0};
        std::memcpy(f.b1, ib[cur] + off, f.k1);
        if (f2) std::memcpy(f.b2, jb[cur] + off2, f.k2);
        pthread_create(&rd, nullptr, fill, &f);
        for (int i = 0; i < nb; i++) pthread_create(&th[i], nullptr, run_block, &b[i]);
        for (int i = 0; i < nb; i++) pthread_join(th[i], nullptr);
        for (int i = 0; i < nb; i++) { std::fwrite(b[i].out, 1, b[i].outlen, fo); b[i].outlen = 0; }
        pthread_join(rd, nullptr);
        cur ^= 1;
        ilen = f.l1;
        jlen = f.l2;
    }
    for (int i = 0; i < t; i++) { std::free(b[i].out); std::free(b[i].buf); }
    std::free(b);
    std::free(th);
    std::free(ib[0]);
    std::free(ib[1]);
    std::free(jb[0]);
    std::free(jb[1]);
    return nrec;
}
