#include <pthread.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "core.h"
#include "count.h"
#include "io.h"

static const uint32_t MM = 3, MARGIN = 10, NS = 65536, M31 = (1u << 31) - 1;

struct cctx { rdcore* r; uint64_t *tab, *pr, *hs; uint16_t* cn; uint32_t *sc, *deg; std::vector<uint64_t> ev[64]; uint64_t z, np, K; int w, t; };
struct job { cctx* x; int i; };

static uint64_t home(const cctx* x, uint64_t h) { return (h >> 48) * x->z + (((h >> 32) & 0xffff) * x->z >> 16); }

static uint64_t phome(const cctx* x, uint64_t k) { return (core_hash(&k, 32) >> 32) * x->np >> 32; }

static uint64_t add(cctx* x, uint64_t h) {
    uint64_t s = (h >> 48) * x->z, tag = h << 16 >> 16;
    for (uint64_t i = home(x, h);; i = i + 1 == s + x->z ? s : i + 1) {
        uint64_t v = __atomic_load_n(&x->tab[i], __ATOMIC_RELAXED);
        while (!v || v >> 16 == tag) {
            if ((v & 0xffff) == 0xffff) return i;
            if (__atomic_compare_exchange_n(&x->tab[i], &v, v ? v + 1 : tag << 16 | 1, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return i;
        }
    }
}

static void padd(cctx* x, uint64_t k) {
    for (uint64_t i = phome(x, k);; i = i + 1 == x->np ? 0 : i + 1) {
        uint64_t v = 0;
        if (__atomic_compare_exchange_n(&x->pr[i], &v, k, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED) || v == k) return;
    }
}

static void chain(cctx* x, rdcore* r, uint32_t lo, uint32_t hi, uint64_t* h, uint64_t* q, uint64_t* wd) {
    uint32_t m = 0, e = 0;
    for (uint32_t i = lo; i < hi; i++)
        if (r->a[i] >= MARGIN && r->b[i] + MARGIN <= r->len) h[m++] = core_hash(wd, core_cut(r, i, wd));
    for (uint32_t i = 0; i < m; i++) {
        if (i + 16 < m) __builtin_prefetch(&x->tab[home(x, h[i + 16])], 1);
        h[i] = add(x, h[i]);
    }
    for (uint32_t i = 1; i < m; i++) if (h[i - 1] != h[i]) q[e++] = 1ull << 63 | h[i - 1] << 31 | h[i];
    for (uint32_t i = 0; i < e; i++) {
        if (i + 16 < e) __builtin_prefetch(&x->pr[phome(x, q[i + 16])], 1);
        padd(x, q[i]);
    }
}

static void cnt_work(block* b) {
    cctx* x = (cctx*)b->ctx;
    rdcore* r = &x->r[b->id];
    size_t o2 = 0;
    for (size_t o = 0; o < b->inlen; o += svb_size(b->in + o, b->inlen - o)) {
        const uint8_t* u = (const uint8_t*)b->in2 + o2;
        collapse(r, (const uint8_t*)b->in + o);
        load_cores(r, u, x->w);
        bgrow(b, 8 * (2 * r->n + r->k / 32 + 2));
        uint64_t* h = (uint64_t*)b->buf;
        chain(x, r, 0, r->nf, h, h + r->n, h + 2 * r->n);
        chain(x, r, r->nf, r->n, h, h + r->n, h + 2 * r->n);
        o2 += exc_size((const char*)u, b->inlen2 - o2);
    }
}

static void* shard_count(void* a) {
    job* j = (job*)a;
    cctx* x = j->x;
    for (uint64_t s = NS * j->i / x->t; s < NS * (j->i + 1) / x->t; s++)
        for (uint64_t k = s * x->z; k < (s + 1) * x->z; k++) x->sc[s + 1] += (x->tab[k] & 0xffff) >= MM;
    return nullptr;
}

static void* shard_emit(void* a) {
    job* j = (job*)a;
    cctx* x = j->x;
    uint64_t* v = (uint64_t*)std::malloc(x->z * 8);
    for (uint64_t s = NS * j->i / x->t; s < NS * (j->i + 1) / x->t; s++) {
        uint64_t m = 0, *t = x->tab + s * x->z;
        for (uint64_t k = 0; k < x->z; k++) {
            if ((t[k] & 0xffff) >= MM) v[m++] = t[k] >> 16 << 16 | k;
            else t[k] = 0;
        }
        std::sort(v, v + m);
        for (uint64_t i = 0; i < m; i++) {
            uint64_t k = v[i] & 0xffff, rk = x->sc[s] + i;
            x->hs[rk] = s << 48 | v[i] >> 16;
            x->cn[rk] = (uint16_t)t[k];
            t[k] = rk + 1;
        }
    }
    std::free(v);
    return nullptr;
}

static void* edge_pick(void* a) {
    job* j = (job*)a;
    cctx* x = j->x;
    uint64_t lo = NS * j->i / x->t * x->z, hi = NS * (j->i + 1) / x->t * x->z;
    std::vector<uint64_t>& ev = x->ev[j->i];
    for (uint64_t i = 0; i < x->np; i++) {
        uint64_t k = x->pr[i], s1 = k >> 31 & M31;
        if (!k || s1 < lo || s1 >= hi) continue;
        uint64_t r1 = x->tab[s1], r2 = x->tab[k & M31];
        if (r1 && r2) { ev.push_back((r1 - 1) << 32 | (r2 - 1)); x->deg[r1 - 1]++; }
    }
    std::sort(ev.begin(), ev.end());
    return nullptr;
}

static void par(cctx* x, void* (*f)(void*)) {
    pthread_t th[64];
    job j[64];
    for (int i = 0; i < x->t; i++) { j[i] = {x, i}; pthread_create(&th[i], nullptr, f, &j[i]); }
    for (int i = 0; i < x->t; i++) pthread_join(th[i], nullptr);
}

static int nbits(uint64_t v) { return v ? 64 - __builtin_clzll(v) : 1; }

int cmd_count(const char* in, const char* out, int t) {
    std::string path, epath = std::string(in) + ".exc";
    if (!out) { path = std::string(in) + ".cnt"; out = path.c_str(); }
    FILE* fi = std::fopen(in, "rb");
    FILE* fe = std::fopen(epath.c_str(), "rb");
    bin_header h, e;
    std::fread(&h, 1, 32, fi);
    std::fread(&e, 1, 32, fe);
    cctx x = {};
    x.r = (rdcore*)std::calloc(t, sizeof(rdcore));
    x.z = (e.total / 24 + 120000000) / NS + 64;
    x.np = NS * x.z * 5 / 4;
    x.w = e.resv[1];
    x.t = t;
    x.tab = (uint64_t*)std::calloc(NS * x.z, 8);
    x.pr = (uint64_t*)std::calloc(x.np, 8);
    FILE* fo = open_tmp(out);
    stream_run(fi, fe, fo, 32, t, svb_size, exc_size, cnt_work, &x);
    std::fclose(fi);
    std::fclose(fe);
    for (int i = 0; i < t; i++) free_core(&x.r[i]);
    x.sc = (uint32_t*)std::calloc(NS + 1, 4);
    par(&x, shard_count);
    for (uint32_t s = 0; s < NS; s++) x.sc[s + 1] += x.sc[s];
    x.K = x.sc[NS];
    x.hs = (uint64_t*)std::malloc(x.K * 8);
    x.cn = (uint16_t*)std::malloc(x.K * 2);
    par(&x, shard_emit);
    x.deg = (uint32_t*)std::calloc(x.K, 4);
    par(&x, edge_pick);
    std::free(x.tab);
    std::free(x.pr);
    uint64_t E = 0, md = 0;
    for (uint64_t u = 0; u < x.K; u++) { E += x.deg[u]; md = std::max<uint64_t>(md, x.deg[u]); }
    int rb = nbits(x.K - 1), db = nbits(md);
    std::fwrite(x.hs, 8, x.K, fo);
    std::fwrite(x.cn, 2, x.K, fo);
    block ob = {};
    bits bw = {};
    for (uint64_t u = 0; u < x.K; u++) wbits(&ob, &bw, x.deg[u], db);
    for (int i = 0; i < t; i++) for (uint64_t e : x.ev[i]) wbits(&ob, &bw, (uint32_t)e, rb);
    wflush(&ob, &bw);
    std::fwrite(ob.out, 1, ob.outlen, fo);
    bin_header o = {{'S', 'V', 'B', 'C'}, 1, x.K, E, {(uint8_t)rb, (uint8_t)db, MM, MARGIN}};
    std::fseek(fo, 0, SEEK_SET);
    std::fwrite(&o, 1, 32, fo);
    commit_file(fo, out);
    std::free(ob.out);
    std::free(x.hs);
    std::free(x.cn);
    std::free(x.sc);
    std::free(x.deg);
    std::free(x.r);
    return 0;
}
