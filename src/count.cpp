#include <unistd.h>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>
#include "cmd.h"
#include "core.h"
#include "io.h"

static const uint32_t MM = 3; // occurrences a core needs to become a node
static const uint32_t MARGIN = 10; // cores closer than this to either read end are not counted
static const uint32_t NS = 65536; // table shards, one per value of the hash's top 16 bits
static const uint64_t ZDIV = 24; // table size estimate: one slot per ZDIV core occurrences of the .exc
static const uint64_t ZADD = 120000000; // plus this many slots
static const uint64_t ZPAD = 64; // plus this many slots per shard
static const uint64_t M31 = (1u << 31) - 1; // a table slot in a pair key

struct cth { rdcore r; std::vector<uint64_t> h; std::vector<std::pair<uint64_t, uint16_t>> ev; };

struct cctx { cth* c; uint64_t *tab, *pr; uint16_t* pc; uint32_t* sc; std::vector<uint32_t> deg; uint64_t z, np, K; int w, t, fd; };

static uint64_t home(const cctx* x, uint64_t h) { return (h >> 48) * x->z + (((h >> 32) & 0xffff) * x->z >> 16); }

static uint64_t phome(const cctx* x, uint64_t k) { return (core_hash(&k, 32) >> 32) * x->np >> 32; }

static uint64_t shard(const cctx* x, int i) { return NS * i / x->t; }

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
        if (__atomic_compare_exchange_n(&x->pr[i], &v, k, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED) || v == k) {
            if (x->pc[i] < 0xffff) __atomic_fetch_add(&x->pc[i], 1, __ATOMIC_RELAXED);
            return;
        }
    }
}

static void chain(cctx* x, rdcore* r, uint32_t lo, uint32_t hi, uint64_t* h, uint64_t* q, uint64_t* wd) {
    uint32_t m = 0, e = 0;
    for (uint32_t i = lo; i < hi; i++)
        if (r->a[i] >= MARGIN && r->b[i] + MARGIN <= r->len) h[m++] = core_hash(wd, core_cut(r, i, wd));
    for (uint32_t i = 0; i < m; i++) {
        if (i + PF < m) __builtin_prefetch(&x->tab[home(x, h[i + PF])], 1);
        h[i] = add(x, h[i]);
    }
    for (uint32_t i = 1; i < m; i++) if (h[i - 1] != h[i]) q[e++] = 1ull << 63 | h[i - 1] << 31 | h[i];
    for (uint32_t i = 0; i < e; i++) {
        if (i + PF < e) __builtin_prefetch(&x->pr[phome(x, q[i + PF])], 1);
        padd(x, q[i]);
    }
}

static void cnt_work(block* b) {
    cctx* x = (cctx*)b->ctx;
    cth* c = &x->c[b->id];
    for (size_t o = 0, e = 0; o < b->inlen[0];) {
        next_read(&c->r, b, &o, &e, x->w);
        uint32_t n = c->r.n;
        if (c->h.size() < 2 * n + c->r.k / 32 + 2) c->h.resize(4 * n + c->r.k / 16 + 4);
        chain(x, &c->r, 0, c->r.nf, c->h.data(), c->h.data() + n, c->h.data() + 2 * n);
        chain(x, &c->r, c->r.nf, n, c->h.data(), c->h.data() + n, c->h.data() + 2 * n);
    }
}

static void shard_count(void* a, int i) {
    cctx* x = (cctx*)a;
    for (uint64_t s = shard(x, i); s < shard(x, i + 1); s++)
        for (uint64_t k = s * x->z; k < (s + 1) * x->z; k++) x->sc[s + 1] += (x->tab[k] & 0xffff) >= MM;
}

static void shard_emit(void* a, int i) {
    cctx* x = (cctx*)a;
    std::vector<uint64_t> v(x->z);
    std::vector<uint16_t> n(x->z);
    for (uint64_t s = shard(x, i); s < shard(x, i + 1); s++) {
        uint64_t m = 0, *t = x->tab + s * x->z;
        for (uint64_t k = 0; k < x->z; k++) {
            if ((t[k] & 0xffff) >= MM) v[m++] = t[k] >> 16 << 16 | k;
            else t[k] = 0;
        }
        std::sort(v.begin(), v.begin() + m);
        for (uint64_t j = 0; j < m; j++) {
            uint64_t k = v[j] & 0xffff;
            n[j] = (uint16_t)t[k];
            v[j] = s << 48 | v[j] >> 16;
            t[k] = x->sc[s] + j + 1;
        }
        pwrite(x->fd, v.data(), m * 8, sizeof(bin_header) + x->sc[s] * 8);
        pwrite(x->fd, n.data(), m * 2, sizeof(bin_header) + x->K * 8 + x->sc[s] * 2);
    }
}

static void edge_pick(void* a, int i) {
    cctx* x = (cctx*)a;
    std::vector<std::pair<uint64_t, uint16_t>>& ev = x->c[i].ev;
    for (uint64_t j = 0; j < x->np; j++) {
        uint64_t k = x->pr[j], s1 = k >> 31 & M31;
        if (!k || s1 < shard(x, i) * x->z || s1 >= shard(x, i + 1) * x->z) continue;
        uint64_t r1 = x->tab[s1], r2 = x->tab[k & M31];
        if (r1 && r2) { ev.push_back({(r1 - 1) << 32 | (r2 - 1), x->pc[j]}); x->deg[r1 - 1]++; }
    }
    std::sort(ev.begin(), ev.end());
}

static int nbits(uint64_t v) { return v ? 64 - __builtin_clzll(v) : 1; }

int cmd_count(const char* in, const char* out, int t) {
    bin_header h, e;
    std::vector<cth> c(t);
    std::vector<uint32_t> sc(NS + 1);
    std::vector<char> ob;
    FILE* fi = open_in(in, &h), *fe = open_in((std::string(in) + ".exc").c_str(), &e), *fo = open_tmp(out);
    cctx x = {c.data(), nullptr, nullptr, nullptr, sc.data(), {}, (e.total / ZDIV + ZADD) / NS + ZPAD, 0, 0, e.resv[1], t, fileno(fo)};
    x.np = NS * x.z * 5 / 4;
    x.tab = (uint64_t*)std::calloc(NS * x.z, 8);
    x.pr = (uint64_t*)std::calloc(x.np, 8);
    x.pc = (uint16_t*)std::calloc(x.np, 2);
    stream_run(fi, fe, fo, t, cnt_work, &x);
    par(t, shard_count, &x);
    for (uint32_t s = 0; s < NS; s++) sc[s + 1] += sc[s];
    x.K = sc[NS];
    x.deg.resize(x.K);
    par(t, shard_emit, &x);
    par(t, edge_pick, &x);
    std::free(x.tab);
    std::free(x.pr);
    std::free(x.pc);
    uint64_t E = 0, md = 0;
    for (uint32_t d : x.deg) { E += d; md = std::max<uint64_t>(md, d); }
    bits bw = {&ob, 0, 0};
    int rb = nbits(x.K - 1), db = nbits(md);
    for (uint32_t d : x.deg) wbits(&bw, d, db);
    std::vector<uint16_t> ec;
    for (cth& q : c) for (auto& v : q.ev) { wbits(&bw, (uint32_t)v.first, rb); ec.push_back(v.second); }
    wflush(&bw);
    std::fseek(fo, sizeof(bin_header) + x.K * 10, SEEK_SET);
    std::fwrite(ec.data(), 2, E, fo);
    std::fwrite(ob.data(), 1, ob.size(), fo);
    bin_header o = header('C', 1, x.K, E);
    o.resv[0] = rb; o.resv[1] = db; o.resv[2] = MM; o.resv[3] = MARGIN;
    finish(fo, out, &o);
    return 0;
}
