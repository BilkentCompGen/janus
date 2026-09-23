#include <immintrin.h>
#include <cstdlib>
#include <cstring>
#include "core.h"
#include "io.h"

__extension__ typedef unsigned __int128 u128;

static uint64_t be64(const uint8_t* p) {
    uint64_t x;
    std::memcpy(&x, p, 8);
    return __builtin_bswap64(x);
}

static uint32_t rbits(const uint8_t* p, uint64_t* pos, int w) {
    uint32_t v = (uint32_t)((be64(p + (*pos >> 3)) << (*pos & 7)) >> (64 - w));
    *pos += w;
    return v;
}

static uint64_t rev2(uint64_t x) {
    x = __builtin_bswap64(x);
    x = ((x >> 4) & 0x0f0f0f0f0f0f0f0full) | ((x & 0x0f0f0f0f0f0f0f0full) << 4);
    return ((x >> 2) & 0x3333333333333333ull) | ((x & 0x3333333333333333ull) << 2);
}

static uint32_t rank(const rdcore* r, uint32_t x) {
    uint32_t i = x & 31;
    return r->cum[x >> 5] + __builtin_popcountll(r->mk[x >> 5] & (i ? ~0ull << (64 - 2 * i) : 0));
}

static void rgrow(rdcore* r, uint32_t len) {
    if (len <= r->cap) return;
    r->cap = 2 * len;
    size_t nw = r->cap / 32 + 4;
    r->raw = (uint8_t*)std::realloc(r->raw, nw * 8);
    r->hs = (uint64_t*)std::realloc(r->hs, nw * 8);
    r->rs = (uint64_t*)std::realloc(r->rs, nw * 8);
    r->mk = (uint64_t*)std::realloc(r->mk, nw * 8);
    r->cum = (uint32_t*)std::realloc(r->cum, nw * 4);
}

size_t exc_size(const char* p, size_t n) {
    if (n < 6) return 0;
    size_t sz = (uint8_t)p[4] << 8 | (uint8_t)p[5];
    return sz <= n ? sz : 0;
}

void collapse(rdcore* r, const uint8_t* rec) {
    uint32_t tag;
    std::memcpy(&tag, rec, 4);
    uint32_t len = tag & 0x7fffffff, nw = (len + 31) / 32, k = 0;
    rgrow(r, len);
    std::memset(r->raw, 0, nw * 8 + 8);
    if (tag >> 31) for (uint32_t i = 0; i < len; i++) r->raw[i >> 2] |= CODE.t[rec[4 + i] & 0xdf] << (6 - 2 * (i & 3));
    else std::memcpy(r->raw, rec + 4, (len + 3) / 4);
    std::memset(r->hs, 0, nw * 8 + 16);
    uint64_t prev = ~(uint64_t)r->raw[0] >> 6 & 3;
    for (uint32_t j = 0; j < nw; j++) {
        uint64_t x = be64(r->raw + 8 * j), d = x ^ (x >> 2 | prev << 62), s = (d | d >> 1) & 0x5555555555555555ull;
        if (len - 32 * j < 32) s &= ~0ull << (64 - 2 * (len - 32 * j));
        prev = x & 3;
        r->mk[j] = s;
        r->cum[j] = k;
        uint32_t n = __builtin_popcountll(s), o = 2 * (k & 31);
        uint64_t v = n ? _pext_u64(x, s | s << 1) << (64 - 2 * n) : 0;
        r->hs[k >> 5] |= v >> o;
        if (o) r->hs[(k >> 5) + 1] |= v << (64 - o);
        k += n;
    }
    r->mk[nw] = 0;
    r->cum[nw] = k;
    r->len = len;
    r->k = k;
}

void load_cores(rdcore* r, const uint8_t* unit, int w) {
    uint32_t nw = (r->k + 31) / 32, m = (1u << w) - 1, prev = 0;
    for (uint32_t j = 0; j < nw; j++) r->rs[j] = rev2(~r->hs[nw - 1 - j]);
    r->rs[nw] = 0;
    r->pad = 32 * nw - r->k;
    r->n = unit[0] << 8 | unit[1];
    r->nf = unit[2] << 8 | unit[3];
    size_t ub = unit[4] << 8 | unit[5];
    if (r->n + 1 > r->ccap || ub + 8 > 8 * (size_t)r->ccap) {
        r->ccap = 2 * (r->n + ub) + 16;
        r->a = (uint32_t*)std::realloc(r->a, r->ccap * 8);
        r->b = r->a + r->ccap;
        r->u = (uint8_t*)std::realloc(r->u, r->ccap * 8);
    }
    std::memcpy(r->u, unit, ub);
    std::memset(r->u + ub, 0, 8);
    uint64_t pos = 48;
    for (uint32_t i = 0; i < r->n; i++) {
        uint32_t d = rbits(r->u, &pos, w), l;
        if (d == m) { d = rbits(r->u, &pos, 32); l = rbits(r->u, &pos, 32); }
        else l = rbits(r->u, &pos, w);
        if (i == r->nf) prev = 0;
        r->a[i] = prev + d;
        r->b[i] = r->a[i] + l;
        prev = r->a[i];
    }
}

uint32_t core_cut(const rdcore* r, uint32_t i, uint64_t* out) {
    uint32_t st, en;
    const uint64_t* W;
    if (i < r->nf) { st = rank(r, r->a[i]); en = rank(r, r->b[i]); W = r->hs; }
    else { st = r->k - rank(r, r->len - r->a[i]) + r->pad; en = r->k - rank(r, r->len - r->b[i]) + r->pad; W = r->rs; }
    uint32_t l = en - st, n = (l + 31) / 32, q = st >> 5, s = 2 * (st & 31);
    for (uint32_t j = 0; j < n; j++) out[j] = s ? W[q + j] << s | W[q + j + 1] >> (64 - s) : W[q + j];
    if (l & 31) out[n - 1] &= ~0ull << (64 - 2 * (l & 31));
    return l;
}

static uint64_t mix(uint64_t a, uint64_t b) {
    u128 p = (u128)a * b;
    return (uint64_t)p ^ (uint64_t)(p >> 64);
}

uint64_t core_hash(const uint64_t* w, uint32_t k) {
    uint64_t h = mix(k ^ 0xa0761d6478bd642full, 0xe7037ed1a0b428dbull);
    for (uint32_t i = 0; i < (k + 31) / 32; i++) h = mix(h ^ w[i], 0x9e3779b97f4a7c15ull);
    return mix(h ^ k, 0xd6e8feb86659fd93ull);
}

void free_core(rdcore* r) {
    std::free(r->raw);
    std::free(r->hs);
    std::free(r->rs);
    std::free(r->mk);
    std::free(r->cum);
    std::free(r->a);
    std::free(r->u);
}
