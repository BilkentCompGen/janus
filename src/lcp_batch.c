// Batch LCP parsing for janus: the level-1 parse, the deterministic coin tossing and the level-n parse of
// lcptools (biolcp branch) on struct-of-arrays buffers, scanning 64 positions per word (AVX-512 where available).
// Same cores as lcptools' init_lps/lps_deepen, about 2.7 times faster; built for lcptools configured with
// LABEL=0 POS=64 DCT=1 CORE=var ALPHABET=dna (see the Makefile).
#include <stdlib.h>
#include <string.h>
#if defined(__AVX512F__)
#include <immintrin.h>
#endif
#include "lcp_batch.h"

#if LCP_FIXED_CORE || LCP_LABEL_BITS != 0 || LCP_POS_BITS == 0 || LCP_DCT_ITERATION_COUNT != 1 || LCP_ALPHABET_PROTEIN
#error "lcp_batch.c needs lcptools built with LABEL=0 POS=64 DCT=1 CORE=var ALPHABET=dna"
#endif
#ifndef minimum
#define minimum(a, b) ((a) < (b) ? (a) : (b))
#endif

#if !defined(__AVX512BW__)
static const uint8_t PARSE1_ACGT[256] = { ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1 };
#endif

static inline size_t ones_from(const uint64_t *e, size_t nw, size_t p) {
    size_t cnt = 0, w = p >> 6;
    unsigned b = (unsigned)(p & 63);
    while (w < nw) {
        const uint64_t bits = e[w] >> b;
        if (bits == (~0ull >> b)) { cnt += 64 - b; w++; b = 0; continue; }
        return cnt + (size_t)__builtin_ctzll(~bits);
    }
    return cnt;
}


static inline uint32_t bl32(uint32_t o) {
    uint32_t r = 0, t;
    t = (o >> 16) != 0u; r += t << 4; o >>= t << 4;
    t = (o >> 8) != 0u;  r += t << 3; o >>= t << 3;
    t = (o >> 4) != 0u;  r += t << 2; o >>= t << 2;
    t = (o >> 2) != 0u;  r += t << 1; o >>= t << 1;
    t = (o >> 1) != 0u;  r += t;      o >>= t;
    return r + o;
}

#if defined(__AVX512CD__)
static inline uint32_t ctz4(uint32_t x) { return (uint32_t)__builtin_ctz(x | 16u); }
static inline uint32_t ctz2(uint32_t y) { return (uint32_t)__builtin_ctz(y | 4u); }
static inline uint32_t bitlen32(uint32_t o) { return o ? 32u - (uint32_t)__builtin_clz(o) : 0u; }
#else
static inline uint32_t ctz4(uint32_t x) {
    return ((x & 1u) == 0u) + ((x & 3u) == 0u) + ((x & 7u) == 0u) + ((x & 15u) == 0u);
}
static inline uint32_t ctz2(uint32_t y) {
    return ((y & 1u) == 0u) + ((y & 3u) == 0u);
}
static inline uint32_t bitlen32(uint32_t o) { return bl32(o); }
#endif

static inline uint32_t dct1_out32(uint32_t L, uint32_t R, uint32_t Rsz) {
    const uint32_t x = (L ^ R) & 15u;
    const uint32_t L2 = (L >> 2) & 3u, L1 = (L >> 4) & 3u;
    const uint32_t R2 = (R >> 2) & 3u, R1 = (R >> 4) & 3u;
    const uint32_t Lm = L >> LCP_LENGTH_SHIFT, Rm = R >> LCP_LENGTH_SHIFT;
    const uint32_t kx = ctz4(x);
    const uint32_t o32 = 2u * kx + ((R >> kx) & 1u);
    const uint32_t lt = 0u - (uint32_t)(Lm < Rm);
    const uint32_t a2 = (L1 & lt) | (L2 & ~lt);
    const uint32_t b2 = (R2 & lt) | (R1 & ~lt);
    const uint32_t mn = (Lm & lt) | (Rm & ~lt);
    const uint32_t km = ctz2(a2 ^ b2);
    const uint32_t om = 4u * (mn + 1u) + 2u * km + ((b2 >> km) & 1u);
    const uint32_t k1 = ctz2(L1 ^ R1);
    const uint32_t o1 = 4u * (Lm + 1u) + 2u * k1 + ((R1 >> k1) & 1u);
    const uint32_t o0 = 2u * Rsz;
    const uint32_t m1 = 0u - (uint32_t)(L1 != R1);
    const uint32_t mm = 0u - (uint32_t)(Lm != Rm);
    const uint32_t mx = 0u - (uint32_t)(x != 0u);
    uint32_t out = (o1 & m1) | (o0 & ~m1);
    out = (om & mm) | (out & ~mm);
    out = (o32 & mx) | (out & ~mx);
    return out;
}

#if defined(__AVX512CD__)
static inline uint32_t dctu_out(uint64_t L, uint32_t Lsz, uint64_t R, uint32_t Rsz) {
    const uint64_t bound = minimum(Rsz, minimum(Lsz, 64u));
    const uint64_t x = L ^ R;
    const uint64_t eq = -(uint64_t)(x == 0);
    const uint64_t cz = (uint64_t)__builtin_ctzll(x | eq);
    const uint64_t idx = (bound & eq) | (minimum(cz, bound) & ~eq);
    return (uint32_t)(2ull * idx + ((R >> idx) & 1ull));
}
#else
static inline uint32_t dctu_out(uint64_t L, uint32_t Lsz, uint64_t R, uint32_t Rsz) {
    const uint64_t bound = minimum(Rsz, minimum(Lsz, 64u));
    uint64_t idx = bound;
    if (L != R) idx = minimum((uint64_t)__builtin_ctzll(L ^ R), bound);
    return (uint32_t)(2ull * idx + ((R >> idx) & 1ull));
}
#endif

static inline uint32_t bitlen_min2_32(uint32_t x) {
    const uint32_t bl = bitlen32(x);
    return bl < 2u ? 2u : bl;
}

static inline uint32_t bitlen_min2_clz(uint32_t x) {
    const uint32_t bl = x ? 32u - (uint32_t)__builtin_clz(x) : 0u;
    return bl < 2u ? 2u : bl;
}

static inline void soa1_core(struct lcp_soa *o, size_t idx, const unsigned char *s, size_t seg_begin, uint32_t len,
                             uint64_t start, uint64_t end, size_t i_mid, size_t i_last) {
    o->rep32[idx] = ((len - 2u) << LCP_LENGTH_SHIFT)
        | ((uint32_t)alphabet[s[seg_begin]] << (2u * LCP_SYMBOL_BITS))
        | ((uint32_t)alphabet[s[i_mid]] << LCP_SYMBOL_BITS)
        | ((uint32_t)alphabet[s[i_last]]);
    o->sz[idx] = LCP_SYMBOL_BITS * len;
    o->st[idx] = start;
    o->en[idx] = end;
}

int parse1_soa(const char *begin, const char *end, struct lcp_soa *o, uint64_t offset) {
    const size_t n = (size_t)(end - begin);
    const unsigned char *s = (const unsigned char *)begin;
    if (n >= 4) {
        const size_t nw = (n + 63) / 64;
        uint64_t dstack[512];
        uint64_t *d = nw <= 512 ? dstack : (uint64_t *)malloc(nw * sizeof(uint64_t));
        unsigned bad = 0;
#if defined(__AVX512BW__)
        {
            const __m512i cA = _mm512_set1_epi8('A'), cC = _mm512_set1_epi8('C');
            const __m512i cG = _mm512_set1_epi8('G'), cT = _mm512_set1_epi8('T');
            __mmask64 badm = 0;
            for (size_t w = 0; w < nw; w++) {
                const size_t base = w * 64;
                const size_t left = n - base;
                const __mmask64 va = left >= 64 ? ~0ull : ((1ull << left) - 1ull);
                const __mmask64 vb = left - 1 >= 64 ? ~0ull : ((1ull << (left - 1)) - 1ull);
                const __m512i a = _mm512_maskz_loadu_epi8(va, s + base);
                const __m512i b = _mm512_maskz_loadu_epi8(vb, s + base + 1);
                d[w] = _mm512_mask_cmplt_epu8_mask(vb, a, b);
                badm |= _mm512_mask_cmpeq_epi8_mask(vb, a, b);
                const __mmask64 ok = _mm512_cmpeq_epi8_mask(a, cA) | _mm512_cmpeq_epi8_mask(a, cC)
                                   | _mm512_cmpeq_epi8_mask(a, cG) | _mm512_cmpeq_epi8_mask(a, cT);
                badm |= va & ~ok;
            }
            bad = badm != 0;
        }
#else
        memset(d, 0, nw * sizeof(uint64_t));
        bad = 1u - PARSE1_ACGT[s[n - 1]];
        for (size_t j = 0; j + 1 < n; j++) {
            const unsigned a = s[j], b = s[j + 1];
            d[j >> 6] |= (uint64_t)(a < b) << (j & 63);
            bad |= (unsigned)(a == b) | (1u - PARSE1_ACGT[a]);
        }
#endif
        if (!bad) {
            size_t idx = 0, it2 = n;
            for (size_t w = 0; w < nw; w++) {
                const uint64_t D = d[w];
                const uint64_t Dn = w + 1 < nw ? d[w + 1] : 0;
                const uint64_t Dp = w ? d[w - 1] : 0;
                const uint64_t R1 = (D >> 1) | (Dn << 63);
                const uint64_t R2 = (D >> 2) | (Dn << 62);
                const uint64_t L1 = (D << 1) | (Dp >> 63);
                uint64_t minm = ~D & R1;
                uint64_t maxm = D & ~R1 & L1 & ~R2;
                const size_t base = w * 64;
                const size_t lim_min = n >= 2 && n - 2 > base ? n - 2 - base : 0;
                const size_t lim_max = n >= 3 && n - 3 > base ? n - 3 - base : 0;
                if (lim_min < 64) minm &= lim_min ? (1ull << lim_min) - 1ull : 0ull;
                if (lim_max < 64) maxm &= lim_max ? (1ull << lim_max) - 1ull : 0ull;
                if (w == 0) maxm &= ~1ull;
                uint64_t m = minm | maxm;
                while (m) {
                    const size_t i = base + (size_t)__builtin_ctzll(m);
                    m &= m - 1;
                    const size_t g = it2 < i;
                    const size_t pm = i - 1 + (i == 0);
                    soa1_core(o, idx, s, it2 - 1, (uint32_t)(i - it2 + 2), it2 - 1 + offset, i + 1 + offset, pm, i);
                    idx += g;
                    soa1_core(o, idx, s, i, 3, i + offset, i + 3 + offset, i + 1, i + 2);
                    idx++;
                    it2 = i + 3;
                }
            }
            if (d != dstack) free(d);
            return (int)idx;
        }
        if (d != dstack) free(d);
    }
    return -1;
}

void lcp_dct_soa(const struct lcp_soa *in, int n, int level, uint32_t *lab, uint32_t *lsz) {
    const uint32_t *sz = in->sz;
    if (level == 1) {
        const uint32_t *rep = in->rep32;
        for (int i = 1; i < n; i++) {
            const uint32_t o = dct1_out32(rep[i - 1], rep[i], sz[i]);
            lab[i] = o; lsz[i] = bitlen_min2_32(o);
        }
    } else {
        const uint64_t *rep = in->rep;
        for (int i = 1; i < n; i++) {
            const uint32_t o = dctu_out(rep[i - 1], sz[i - 1], rep[i], sz[i]);
            lab[i] = o; lsz[i] = bitlen_min2_clz(o);
        }
    }
}

static inline void soa3_core(struct lcp_soa *o, size_t idx, const uint32_t *LB, const uint32_t *LS,
                             const uint64_t *ST, const uint64_t *EN, size_t b, size_t distance) {
    o->st[idx] = ST[b];
    o->en[idx] = EN[b + distance - 1];
    uint32_t size = 0;
    for (size_t t = b; t < b + distance; t++) size += LS[t];
    uint64_t rep = 0;
    unsigned index = 0;
    for (size_t t = b + distance; t-- > b && index + LS[t] <= 64; ) {
        rep |= (uint64_t)LB[t] << index;
        index += LS[t];
    }
    o->rep[idx] = LCP_LEVEL1_MASK & rep;
    o->sz[idx] = minimum(size, 63u);
}

static inline void soa3_three(struct lcp_soa *o, size_t idx, const uint32_t *LB, const uint32_t *LS,
                              const uint64_t *ST, const uint64_t *EN, size_t b) {
    o->st[idx] = ST[b];
    o->en[idx] = EN[b + 2];
    const uint64_t r0 = LB[b], r1 = LB[b + 1], r2 = LB[b + 2];
    const unsigned s0 = LS[b], s1 = LS[b + 1], s2 = LS[b + 2];
    const unsigned c2 = s2 <= 64u;
    const unsigned c1 = c2 & (s2 + s1 <= 64u);
    const unsigned c0 = c1 & (s2 + s1 + s0 <= 64u);
    const uint64_t rep = (r2 & (0ull - (uint64_t)c2))
                       | ((r1 << (s2 & 63u)) & (0ull - (uint64_t)c1))
                       | ((r0 << ((s2 + s1) & 63u)) & (0ull - (uint64_t)c0));
    o->rep[idx] = LCP_LEVEL1_MASK & rep;
    o->sz[idx] = minimum(s0 + s1 + s2, 63u);
}

int parse3_soa(const uint32_t *lab, const uint32_t *lsz, const struct lcp_soa *in, int n, struct lcp_soa *o) {
    if (n < 4) return 0;
    const size_t m = (size_t)(n - 1);
    const uint32_t *LB = lab + 1;
    const uint32_t *LS = lsz + 1;
    const uint64_t *ST = in->st;
    const uint64_t *EN = in->en + 1;
    const size_t nw = (m + 63) / 64;
    uint64_t dst[512], est[512];
    uint64_t *d = nw <= 512 ? dst : (uint64_t *)malloc(nw * sizeof(uint64_t));
    uint64_t *e = nw <= 512 ? est : (uint64_t *)malloc(nw * sizeof(uint64_t));
    memset(d, 0, nw * sizeof(uint64_t));
    memset(e, 0, nw * sizeof(uint64_t));
#if defined(__AVX512F__)
    for (size_t j = 0; j + 1 < m; j += 16) {
        const size_t left = m - 1 - j;
        const __mmask16 k = left >= 16 ? (__mmask16)0xFFFF : (__mmask16)((1u << left) - 1u);
        const __m512i a = _mm512_maskz_loadu_epi32(k, LB + j);
        const __m512i b = _mm512_maskz_loadu_epi32(k, LB + j + 1);
        d[j >> 6] |= (uint64_t)_mm512_mask_cmplt_epu32_mask(k, a, b) << (j & 63);
        e[j >> 6] |= (uint64_t)_mm512_mask_cmpeq_epi32_mask(k, a, b) << (j & 63);
    }
#else
    for (size_t j = 0; j + 1 < m; j++) {
        const uint32_t a = LB[j], b = LB[j + 1];
        d[j >> 6] |= (uint64_t)(a < b) << (j & 63);
        e[j >> 6] |= (uint64_t)(a == b) << (j & 63);
    }
#endif
    size_t idx = 0, it2 = m;
    for (size_t w = 0; w < nw; w++) {
        const uint64_t D = d[w], E = e[w];
        const uint64_t Dn = w + 1 < nw ? d[w + 1] : 0, En = w + 1 < nw ? e[w + 1] : 0;
        const uint64_t Dp = w ? d[w - 1] : 0, Ep = w ? e[w - 1] : 0;
        const uint64_t R1d = (D >> 1) | (Dn << 63), R1e = (E >> 1) | (En << 63);
        const uint64_t R2d = (D >> 2) | (Dn << 62);
        const uint64_t L1d = (D << 1) | (Dp >> 63), L1e = (E << 1) | (Ep >> 63);
        const uint64_t G = ~(D | E), GR1 = ~(R1d | R1e), GL1 = ~(L1d | L1e);
        uint64_t minm = G & R1d;
        uint64_t maxm = D & GR1 & ~GL1 & ~R2d;
        uint64_t runm = ~E & R1e;
        const size_t base = w * 64;
        const size_t lim_min = m >= 2 && m - 2 > base ? m - 2 - base : 0;
        const size_t lim_max = m >= 3 && m - 3 > base ? m - 3 - base : 0;
        if (lim_min < 64) { const uint64_t k = lim_min ? (1ull << lim_min) - 1ull : 0ull; minm &= k; runm &= k; }
        if (lim_max < 64) maxm &= lim_max ? (1ull << lim_max) - 1ull : 0ull;
        if (w == 0) maxm &= ~1ull;
        uint64_t mk = minm | maxm | runm;
        while (mk) {
            const size_t i = base + (size_t)__builtin_ctzll(mk);
            const uint64_t bit = mk & (~mk + 1ull);
            mk &= mk - 1;
            size_t len = 3;
            if (runm & bit) {
                const size_t mc = 1 + ones_from(e, nw, i + 1);
                if (i + 1 + mc >= m) continue;
                len = 2 + mc;
            }
            const size_t g = it2 < i;
            soa3_core(o, idx, LB, LS, ST, EN, it2 - 1, g ? i - it2 + 2 : 1);
            idx += g;
            if (len == 3) soa3_three(o, idx, LB, LS, ST, EN, i);
            else soa3_core(o, idx, LB, LS, ST, EN, i, len);
            idx++;
            it2 = i + len;
        }
    }
    if (d != dst) free(d);
    if (e != est) free(e);
    return (int)idx;
}
