#include "lps.h"
#if defined(__AVX512F__)
#include <immintrin.h>
#endif

#define minimum(a, b) ((a) < (b) ? (a) : (b))

#if LCP_FIXED_CORE
#include "emit_fixed.h"
#else
#include "emit_var.h"
#endif


/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
#if !LCP_FIXED_CORE
#if !defined(__AVX512BW__)
static const uint8_t PARSE1_ACGT[256] = { ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1 };
#endif

static inline void parse1_fast_core(struct core *cr, const unsigned char *s, size_t seg_begin, uint64_t len,
                                    uint64_t start, uint64_t end, size_t i_mid, size_t i_last) {
    cr->bit_rep = LCP_LEVEL1_TAG
        | ((len - 2) << LCP_LENGTH_SHIFT)
        | ((uint64_t)alphabet[s[seg_begin]] << (2u * LCP_SYMBOL_BITS))
        | ((uint64_t)alphabet[s[i_mid]] << LCP_SYMBOL_BITS)
        | ((uint64_t)alphabet[s[i_last]]);
#if LCP_LABEL_BITS != 0
    cr->label = (lcp_label)cr->bit_rep;
#endif
    cr->bit_size = (ubit_size)(LCP_SYMBOL_BITS * len);
#if LCP_POS_BITS != 0
    cr->start = (lcp_pos)start;
    cr->end = (lcp_pos)end;
#else
    (void)start; (void)end;
#endif
}

static int parse1_fast(const char *begin, const char *end, struct core *cores, uint64_t offset) {
    const size_t n = (size_t)(end - begin);
    const unsigned char *s = (const unsigned char *)begin;
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
    if (bad) { if (d != dstack) free(d); return -1; }

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
            parse1_fast_core(&cores[idx], s, it2 - 1, (uint64_t)(i - it2 + 2), it2 - 1 + offset, i + 1 + offset, pm, i);
            idx += g;
            parse1_fast_core(&cores[idx], s, i, 3, i + offset, i + 3 + offset, i + 1, i + 2);
            idx++;
            it2 = i + 3;
        }
    }
    if (d != dstack) free(d);
    return (int)idx;
}
#endif

int parse1(const char *begin, const char *end, struct core *cores, uint64_t offset) {

#if !LCP_FIXED_CORE
    if (end - begin >= 4) {
        const int r = parse1_fast(begin, end, cores, offset);
        if (r >= 0) return r;
    }
#endif

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;
    int last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == -1) {
            last_invalid_char_index = it1 - begin;
            continue;
        }

        if (alphabet[(unsigned char)*it1] == alphabet[(unsigned char)*(it1+1)]) {
            continue;
        }

        // check for RINT core
        if (alphabet[(unsigned char)*(it1+1)] == alphabet[(unsigned char)*(it1+2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 + 2;
            while (temp < end && alphabet[(unsigned char)*(temp-1)] == alphabet[(unsigned char)*temp]) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                    core_index = emit_core1(cores, core_index, begin, offset, it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                core_index = emit_core1(cores, core_index, begin, offset, it1, 2+middle_count, it1-begin+offset, it2-begin+offset);

                continue;
            }
        }

        if (alphabet[(unsigned char)*it1] > alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] < alphabet[(unsigned char)*(it1+2)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                core_index = emit_core1(cores, core_index, begin, offset, it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            alphabet[(unsigned char)*it1] < alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] > alphabet[(unsigned char)*(it1+2)] &&
            alphabet[(unsigned char)*(it1-1)] <= alphabet[(unsigned char)*(it1)] &&
            alphabet[(unsigned char)*(it1+2)] >= alphabet[(unsigned char)*(it1+3)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                core_index = emit_core1(cores, core_index, begin, offset, it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

#if !LCP_ALPHABET_PROTEIN

/**
 * @brief Parses a sequence to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores using complement alphabet.
 *
 * This function iterates over a sequence defined by iterators `begin` and `end` and identifies key
 * segments, called "cores," that represent the (LCP) regions. By analyzing
 * character relationships in the sequence (such as equality or relative order based on complement), 
 * it builds and stores these cores for further processing in the LCP framework.
 *
 * Positions are reverse-complement indices: a character at address p sits at
 * index (end - p - 1). Every core here therefore reports
 *
 *     start = end - seg_begin - 1        end = start + seg_len
 *
 * where seg_begin is the segment's highest address.
 *
 * @param begin Iterator pointing to the beginning of the sequence to parse.
 * @param end Iterator pointing to the end of the sequence to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @param offset The distance measure where the indecies of the core will be shifted by.
 * @return Size of the cores identified in the given string.
 */
int parse2(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = end - 1;
    const char *it2 = begin - 1;
    int core_index = 0;

    // find lcp cores
    for (; begin <= it1 - 2; it1--) {

        // skip invalid character
        if (rc_alphabet[(unsigned char)*it1] == rc_alphabet[(unsigned char)*(it1-1)]) {
            continue;
        }

        // check for RINT core
        if (rc_alphabet[(unsigned char)*(it1-1)] == rc_alphabet[(unsigned char)*(it1-2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 - 2;
            while (begin <= temp && rc_alphabet[(unsigned char)*(temp+1)] == rc_alphabet[(unsigned char)*temp]) {
                temp--;
                middle_count++;
            }
            if (begin <= temp) {
                // check if there is any SSEQ cores left behind
                if (it1 < it2) {
                    core_index = emit_core2(cores, core_index, end, offset, it2+1, it2-it1+2, end-it2-2+offset, end-it1+offset);
                }

                // create RINT core
                it2 = it1 - 2 - middle_count;
                core_index = emit_core2(cores, core_index, end, offset, it1, 2+middle_count, end-it1-1+offset, end-it2-1+offset);

                continue;
            }
        }

        if (rc_alphabet[(unsigned char)*it1] > rc_alphabet[(unsigned char)*(it1-1)] &&
            rc_alphabet[(unsigned char)*(it1-1)] < rc_alphabet[(unsigned char)*(it1-2)]) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                core_index = emit_core2(cores, core_index, end, offset, it2+1, it2-it1+2, end-it2-2+offset, end-it1+offset);
            }

            // create LMIN core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (begin <= it1-3 &&
            rc_alphabet[(unsigned char)*it1] < rc_alphabet[(unsigned char)*(it1-1)] &&
            rc_alphabet[(unsigned char)*(it1-1)] > rc_alphabet[(unsigned char)*(it1-2)] &&
            rc_alphabet[(unsigned char)*(it1+1)] <= rc_alphabet[(unsigned char)*(it1)] &&
            rc_alphabet[(unsigned char)*(it1-2)] >= rc_alphabet[(unsigned char)*(it1-3)]) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                core_index = emit_core2(cores, core_index, end, offset, it2+1, it2-it1+2, end-it2-2+offset, end-it1+offset);
            }

            // create LMAX core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

#endif /* !LCP_ALPHABET_PROTEIN */

/**
 * @brief Parses a array of cores to extract Locally Consisted Parsing (LCP) cores and stores them in a 
 * array of cores.
 *
 * This function iterates over a array of `core` structures defined by iterators `begin` and `end` and 
 * identifies key segments, called "cores," that represent the (LCP) regions. By analyzing
 * `core` structure relationships in the array (such as equality or relative order), it builds and stores
 * these cores for further processing in the LCP framework.
 *
 * @param begin Iterator pointing to the beginning of the `core` array to parse.
 * @param end Iterator pointing to the end of the `core` array to parse.
 * @param cores Pointer to a array where the identified LCP cores will be stored.
 * @return Size of the cores identified in the given string.
 */
#if LCP_FIXED_CORE
static int parse3_general(struct core *begin, struct core *end, struct core *cores);
#endif

#if !LCP_FIXED_CORE
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

#if LCP_LABEL_BITS == 0
static inline void core3_any(struct core *cr, const struct core *b, lcp_pos distance) {
#if LCP_POS_BITS != 0
    cr->start = b->start;
    cr->end = (b + distance - 1)->end;
#endif
    ubit_size size = 0;
    for (const struct core *it = b; it < b + distance; it++) size += it->bit_size;
    uint64_t rep = 0;
    unsigned index = 0;
    for (const struct core *it = b + distance - 1; b <= it && index + it->bit_size <= 64; it--) {
        rep |= it->bit_rep << index;
        index += it->bit_size;
    }
    cr->bit_rep = 0x7FFFFFFFFFFFFFFF & rep;
    cr->bit_size = minimum(size, 63);
}
static inline void core3_three(struct core *cr, const struct core *b) {
#if LCP_POS_BITS != 0
    cr->start = b[0].start;
    cr->end = b[2].end;
#endif
    const uint64_t r0 = b[0].bit_rep, r1 = b[1].bit_rep, r2 = b[2].bit_rep;
    const unsigned s0 = b[0].bit_size, s1 = b[1].bit_size, s2 = b[2].bit_size;
    const unsigned c2 = s2 <= 64u;
    const unsigned c1 = c2 & (s2 + s1 <= 64u);
    const unsigned c0 = c1 & (s2 + s1 + s0 <= 64u);
    const uint64_t rep = (r2 & (0ull - (uint64_t)c2))
                       | ((r1 << (s2 & 63u)) & (0ull - (uint64_t)c1))
                       | ((r0 << ((s2 + s1) & 63u)) & (0ull - (uint64_t)c0));
    cr->bit_rep = 0x7FFFFFFFFFFFFFFF & rep;
    cr->bit_size = minimum(s0 + s1 + s2, 63u);
}
#else
#define core3_any(cr, b, distance) init_core3((cr), (struct core *)(b), (distance))
#define core3_three(cr, b) init_core3((cr), (struct core *)(b), 3)
#endif

static int parse3_fast(struct core *begin, struct core *end, struct core *cores) {
    const size_t n = (size_t)(end - begin);
    if (n < 3) return 0;
    const size_t nw = (n + 63) / 64;
    uint64_t dst[512], est[512];
    uint64_t *d = nw <= 512 ? dst : (uint64_t *)malloc(nw * sizeof(uint64_t));
    uint64_t *e = nw <= 512 ? est : (uint64_t *)malloc(nw * sizeof(uint64_t));
    memset(d, 0, nw * sizeof(uint64_t));
    memset(e, 0, nw * sizeof(uint64_t));
#if defined(__AVX512F__)
    {
        const long long st = (long long)(sizeof(struct core) / 8);
        const __m512i vidx = _mm512_set_epi64(7 * st, 6 * st, 5 * st, 4 * st, 3 * st, 2 * st, st, 0);
        const __m512i z = _mm512_setzero_si512();
        for (size_t j = 0; j + 1 < n; j += 8) {
            const size_t left = n - 1 - j;
            const __mmask8 k = left >= 8 ? (__mmask8)0xFF : (__mmask8)((1u << left) - 1u);
            const __m512i a = _mm512_mask_i64gather_epi64(z, k, vidx, &begin[j].bit_rep, 8);
            const __m512i b = _mm512_mask_i64gather_epi64(z, k, vidx, &begin[j + 1].bit_rep, 8);
            d[j >> 6] |= (uint64_t)_mm512_mask_cmplt_epu64_mask(k, a, b) << (j & 63);
            e[j >> 6] |= (uint64_t)_mm512_mask_cmpeq_epi64_mask(k, a, b) << (j & 63);
        }
    }
#else
    for (size_t j = 0; j + 1 < n; j++) {
        const uint64_t a = begin[j].bit_rep, b = begin[j + 1].bit_rep;
        d[j >> 6] |= (uint64_t)(a < b) << (j & 63);
        e[j >> 6] |= (uint64_t)(a == b) << (j & 63);
    }
#endif
    size_t idx = 0, it2 = n;
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
        const size_t lim_min = n >= 2 && n - 2 > base ? n - 2 - base : 0;
        const size_t lim_max = n >= 3 && n - 3 > base ? n - 3 - base : 0;
        if (lim_min < 64) { const uint64_t k = lim_min ? (1ull << lim_min) - 1ull : 0ull; minm &= k; runm &= k; }
        if (lim_max < 64) maxm &= lim_max ? (1ull << lim_max) - 1ull : 0ull;
        if (w == 0) maxm &= ~1ull;
        uint64_t m = minm | maxm | runm;
        while (m) {
            const size_t i = base + (size_t)__builtin_ctzll(m);
            const uint64_t bit = m & (~m + 1ull);
            m &= m - 1;
            size_t len = 3;
            if (runm & bit) {
                const size_t mc = 1 + ones_from(e, nw, i + 1);
                if (i + 1 + mc >= n) continue;
                len = 2 + mc;
            }
            const size_t g = it2 < i;
            core3_any(&cores[idx], begin + (it2 - 1), (lcp_pos)(g ? i - it2 + 2 : 1));
            idx += g;
            if (len == 3) core3_three(&cores[idx], begin + i);
            else core3_any(&cores[idx], begin + i, (lcp_pos)len);
            idx++;
            it2 = i + len;
        }
    }
    if (d != dst) free(d);
    if (e != est) free(e);
    return (int)idx;
}
#endif

int parse3(struct core *begin, struct core *end, struct core *cores) {
#if !LCP_FIXED_CORE
    return parse3_fast(begin, end, cores);
#else
    return parse3_general(begin, end, cores);
#endif
}

#if LCP_FIXED_CORE
static int parse3_general(struct core *begin, struct core *end, struct core *cores) {

    struct core *it1 = begin;
    struct core *it2 = end;
    int core_index = 0;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (core_eq(it1, it1+1)) {
            continue;
        }

        // check for RINT core
        if (core_eq(it1+1, it1+2)) {

            // count middle characters
            uint32_t middle_count = 1;
            struct core *temp = it1 + 2;
            while (temp < end && core_eq(temp-1, temp)) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    core_index = emit_core3(cores, core_index, it2-1, it1-it2+2);
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                core_index = emit_core3(cores, core_index, it1, it2-it1);

                continue;
            }
        }

        // check for LMIN
        if (core_gt(it1, it1+1) && core_lt(it1+1, it1+2)) {
            
            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                core_index = emit_core3(cores, core_index, it2-1, it1-it2+2);
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            core_lt(it1, it1+1) &&
            core_gt(it1+1, it1+2) &&
            core_leq(it1-1, it1) &&
            core_geq(it1+2, it1+3)) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                core_index = emit_core3(cores, core_index, it2-1, it1-it2+2);
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }
    }
    return core_index;
}
#endif

void init_lps(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/LCP_CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, 0);
}

void init_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/LCP_CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, offset);
}

#if !LCP_ALPHABET_PROTEIN

void init_lps2(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/LCP_CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse2(str, str+len, lps_ptr->cores, 0);
}

#endif /* !LCP_ALPHABET_PROTEIN */

void init_lps3(struct lps *lps_ptr, FILE *in) {
    // read the level from the binary file
    if (fread(&(lps_ptr->level), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading level from file\n");
        exit(EXIT_FAILURE);
    }

    // read the size (number of cores)
    if(fread(&(lps_ptr->size), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading size from file\n");
        exit(EXIT_FAILURE);
    }

    lps_ptr->cores = NULL;

    if (lps_ptr->size) {
        // allocate memory for the cores array
        lps_ptr->cores = (struct core *)malloc(lps_ptr->size * sizeof(struct core));
        if (fread(lps_ptr->cores, lps_ptr->size * sizeof(struct core), 1, in) != 1) {
            fprintf(stderr, "Error reading cores from file\n");
            exit(EXIT_FAILURE);
        }
    }
}

#if LCP_POS_BITS != 0

void init_lps4(struct lps *lps_ptr, const char *str, int len, int lcp_level, int chunk_size) {

    if (lcp_level < 1)
        return;

    lps_ptr->level = 1;
    lps_ptr->size = 0; 
    int estimated_size = (int)(len / pow((double)LCP_CONSTANT_FACTOR, lcp_level));
    lps_ptr->cores = (struct core *)malloc(estimated_size*sizeof(struct core));

    int str_index = 0, core_index = 0;

    {
        int str_len = minimum(chunk_size, len);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str, str_len, 0);
        lps_deepen(&temp_lps, lcp_level);

        if (temp_lps.size) {
            memcpy(lps_ptr->cores, temp_lps.cores, (temp_lps.size)*sizeof(struct core));
            core_index = (temp_lps.size);
            lps_ptr->size = (temp_lps.size);
            if (temp_lps.size>1)
                str_index = lps_ptr->cores[core_index-2].start;
            else 
                str_index = lps_ptr->cores[core_index-1].start;
        }
        free(temp_lps.cores);
    }

    while (str_index < len) {
        int str_len = minimum(chunk_size, len-str_index);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str+str_index, str_len, str_index);
        lps_deepen(&temp_lps, lcp_level);

        if (1<temp_lps.size) {
            int overlap = 2;
            while (0<overlap) {
                if (lps_ptr->cores[core_index-overlap].start == temp_lps.cores[0].start)
                    break;
                overlap--;
            }
            memcpy(lps_ptr->cores+core_index, temp_lps.cores+overlap, (temp_lps.size-overlap)*sizeof(struct core));
            core_index += (temp_lps.size-overlap);
            lps_ptr->size += (temp_lps.size-overlap);

            if ((uint64_t)str_index < lps_ptr->cores[core_index-2].start) {
                str_index = lps_ptr->cores[core_index-2].start;
                free(temp_lps.cores);
                continue;
            } 
        }
        
        // find next start point
        for(int i=str_index+str_len-1; str_index <= i; i--) {
            if (alphabet[(unsigned char)*(str+i)] == -1) {
                str_index = i+1;
                break;
            }
        }
        if (alphabet[(unsigned char)*(str+str_index)] != -1) { // all of the characters are valid, so not valid cores found
            str_index += str_len;
        }
        
        free(temp_lps.cores);
    }

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));
}

#endif

void free_lps(struct lps *lps_ptr) {
    free(lps_ptr->cores);
    lps_ptr->size = 0;
}

void write_lps(struct lps *lps_ptr, FILE *out) {
    // write the level field
    fwrite(&(lps_ptr->level), sizeof(int), 1, out);

    // write the size (number of cores)
    fwrite(&(lps_ptr->size), sizeof(int), 1, out);

    // write each core object iteratively
    if (lps_ptr->size) {
        fwrite(lps_ptr->cores, lps_ptr->size*sizeof(struct core), 1, out);
    }
}

int64_t lps_memsize(const struct lps *lps_ptr) {
    return sizeof(struct lps) + lps_ptr->size * sizeof(struct core);
}

/**
 * @brief Performs Deterministic Coin Tossing (DCT) compression on binary sequences.
 *
 * This function is a central part of the LCP (Locally Consisted Parsing) algorithm. It identifies differences
 * between consecutive binary strings, compressing the information by focusing on the position and value of
 * the first divergent bit from the right-end of the strings. This difference is used to generate a compact
 * 'core' that encapsulates the unique elements of each sequence in a smaller binary form.
 *
 * This compression significantly reduces redundant information, making further analysis of the sequences
 * within the LCP framework more efficient and manageable.
 *
 * @return 0 if dct is performed, -1 if no enough cores are available for dct.
 */
int lcp_dct(struct lps *lps_ptr) {

    // at least 2 cores are needed for compression
    if (lps_ptr->size < LCP_DCT_ITERATION_COUNT + 1) {
        return -1;
    }

    if (lps_ptr->level == 1) {
        for (uint64_t dct_index = 0; dct_index < LCP_DCT_ITERATION_COUNT; dct_index++) {
            struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

            for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
                core_compress_level1(it_left, it_right);
            }
        }
    } else {
        for (uint64_t dct_index = 0; dct_index < LCP_DCT_ITERATION_COUNT; dct_index++) {
            struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

            for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
                core_compress_upper(it_left, it_right);
            }
        }
    }

    return 0;
}

int lps_deepen1(struct lps *lps_ptr) {

    // compress cores
    if (lcp_dct(lps_ptr) < 0) {
        lps_ptr->size = 0;
        lps_ptr->level++;
        return 0;
    }

    // find new cores
#if !LCP_FIXED_CORE
    int new_size = parse3(lps_ptr->cores + LCP_DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, lps_ptr->cores);
#else
    struct core *new_cores = (struct core *)malloc((lps_ptr->size - LCP_DCT_ITERATION_COUNT) * sizeof(struct core));
    if (new_cores == NULL) {
        return 0;
    }

    int new_size = parse3(lps_ptr->cores + LCP_DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, new_cores);

    // replace the previous level; struct core owns no heap memory
    free(lps_ptr->cores);
    lps_ptr->cores = new_cores;
#endif
    lps_ptr->size = new_size;

    lps_ptr->level++;

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));

    return 1;
}

int lps_deepen(struct lps *lps_ptr, int lcp_level) {

    if (lcp_level <= lps_ptr->level)
        return 0;

    while (lps_ptr->level < lcp_level && lps_deepen1(lps_ptr))
        ;

    return 1;
}

typedef struct {
    struct core *cores;     // destination cores (real array)
    struct core dummy_left; // uncompressed core at begin-1
    int offset_begin;       // inclusive i
    int offset_end;         // exclusive i
    int flags;              // 1 bit (process head) 1 bit (level 1 or not)
} dct_worker_args_t;

/**
 * @brief Worker routine for parallel DCT compression.
 *
 * This function is executed by a single POSIX thread and performs
 * compression on a contiguous subrange of cores. The behavior depends
 * on the `flags` field:
 *
 * - flags == -1 : No work is performed (thread exits immediately).
 * - flags % 2   : Use core_compress_level1().
 * - otherwise   : Use core_compress_upper().
 *
 * If flags > 1, the worker also performs an additional boundary
 * compression using `dummy_left` and the first core in its range.
 * This is required when the thread's assigned range does not begin
 * at the global DCT starting index, ensuring correctness across
 * chunk boundaries.
 *
 * The function assumes that:
 * - Each thread operates on a disjoint range of cores.
 * - Any required boundary state is provided through `dummy_left`.
 *
 * @param argp Pointer to a dct_worker_args_t structure containing
 *             the thread parameters.
 *
 * @return Always returns NULL (required by pthread signature).
 */
static void *dct_worker(void *argp) {
    dct_worker_args_t *dct_params = (dct_worker_args_t *)argp;

    if (dct_params->flags == -1) return NULL;

    if (dct_params->flags % 2) {
        for (int i = dct_params->offset_end - 1; dct_params->offset_begin < i; i--) {
            core_compress_level1(dct_params->cores + (i - 1), dct_params->cores + i);
        }
        if (dct_params->flags > 1) {
            core_compress_level1(&(dct_params->dummy_left), dct_params->cores + dct_params->offset_begin);
        }
    } else {
        for (int i = dct_params->offset_end - 1; dct_params->offset_begin < i; i--) {
            core_compress_upper(dct_params->cores + (i - 1), dct_params->cores + i);
        }
        if (dct_params->flags > 1) {
            core_compress_upper(&(dct_params->dummy_left), dct_params->cores + dct_params->offset_begin);
        }
    }

    return NULL;
}

/**
 * @brief Executes one parallel DCT sweep over the cores array.
 *
 * This function divides the compression work starting at `dct_index`
 * into approximately equal chunks and distributes them among
 * `thread_number` threads. Each thread processes a contiguous
 * subrange of cores using the dct_worker() routine.
 *
 * Work partitioning:
 * - total = lps_ptr->size - dct_index
 * - chunk size is computed using ceiling division.
 * - Each thread receives a non-overlapping [begin, end) interval.
 *
 * Boundary handling:
 * - Threads whose range does not start at dct_index receive a
 *   copy of the left boundary core (`dummy_left`) to ensure
 *   correct cross-boundary compression.
 *
 * Memory management:
 * - Dynamically allocates thread handles and argument arrays.
 * - Joins all threads before returning.
 *
 * @param lps_ptr       Pointer to the LPS structure containing cores.
 * @param dct_index     Starting index for this DCT iteration.
 * @param thread_number Number of worker threads to spawn.
 *
 * @return 0 on success.
 * @return -1 on allocation failure.
 */
static int run_parallel_sweep(struct lps *lps_ptr, int dct_index, int thread_number) {

    const int total = (int)lps_ptr->size - dct_index;   // number of pairs (i, i+1)
    if (total <= 0) return 0;                               // nothing to do
    const int chunk = (total + thread_number - 1) / thread_number;

    pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * thread_number);
    dct_worker_args_t *args = (dct_worker_args_t *)malloc(sizeof(dct_worker_args_t) * thread_number);

    if (!threads || !args) {
        perror("LCP: couldn't allocate thread arguments\n");
        return -1;
    }

    for (int t = 0; t < thread_number; t++) {
        int begin = dct_index + t * chunk;
        int end = ((begin + chunk) <= lps_ptr->size ? (begin + chunk) : lps_ptr->size);

        args[t].cores = lps_ptr->cores;
        args[t].offset_begin = begin;
        args[t].offset_end = end;

        if ((int)lps_ptr->size <= begin || end <= begin) {
            args[t].flags = -1;
        } else if (begin != dct_index) {
            args[t].dummy_left = lps_ptr->cores[begin-1];
            args[t].flags = 2 + (lps_ptr->level == 1 ? 1 : 0);
        } else {
            args[t].flags = (lps_ptr->level == 1 ? 1 : 0);
        }
    }

    for (int t = 0; t < thread_number; t++) {
        pthread_create(&threads[t], NULL, dct_worker, &args[t]);
    }

    for (int t = 0; t < thread_number; t++) {
        pthread_join(threads[t], NULL);
    }

    free(threads);
    free(args);
    return 0;
}

/**
 * @brief Performs multi-iteration parallel Deterministic Coin Tossing (DCT).
 *
 * This function executes LCP_DCT_ITERATION_COUNT consecutive compression
 * sweeps over the cores stored in `lps_ptr`, using parallel execution
 * for each sweep.
 *
 * Each iteration reduces redundancy by compressing adjacent core pairs.
 * The process prepares the core sequence for subsequent parsing stages
 * in the LCP framework.
 *
 * Preconditions:
 * - At least LCP_DCT_ITERATION_COUNT + 1 cores must be available.
 *
 * @param lps_ptr       Pointer to the LPS structure.
 * @param thread_number Number of worker threads per sweep.
 *
 * @return 0 on success.
 * @return -1 if there are not enough cores for DCT.
 * @return Propagates non-zero errors from run_parallel_sweep().
 */
int lcp_dct_parallel(struct lps *lps_ptr, int thread_number) {
    // at least 2 cores are needed for compression
    if (lps_ptr->size < LCP_DCT_ITERATION_COUNT + 1) {
        return -1;
    }

    for (int dct_index = 0; dct_index < LCP_DCT_ITERATION_COUNT; dct_index++) {
        int rc = run_parallel_sweep(lps_ptr, dct_index, thread_number);
        if (rc != 0) {
            return rc;
        }
    }

    return 0;
}

/**
 * @brief Performs a single LPS deepening step using parallel DCT.
 *
 * This function advances the LPS structure by one level:
 *
 * 1. Executes parallel DCT compression.
 * 2. Parses the resulting cores using parse3() to identify new cores.
 * 3. Updates the core array and shrinks it via realloc().
 * 4. Increments the LPS level.
 *
 * If DCT cannot be performed (insufficient cores),
 * the structure is reset (size = 0) and the level is still incremented.
 *
 * @param lps_ptr       Pointer to the LPS structure to update.
 * @param thread_number Number of worker threads used during DCT.
 *
 * @return 1 if deepening produced a valid new level.
 * @return 0 if no further deepening was possible.
 */
int lps_deepen1_parallel(struct lps *lps_ptr, int thread_number) {

    // compress cores
    if (lcp_dct_parallel(lps_ptr, thread_number) < 0) {
        lps_ptr->size = 0;
        lps_ptr->level++;
        return 0;
    }

    // find new cores
#if !LCP_FIXED_CORE
    int new_size = parse3(lps_ptr->cores + LCP_DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, lps_ptr->cores);
#else
    struct core *new_cores = (struct core *)malloc((lps_ptr->size - LCP_DCT_ITERATION_COUNT) * sizeof(struct core));
    if (new_cores == NULL) {
        return 0;
    }

    int new_size = parse3(lps_ptr->cores + LCP_DCT_ITERATION_COUNT, lps_ptr->cores + lps_ptr->size, new_cores);

    // replace the previous level; struct core owns no heap memory
    free(lps_ptr->cores);
    lps_ptr->cores = new_cores;
#endif
    lps_ptr->size = new_size;

    lps_ptr->level++;

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));

    return 1;
}

/**
 * @brief Deepens the LPS structure up to a target LCP level.
 *
 * Repeatedly invokes lps_deepen1_parallel() until either:
 * - The desired `lcp_level` is reached, or
 * - Further deepening is no longer possible.
 *
 * If the current level already satisfies the requested level,
 * no work is performed.
 *
 * @param lps_ptr       Pointer to the LPS structure.
 * @param lcp_level     Target LCP level to reach.
 * @param thread_number Number of worker threads used per deepening step.
 *
 * @return 1 if processing completed (even if stopped early).
 * @return 0 if no deepening was required.
 */
int lps_deepen_parallel(struct lps *lps_ptr, int lcp_level, int thread_number) {

    if (lcp_level <= lps_ptr->level)
        return 0;

    while (lps_ptr->level < lcp_level && lps_deepen1_parallel(lps_ptr, thread_number))
        ;

    return 1;
}

void print_lps(const struct lps *lps_ptr) {
    printf("Level: %d \n", lps_ptr->level);
    for(int i=0; i<lps_ptr->size; i++) {
        print_core(&(lps_ptr->cores[i]));
        printf(" ");
    }
}

int lps_eq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 0;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 0;
        }
    }

    return 1;
}

int lps_neq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 1;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 1;
        }
    }

    return 0;
}
#if !LCP_FIXED_CORE && LCP_LABEL_BITS == 0 && LCP_POS_BITS != 0 && LCP_DCT_ITERATION_COUNT == 1 && !LCP_ALPHABET_PROTEIN

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
#endif
