#include "lps.h"
#include <stdio.h>
#include <stdlib.h>
int lcp_dct(struct lps *lps_ptr);
int main(void) {
    enum { N = 1 << 22 };
    struct core *c = malloc(N * sizeof *c);
    uint64_t *rep = malloc(N * 8), *st = malloc(N * 8), *en = malloc(N * 8);
    uint32_t *rep32 = malloc(N * 4), *sz = malloc(N * 4), *lab = malloc(N * 4), *lsz = malloc(N * 4);
    srand(7);
    long bad1 = 0, badu = 0, ub = 0, cmp1 = 0;
    for (int i = 0; i < N; i++) {
        const uint64_t f = rand() & 3, m = rand() & 3, l = rand() & 3, cnt = rand() % 6, len = cnt + 2;
        c[i].bit_rep = LCP_LEVEL1_TAG | (cnt << LCP_LENGTH_SHIFT) | (f << 4) | (m << 2) | l;
        c[i].bit_size = (ubit_size)(2 * len) + (rand() % 3 == 0 ? (rand() % 80) : 0);
        c[i].start = i; c[i].end = i + len;
        rep32[i] = (uint32_t)(c[i].bit_rep & LCP_LEVEL1_MASK); sz[i] = c[i].bit_size; st[i] = i; en[i] = i + len;
    }
    struct lcp_soa in = { rep, rep32, st, en, sz };
    for (int i = 1; i < N; i++) {
        const uint32_t L = rep32[i - 1], R = rep32[i];
        const uint32_t x = (L ^ R) & 15u, L2 = (L >> 2) & 3, L1 = (L >> 4) & 3, R2 = (R >> 2) & 3, R1 = (R >> 4) & 3;
        const uint32_t Lm = L >> 6, Rm = R >> 6;
        if (x == 0 && Lm != Rm) { const uint32_t a2 = Lm < Rm ? L1 : L2, b2 = Lm < Rm ? R2 : R1; if (a2 == b2) { ub++; rep32[i] = 0; } }
    }
    struct lps p = { 1, N, c };
    lcp_dct(&p);
    { int shown = 0; for (int i = 1; i < N && shown < 6; i++) if (rep32[i] == 0 && rep32[i - 1] != 0) { const uint32_t L = rep32[i - 1], R = (uint32_t)(c[i].bit_rep); (void)R; printf("UB pair i=%d: L=0x%x Lm=%u  lib out rep=%llu (bits %u)\n", i, L, L >> 6, (unsigned long long)c[i].bit_rep, c[i].bit_size); shown++; } }
    lcp_dct_soa(&in, N, 1, lab, lsz);
    for (int i = 1; i < N; i++) {
        if (rep32[i] == 0 || rep32[i - 1] == 0) continue;
        cmp1++;
        if (c[i].bit_rep != lab[i] || c[i].bit_size != lsz[i]) {
            if (bad1 < 5) printf("L1 mismatch i=%d lib rep=%llu sz=%u  soa rep=%u sz=%u\n", i, (unsigned long long)c[i].bit_rep, c[i].bit_size, lab[i], lsz[i]);
            bad1++;
        }
    }
    for (int i = 0; i < N; i++) {
        const uint32_t s = 2 + rand() % 62;
        const uint64_t r = (((uint64_t)rand() << 31) ^ (uint64_t)rand() ^ ((uint64_t)rand() << 47)) & ((1ull << s) - 1ull);
        c[i].bit_rep = rand() % 5 == 0 && i ? c[i - 1].bit_rep : r; c[i].bit_size = rand() % 5 == 0 && i ? c[i - 1].bit_size : s;
        rep[i] = c[i].bit_rep; sz[i] = c[i].bit_size;
    }
    p.level = 2;
    lcp_dct(&p);
    lcp_dct_soa(&in, N, 2, lab, lsz);
    for (int i = 1; i < N; i++) if (c[i].bit_rep != lab[i] || c[i].bit_size != lsz[i]) {
        if (badu < 5) printf("LU mismatch i=%d lib rep=%llu sz=%u  soa rep=%u sz=%u\n", i, (unsigned long long)c[i].bit_rep, c[i].bit_size, lab[i], lsz[i]);
        badu++;
    }
    printf("dct test: level1 compared %ld pairs (skipped %ld ctz(0) pairs), mismatches %ld; upper %d pairs, mismatches %ld\n", cmp1, ub, bad1, N - 1, badu);
    return 0;
}
