#pragma once
#include <stdint.h>
#include "lps.h"

#ifdef __cplusplus
extern "C" {
#endif

// one LCP level of a read as parallel arrays: 32-bit (level 1) or 64-bit representations, positions and bit sizes
struct lcp_soa {
    uint64_t *rep;
    uint32_t *rep32;
    uint64_t *st;
    uint64_t *en;
    uint32_t *sz;
};

// level-1 cores of [begin, end); -1 if the sequence has a non-ACGT base or a homopolymer run
int parse1_soa(const char *begin, const char *end, struct lcp_soa *o, uint64_t offset);
// the labels and label sizes of the cores of one level, by deterministic coin tossing
void lcp_dct_soa(const struct lcp_soa *in, int n, int level, uint32_t *lab, uint32_t *lsz);
// the cores of the next level from the labels
int parse3_soa(const uint32_t *lab, const uint32_t *lsz, const struct lcp_soa *in, int n, struct lcp_soa *o);

#ifdef __cplusplus
}
#endif
