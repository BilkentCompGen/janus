#pragma once
#include <cstdint>
#include <vector>
#include "keys.h"

static const uint32_t CL = 3, DN = 32, COV = 8;

struct clus { int32_t d, bd; uint32_t n, lo, hi, a, b; uint64_t best; };

struct site { uint64_t read, key; uint32_t s, rc, fa, fb, a0, b0, span, tonly, shadow, nbr, tspan, nspan, tsame; int32_t lest, dmin, dmax; std::vector<uint32_t> la, ra, lp, rp, cov; clus cl[CL]; int32_t dn[DN]; uint32_t pn[DN]; uint8_t lk; };

void measure(std::vector<site>& w, const keys* x, const char* reads, int t);
