#include "placement.h"
#include <cmath>

namespace {
inline int64_t iabs64(int64_t x) { return x < 0 ? -x : x; }
// expected reference pos of core i given an anchor: fwd strand ref advances with read, rc retreats.
inline int64_t expect_pos(uint64_t ap, uint8_t astrand, int32_t ao, int32_t io) {
    return astrand == 0 ? (int64_t)ap + (io - ao) : (int64_t)ap - (io - ao);
}
}  // namespace

void place(const std::vector<Span>& eff, const posmap::Map& pm, std::vector<uint64_t>& locus) {
    const size_t n = eff.size();
    locus.assign(n, posmap::NOVEL);
    if (n == 0) return;

    std::vector<const posmap::Rec*> cand(n);
    std::vector<uint32_t> cnt(n);
    for (size_t i = 0; i < n; i++) cand[i] = pm.find(eff[i].sig, cnt[i]);

    // anchors (exactly one reference locus) -> pinned
    for (size_t i = 0; i < n; i++)
        if (cnt[i] == 1) locus[i] = posmap::pack_locus(cand[i][0].contig, cand[i][0].pos, cand[i][0].strand);

    // nearest-anchor arrays (O(n)): nl[i]=nearest anchor at index < i, nr[i]=nearest at index > i
    std::vector<size_t> nl(n, n), nr(n, n);
    { size_t last = n; for (size_t i = 0; i < n; i++) { nl[i] = last; if (cnt[i] == 1) last = i; } }
    { size_t last = n; for (size_t i = n; i-- > 0;) { nr[i] = last; if (cnt[i] == 1) last = i; } }

    // repeats (>=2 loci) -> the candidate colinear with the nearest anchor
    for (size_t i = 0; i < n; i++) {
        if (cnt[i] < 2) continue;
        size_t a = n;
        const size_t L = nl[i], R = nr[i];
        if (L != n && R != n) a = (i - L) <= (R - i) ? L : R;
        else if (L != n) a = L;
        else if (R != n) a = R;
        if (a == n) { locus[i] = posmap::UNRESOLVED; continue; }   // no anchor in this read -> do not guess

        const posmap::Rec& A = cand[a][0];
        const int64_t exp = expect_pos(A.pos, A.strand, (int32_t)eff[a].st, (int32_t)eff[i].st);
        const int64_t band = (int64_t)(0.15 * iabs64((int64_t)eff[i].st - (int64_t)eff[a].st) + 50);
        uint64_t bestL = UINT64_MAX; int64_t bestres = -1;
        for (uint32_t k = 0; k < cnt[i]; k++) {
            const posmap::Rec& c = cand[i][k];
            if (c.contig != A.contig || c.strand != A.strand) continue;
            const int64_t res = iabs64((int64_t)c.pos - exp);
            if (res <= band && (bestres < 0 || res < bestres)) {
                bestres = res; bestL = posmap::pack_locus(c.contig, c.pos, c.strand);
            }
        }
        // Split ONLY on a confident, colinear match. Guessing a copy scatters the repeat's edges and is
        // what tangled the graph; an unresolved copy collapses to one node per sig instead.
        locus[i] = (bestres >= 0) ? bestL : posmap::UNRESOLVED;
    }
}
