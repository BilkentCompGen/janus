#include "svp.h"
#include <algorithm>
#include <unordered_map>

void unique_anchors(const std::vector<uint64_t>& nsig, const std::vector<uint64_t>& tsig, std::vector<Anchor>& out) {
    std::unordered_map<uint64_t, int32_t> nmap, tmap;
    nmap.reserve(nsig.size() * 2);
    tmap.reserve(tsig.size() * 2);
    for (size_t i = 0; i < nsig.size(); i++) {
        int32_t& v = nmap[nsig[i]];
        v = (v == 0) ? (int32_t)(i + 1) : -1;
    }
    for (size_t i = 0; i < tsig.size(); i++) {
        int32_t& v = tmap[tsig[i]];
        v = (v == 0) ? (int32_t)(i + 1) : -1;
    }
    out.clear();
    for (size_t i = 0; i < nsig.size(); i++) {
        const int32_t nv = nmap[nsig[i]];
        if (nv != (int32_t)(i + 1)) continue;
        const auto it = tmap.find(nsig[i]);
        if (it == tmap.end() || it->second <= 0) continue;
        out.push_back({(uint32_t)(nv - 1), (uint32_t)(it->second - 1)});
    }
}

void chain(const std::vector<Anchor>& anchors, std::vector<int>& kept) {
    const int K = (int)anchors.size();
    std::vector<int> tails, idx, parent(K, -1);
    for (int i = 0; i < K; i++) {
        const int x = (int)anchors[i].t_idx;
        const auto it = std::lower_bound(tails.begin(), tails.end(), x);
        const int p = (int)(it - tails.begin());
        parent[i] = p ? idx[p - 1] : -1;
        if (p == (int)tails.size()) { tails.push_back(x); idx.push_back(i); }
        else { tails[p] = x; idx[p] = i; }
    }
    kept.clear();
    for (int k = idx.empty() ? -1 : idx.back(); k != -1; k = parent[k]) kept.push_back(k);
    std::reverse(kept.begin(), kept.end());
}

static std::pair<uint64_t, uint64_t> ext_r(const std::string& N, uint64_t ne, uint64_t nmax, const std::string& T, uint64_t te, uint64_t tmax) {
    while (ne < nmax && te < tmax && N[ne] == T[te]) { ne++; te++; }
    return {ne, te};
}

static std::pair<uint64_t, uint64_t> ext_l(const std::string& N, uint64_t ns, uint64_t nmin, const std::string& T, uint64_t ts, uint64_t tmin) {
    while (ns > nmin && ts > tmin && N[ns - 1] == T[ts - 1]) { ns--; ts--; }
    return {ns, ts};
}

void emit_bubbles(const std::string& N, const std::vector<Core>& ncs, const std::string& T, const std::vector<Core>& tcs, const std::vector<Anchor>& anchors, const std::vector<int>& kept, std::vector<Bubble>& out) {
    out.clear();
    const int K = (int)kept.size();
    for (int i = -1; i < K; i++) {
        const uint64_t ne_start = (i == -1) ? 0u : ncs[anchors[kept[i]].n_idx].end;
        const uint64_t te_start = (i == -1) ? 0u : tcs[anchors[kept[i]].t_idx].end;
        const uint64_t ne_end = (i + 1 == K) ? N.size() : ncs[anchors[kept[i + 1]].n_idx].start;
        const uint64_t te_end = (i + 1 == K) ? T.size() : tcs[anchors[kept[i + 1]].t_idx].start;
        auto [ne, te] = ext_r(N, ne_start, ne_end, T, te_start, te_end);
        auto [ns, ts] = ext_l(N, ne_end, ne, T, te_end, te);
        const uint64_t nl = ns - ne, tl = ts - te;
        const char* type = (nl == 1 && tl == 1) ? "SNP" : (nl == 0 && tl > 0) ? "INS" : (nl > 0 && tl == 0) ? "DEL" : "CPLX";
        if (nl > 0 || tl > 0) out.push_back({ne, ns, te, ts, std::string_view(N.data() + ne, nl), std::string_view(T.data() + te, tl), type});
    }
}
