// Repeat-core blacklist, shared by the graph builder (build) and the anchor stage (reference parse).
// A blacklisted sig is a core occurring > minfreq times in the reference (a repeat), built by
// tools/build_blacklist. When active, a maximal run of consecutive blacklisted cores is absorbed into
// its non-blacklisted flanks as two OVERLAPPING merged cores -- the left half [leftflank .. run-end)
// and the right half [run-start .. rightflank) -- which share the run's bytes (a real overlap edge,
// not a zero-overlap bridge). The bare repeat core never becomes a node.
//
// CRITICAL: the graph and the anchor MUST run the SAME merge, or the reference's merged cores will not
// match the graph's merged nodes and anchoring silently fails. That is why this lives in one place.
#include "svp.h"
#include "unordered_dense.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
ankerl::unordered_dense::set<uint64_t> g_bl;
bool g_on = false;
uint64_t g_dropped = 0;
inline bool blk(uint64_t sig) { return g_on && g_bl.count(sig); }
}

void load_blacklist(const char* path) {
    FILE* f = std::fopen(path, "rb");
    if (!f) { std::fprintf(stderr, "[FATAL] SVP_BLACKLIST open %s\n", path); std::abort(); }
    uint64_t n = 0;
    if (std::fread(&n, 8, 1, f) != 1) { std::fprintf(stderr, "[FATAL] SVP_BLACKLIST header\n"); std::abort(); }
    std::vector<uint64_t> buf(n);
    if (n && std::fread(buf.data(), 8, n, f) != n) { std::fprintf(stderr, "[FATAL] SVP_BLACKLIST body\n"); std::abort(); }
    std::fclose(f);
    g_bl.reserve(n * 2);
    for (uint64_t s : buf) g_bl.insert(s);
    g_on = true;
    std::fprintf(stderr, "[blacklist] %lu repeat-core sigs from %s ((A+X)->(X+B) merge)\n", (unsigned long)n, path);
}

bool blacklist_on() { return g_on; }
uint64_t blacklist_dropped() { return g_dropped; }

// (st,en,sig) are extract's cores; keep[] are the margin-surviving indices in read order; out gets the
// effective (merged) spans. A merged span longer than len_cap is dropped (the chain breaks there).
void blacklist_merge(const uint64_t* st, const uint64_t* en, const uint64_t* sig,
                     const std::vector<int>& keep, uint64_t len_cap,
                     std::vector<std::pair<uint32_t, uint32_t>>& out) {
    out.clear();
    const int nk = (int)keep.size();
    int i = 0;
    while (i < nk) {
        const int ci = keep[i];
        if (!blk(sig[ci])) {
            if (i + 1 < nk && blk(sig[keep[i + 1]])) {              // ci = left flank of a blacklisted run
                int j = i + 1; while (j + 1 < nk && blk(sig[keep[j + 1]])) j++;
                const uint64_t Lst = st[keep[i]], runst = st[keep[i + 1]], runen = en[keep[j]];
                const bool hasR = (j + 1 < nk);
                const uint64_t Ren = hasR ? en[keep[j + 1]] : runen;
                if (runen - Lst <= len_cap) out.push_back({(uint32_t)Lst, (uint32_t)runen}); else g_dropped++;   // left half (A+X)
                if (Ren - runst <= len_cap) out.push_back({(uint32_t)runst, (uint32_t)Ren}); else g_dropped++;   // right half (X+B)
                i = hasR ? j + 2 : j + 1;                          // consume run (+ right flank)
            } else { out.push_back({(uint32_t)st[ci], (uint32_t)en[ci]}); i++; }   // plain non-repeat core
        } else {                                                   // run with no left flank (opens inside the repeat)
            int j = i; while (j + 1 < nk && blk(sig[keep[j + 1]])) j++;
            const uint64_t runst = st[keep[i]], runen = en[keep[j]];
            const bool hasR = (j + 1 < nk);
            const uint64_t Ren = hasR ? en[keep[j + 1]] : runen;
            if (Ren - runst <= len_cap) out.push_back({(uint32_t)runst, (uint32_t)Ren}); else g_dropped++;
            i = hasR ? j + 2 : j + 1;
        }
    }
}
