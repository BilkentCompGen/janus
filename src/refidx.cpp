// The reference's two answers, from one index and one lookup per core. See refidx.h for WHY the strand
// vote and the repeat merge belong in the same file: they are the same question asked of the same core.
//
// This replaces src/blacklist.cpp, which knew only about repeats and left the graph builder unable to
// tell a read's strand -- so it inserted both.
//
// The graph builder and the anchor/fasta tools ALL come through here, which is the whole point: the
// parse, the strand choice and the merge must be byte-identical in the build and in the reference scan,
// or the merged reference cores never match the merged graph nodes and anchoring silently fails.
#include "refidx.h"
#include "svp.h"
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
refidx::Index g_ix;
bool g_on = false;
bool g_vote = true;           // SVP_NOVOTE=1: keep the merge, drop the vote (dual-strand A/B)
uint32_t g_pct = 80;          // SVP_VOTE_PCT: a strand needs this % of the DECIDING cores

// SVP_VOTE_ERR=<p>: flip a DECIDED vote with probability p. This is not a model of anything -- it is a
// deliberate FAULT injected into the real vote, so the lab can ask what a given strand-error rate does
// to a topology instead of waiting for a genome to hand us one. Production leaves it at 0 (the rate
// measured there is 1 in 320,000).
uint32_t g_verr = 0;          // parts per million
std::atomic<uint64_t> g_dropped{0}, g_fwd{0}, g_rev{0}, g_und{0}, g_flipped{0};

// Which reads get corrupted is a hash of the READ's own core flags -- not a shared RNG. Two reasons:
// the workers are threads (a shared RNG would need a lock and would not be reproducible), and the SAME
// read must be corrupted the SAME way in the count pass and the edge pass. If it were not, the edges
// would be counted against one strand and the nodes against the other, and the graph would be garbage.
inline bool corrupt(const std::vector<uint8_t>& fl) {
  if (!g_verr) return false;
  uint64_t h = 0x9E3779B97F4A7C15ULL ^ (uint64_t)fl.size();
  for (size_t i = 0; i < fl.size(); i++)
    h = (h ^ ((uint64_t)(fl[i] + 1) * (i + 1))) * 0xBF58476D1CE4E5B9ULL;
  h ^= h >> 31;
  return (uint32_t)(h % 1000000u) < g_verr;
}

inline bool rep(const std::vector<uint8_t>& fl, int i) { return (fl[i] & refidx::REPEAT) != 0; }
}  // namespace

void ref_load(const char* path, int level) {
    if (!g_ix.open(path, level)) { std::fprintf(stderr, "[FATAL] refidx %s\n", path); std::abort(); }
    g_on = true;
    if (std::getenv("SVP_NOVOTE")) g_vote = false;
    if (const char* p = std::getenv("SVP_VOTE_PCT")) g_pct = (uint32_t)std::atoi(p);
    if (const char* p = std::getenv("SVP_VOTE_ERR")) g_verr = (uint32_t)(std::atof(p) * 1e6);
    std::fprintf(stderr,
        "[refidx] %s: %lu cores, %.0f MB | FW %lu  RC %lu  BOTH %lu (abstain)  REPEAT %lu (merged) | vote %s @%u%%\n",
        path, (unsigned long)g_ix.hdr.nentry, g_ix.len / 1e6, (unsigned long)g_ix.hdr.nfw,
        (unsigned long)g_ix.hdr.nrc, (unsigned long)g_ix.hdr.nboth, (unsigned long)g_ix.hdr.nrepeat,
        g_vote ? "on" : "OFF (dual-strand)", g_pct);
    if (g_verr)
        std::fprintf(stderr, "[refidx] FAULT INJECTED: %.4f%% of decided votes are flipped (SVP_VOTE_ERR)\n",
                     g_verr / 1e4);
}

bool     ref_on()      { return g_on; }
uint64_t ref_dropped() { return g_dropped.load(std::memory_order_relaxed); }

void ref_read_stats(uint64_t& fwd, uint64_t& rev, uint64_t& undecided) {
    fwd = g_fwd.load(std::memory_order_relaxed);
    rev = g_rev.load(std::memory_order_relaxed);
    undecided = g_und.load(std::memory_order_relaxed);
}
uint64_t ref_flipped() { return g_flipped.load(std::memory_order_relaxed); }

void ref_flags(const uint64_t* sig, const std::vector<int>& keep, std::vector<uint8_t>& fl) {
    fl.clear();
    fl.reserve(keep.size());
    for (int i : keep) fl.push_back(g_on ? g_ix.flags(sig[i]) : 0);
}

// Only FW_ONLY and RC_ONLY cores get a ballot. BOTH sits in an inverted repeat (same bytes forwards and
// backwards -- no strand information), NEITHER is a read error or an SV junction (not reference at all).
// Both abstain, so the denominator is the DECIDING cores, not the read's cores.
int ref_vote(const std::vector<uint8_t>& fl, bool tally) {
    if (!g_on || !g_vote) return 0;
    uint32_t fw = 0, rc = 0;
    for (uint8_t f : fl) {
        const uint8_t c = f & 3;
        fw += (c == refidx::FW_ONLY);
        rc += (c == refidx::RC_ONLY);
    }
    const uint32_t tot = fw + rc;
    auto hit = [&](std::atomic<uint64_t>& c, int v) { if (tally) c.fetch_add(1, std::memory_order_relaxed); return v; };
    const bool bad = tot && corrupt(fl);              // SVP_VOTE_ERR: a decided vote, deliberately flipped
    if (bad && tally) g_flipped.fetch_add(1, std::memory_order_relaxed);
    if (tot && fw * 100 >= tot * g_pct) return hit(g_fwd, bad ? -1 : +1);
    if (tot && rc * 100 >= tot * g_pct) return hit(g_rev, bad ? +1 : -1);
    return hit(g_und, 0);                            // a split vote, or no reference support at all
}

void ref_merge(const char* data, const uint64_t* st, const uint64_t* en, const uint64_t* sig,
               const std::vector<int>& keep, const std::vector<uint8_t>& fl, uint64_t len_cap,
               std::vector<Span>& out) {
    out.clear();
    const int nk = (int)keep.size();
    // No index: every kept core passes through untouched, and extract's own sig already names it.
    if (!g_on) { for (int i = 0; i < nk; i++) out.push_back({(uint32_t)st[keep[i]], (uint32_t)en[keep[i]], sig[keep[i]]}); return; }

    // A merged span's bytes are not any single core's bytes, so its name has to be hashed afresh.
    auto merged = [&](uint64_t a, uint64_t b) {
        if (b - a > len_cap) { g_dropped.fetch_add(1, std::memory_order_relaxed); return; }
        out.push_back({(uint32_t)a, (uint32_t)b, core_sig(data + a, (int)(b - a))});
    };
    int i = 0;
    while (i < nk) {
        if (!rep(fl, i)) {
            if (i + 1 >= nk || !rep(fl, i + 1)) {    // a plain core: reuse extract's sig, no re-hash
                out.push_back({(uint32_t)st[keep[i]], (uint32_t)en[keep[i]], sig[keep[i]]});
                i++;
                continue;
            }
            int j = i + 1;                            // keep[i] is the left flank A of a repeat run X
            while (j + 1 < nk && rep(fl, j + 1)) j++;
            const bool hasR = (j + 1 < nk);
            merged(st[keep[i]], en[keep[j]]);                                        // (A+X)
            merged(st[keep[i + 1]], hasR ? en[keep[j + 1]] : en[keep[j]]);           // (X+B)
            i = hasR ? j + 2 : j + 1;
        } else {                                      // the read opens inside the run: no left flank
            int j = i;
            while (j + 1 < nk && rep(fl, j + 1)) j++;
            const bool hasR = (j + 1 < nk);
            merged(st[keep[i]], hasR ? en[keep[j + 1]] : en[keep[j]]);
            i = hasR ? j + 2 : j + 1;
        }
    }
}
