// SoA adapter over biolcp (init_lps_offset / lps_deepen). Each core's signature is MurmurHash64A
// over the core's own bytes, computed once here and returned in sigs_out. end is EXCLUSIVE (biolcp
// cores carry inclusive end, +1). LCP_INIT() runs once on first use (process-global alphabet).
extern "C" {
#include "lps.h"
#include "encoding.h"
}
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

// MurmurHash64A (Austin Appleby, public domain). Node identity = this hash over the core's raw
// bytes. One definition, used everywhere a core signature is needed (graph build, truth)
// so all sigs stay consistent; seed fixed at 0 for determinism across runs/tools.
uint64_t core_sig(const char* key, int len) {
    const uint64_t m = 0xc6a4a7935bd1e995ULL;
    const int r = 47;
    uint64_t h = (uint64_t)len * m;                        // seed 0 -> h = 0 ^ (len*m)
    const int nblocks = len / 8;
    for (int i = 0; i < nblocks; i++) {
        uint64_t k; std::memcpy(&k, key + i * 8, 8);       // memcpy: safe on unaligned core offsets
        k *= m; k ^= k >> r; k *= m;
        h ^= k; h *= m;
    }
    const unsigned char* tail = (const unsigned char*)(key + nblocks * 8);
    switch (len & 7) {
        case 7: h ^= (uint64_t)tail[6] << 48; [[fallthrough]];
        case 6: h ^= (uint64_t)tail[5] << 40; [[fallthrough]];
        case 5: h ^= (uint64_t)tail[4] << 32; [[fallthrough]];
        case 4: h ^= (uint64_t)tail[3] << 24; [[fallthrough]];
        case 3: h ^= (uint64_t)tail[2] << 16; [[fallthrough]];
        case 2: h ^= (uint64_t)tail[1] << 8;  [[fallthrough]];
        case 1: h ^= (uint64_t)tail[0];       h *= m;
    }
    h ^= h >> r; h *= m; h ^= h >> r;
    return h;
}

namespace {
struct Scratch { uint64_t* st; uint64_t* en; uint64_t* sig; int cap; };
thread_local Scratch sc{nullptr, nullptr, nullptr, 0};

void grow(int n) {
    if (n <= sc.cap) return;
    int c = n < 1024 ? 1024 : n;
    sc.st  = (uint64_t*)std::realloc(sc.st,  (size_t)c * 8);
    sc.en  = (uint64_t*)std::realloc(sc.en,  (size_t)c * 8);
    sc.sig = (uint64_t*)std::realloc(sc.sig, (size_t)c * 8);
    sc.cap = c;
}

// Homopolymer-compressed copy of the current sequence + the map back to raw offsets.
// hs[i] is the i-th run's base; hmap[i] is where that run STARTS in the raw string, and
// hmap[hn] == len closes the last run, so a compressed span [a,b) is raw [hmap[a], hmap[b]).
thread_local std::string hs;
thread_local std::vector<uint32_t> hmap;

bool init_alphabet() { LCP_INIT(); return true; }
}  // namespace

// Measured on this simulation's reads (scripts/homopolymer_errors.py, 20k reads vs the reference):
// 91.3% of insertion errors merely lengthen a homopolymer run, deletions are 4.6x enriched in runs
// of >= 10, and indels outnumber substitutions 38:1. Every one of those errors changes the bytes of
// the core it falls in, i.e. changes the core's SIGNATURE -- the same locus stops being the same node.
// Collapsing runs before parsing removes that whole error class: TTTT and TTTTT normalise to T.
// The cost is homopolymer-LENGTH resolution, which no SV-scale call depends on (no read error here
// reaches 50 bp, so nothing at SV scale was being carried by run lengths anyway).
// The normalisation extract() applies, one line below: run-length slips -- what long-read indel
// errors overwhelmingly are -- stop changing a locus identity. `map` has out.size()+1 entries,
// map.back() == len. File-local: extract() is the only caller.
static void hpc(const char* s, int len, std::string& out, std::vector<uint32_t>& map) {
    out.clear(); map.clear();
    out.reserve((size_t)len); map.reserve((size_t)len + 1);
    for (int i = 0; i < len; i++) {
        if (i == 0 || s[i] != s[i - 1]) { out.push_back(s[i]); map.push_back((uint32_t)i); }
    }
    map.push_back((uint32_t)len);            // sentinel: closes the final run
}

extern "C" void extract(const char* str, int len, uint64_t offset, int target_level,
                        int* n_out, const uint64_t** starts_out,
                        const uint64_t** ends_out, const uint64_t** sigs_out) {
    static bool once = init_alphabet(); (void)once;
    *n_out = 0; *starts_out = nullptr; *ends_out = nullptr; *sigs_out = nullptr;
    if (target_level < 1 || len < 3) return;

    // Normalise FIRST: parse and hash the homopolymer-compressed sequence, never the raw one.
    // Doing it here rather than at each call site is deliberate -- extract() is the single door
    // every sequence enters through (both graph passes, gt_cores, ref_copynum), so no caller can
    // forget to normalise and end up with signatures that disagree with everyone else's.
    hpc(str, len, hs, hmap);
    const int hn = (int)hs.size();
    if (hn < 3) return;

    struct lps p;
    init_lps_offset(&p, hs.data(), hn, 0);      // parse in compressed space...
    if (target_level > 1) lps_deepen(&p, target_level);
    const int n = p.size;
    if (n <= 0) { free_lps(&p); return; }

    grow(n);
    for (int i = 0; i < n; i++) {
        const uint64_t s = p.cores[i].start, e = p.cores[i].end + 1;
        sc.sig[i] = core_sig(hs.data() + s, (int)(e - s));
        // ...but report RAW coordinates, so every caller's positions, margins and sizes still
        // mean what they meant before. Only the identity of a core changed, not where it is.
        sc.st[i] = offset + hmap[s];
        sc.en[i] = offset + hmap[e];
    }
    free_lps(&p);
    *n_out = n; *starts_out = sc.st; *ends_out = sc.en; *sigs_out = sc.sig;
}
