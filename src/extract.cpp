// SoA adapter over biolcp (init_lps_offset / lps_deepen). Each core's signature is MurmurHash64A
// over the core's own bytes, computed once here and returned in sigs_out. end is EXCLUSIVE (biolcp
// cores carry inclusive end, +1). LCP_INIT() runs once on first use (process-global alphabet).
extern "C" {
#include "lps.h"
#include "encoding.h"
}
#include <cstdint>
#include <cstdlib>
#include <cstring>

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

bool init_alphabet() { LCP_INIT(); return true; }
}  // namespace

extern "C" void extract(const char* str, int len, uint64_t offset, int target_level,
                        int* n_out, const uint64_t** starts_out,
                        const uint64_t** ends_out, const uint64_t** sigs_out) {
    static bool once = init_alphabet(); (void)once;
    *n_out = 0; *starts_out = nullptr; *ends_out = nullptr; *sigs_out = nullptr;
    if (target_level < 1 || len < 3) return;

    struct lps p;
    init_lps_offset(&p, str, len, offset);
    if (target_level > 1) lps_deepen(&p, target_level);
    const int n = p.size;
    if (n <= 0) { free_lps(&p); return; }

    grow(n);
    for (int i = 0; i < n; i++) {
        const uint64_t s = p.cores[i].start, e = p.cores[i].end + 1;
        sc.st[i] = s; sc.en[i] = e;
        sc.sig[i] = core_sig(str + (s - offset), (int)(e - s));
    }
    free_lps(&p);
    *n_out = n; *starts_out = sc.st; *ends_out = sc.en; *sigs_out = sc.sig;
}
