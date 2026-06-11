// SoA adapter: exposes the extract_fast(...) interface that graph.cpp + cores.cpp
// rely on, implemented on top of biolcp's public init_lps_offset / lps_deepen.
//
// Two contracts the rest of svp depends on, preserved here:
//   * end is EXCLUSIVE (one-past-last char). biolcp cores carry inclusive end, so +1.
//   * the per-core identity returned in labels_out is biolcp's native ulabel (a
//     deterministic, offset-invariant content fingerprint) — the signature replacing
//     the old flank-window murmur. start/end still delimit the core bytes so callers
//     can build their ±flank storage/position window.
//
// The alphabet (A/C/G/T coefficients) is a process-global biolcp needs before any
// parse; we LCP_INIT() it exactly once on first use so every caller is safe.
extern "C" {
#include "lps.h"
#include "encoding.h"
}
#include <cstdint>
#include <cstdlib>

namespace {
struct Scratch { uint64_t* st; uint64_t* en; uint32_t* lb; int cap; };
thread_local Scratch sc{nullptr, nullptr, nullptr, 0};

void grow(int n) {
    if (n <= sc.cap) return;
    int c = n < 1024 ? 1024 : n;
    sc.st = (uint64_t*)std::realloc(sc.st, (size_t)c * 8);
    sc.en = (uint64_t*)std::realloc(sc.en, (size_t)c * 8);
    sc.lb = (uint32_t*)std::realloc(sc.lb, (size_t)c * 4);
    sc.cap = c;
}

bool init_alphabet() { LCP_INIT(); return true; }
}  // namespace

extern "C" void extract_fast(const char* str, int len, uint64_t offset, int target_level,
                             int* n_out, const uint64_t** starts_out,
                             const uint64_t** ends_out, const uint32_t** labels_out) {
    static bool once = init_alphabet(); (void)once;   // thread-safe one-time alphabet init
    *n_out = 0; *starts_out = nullptr; *ends_out = nullptr; *labels_out = nullptr;
    if (target_level < 1 || len < 3) return;

    struct lps p;
    init_lps_offset(&p, str, len, offset);
    if (target_level > 1) lps_deepen(&p, target_level);
    const int n = p.size;
    if (n <= 0) { free_lps(&p); return; }

    grow(n);
    for (int i = 0; i < n; i++) {
        sc.st[i] = (uint64_t)p.cores[i].start;
        sc.en[i] = (uint64_t)p.cores[i].end + 1;       // inclusive -> exclusive
        sc.lb[i] = p.cores[i].label;
    }
    free_lps(&p);

    *n_out = n; *starts_out = sc.st; *ends_out = sc.en; *labels_out = sc.lb;
}
