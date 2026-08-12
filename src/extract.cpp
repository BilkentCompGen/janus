// Candidate replacement for src/extract.cpp. Same cores, same signatures, same coordinates.
//
// Two costs removed, and one trap avoided:
//   1. init_lps_offset() malloc'd len/1.5 cores per call, lps_deepen1() realloc'd (shrinking) once
//      per level and free_lps() freed -- six allocator round trips per strand per read. The core
//      buffer is now thread-local and grows monotonically; the loop below is lps_deepen1()'s body
//      with the shrink dropped, so lcptools itself is untouched.
//   2. extract2() parses BOTH strands of a read from one homopolymer compression: run for run,
//      hpc(rc(s)) is rc(hpc(s)), so the reverse strand is a reversal of the compressed string and
//      its run map instead of a second compression of a second raw copy.
//
// The trap: a core's `end` is EXCLUSIVE, and emit() reads [start, end+1) -- one byte past the core,
// which for the LAST core of a read is one byte past the compressed string and one entry past the
// run map. Both reads land on whatever the shipped code left there (a std::string's terminator, and
// a vector slot that was written by an earlier read), and the graph depends on it: `en` decides
// whether that core clears the margin. So hs/hmap keep the shipped types and the shipped
// clear/reserve/push_back sequence, and the reverse strand rewrites those same buffers in place
// with the same values at the same indices, which leaves the same history behind them.
extern "C" {
#include "lps.h"
#include "encoding.h"
// Defined in lps.c, used by lps_deepen1(), but not declared in lps.h. Declared here so the
// deepening chain can run on a caller-owned buffer.
int lcp_dct(struct lps* lps_ptr);
}
#include "snape.h"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// MurmurHash64A (Austin Appleby, public domain). Node identity = this hash over the core's bytes.
// One definition, used everywhere a core signature is needed, seed fixed at 0 for determinism.
uint64_t core_sig(const char* key, int len) {
    const uint64_t m = 0xc6a4a7935bd1e995ULL;
    const int r = 47;
    uint64_t h = (uint64_t)len * m;
    const int nblocks = len / 8;
    for (int i = 0; i < nblocks; i++) {
        uint64_t k; std::memcpy(&k, key + i * 8, 8);
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

// Per-thread scratch: the SoA output and the lps core buffer, both grown monotonically.
struct Scratch {
    uint64_t *st = nullptr, *en = nullptr, *sig = nullptr; int cap = 0;
    struct core* cores = nullptr; int ccap = 0;
};
thread_local Scratch sc, sc2;      // sc: forward strand (and every extract() call). sc2: reverse.

void grow_out(Scratch& s, int n) {
    if (n <= s.cap) return;
    const int c = n < 1024 ? 1024 : n;
    s.st  = (uint64_t*)std::realloc(s.st,  (size_t)c * 8);
    s.en  = (uint64_t*)std::realloc(s.en,  (size_t)c * 8);
    s.sig = (uint64_t*)std::realloc(s.sig, (size_t)c * 8);
    s.cap = c;
}
struct core* grow_cores(Scratch& s, int n) {
    if (n > s.ccap) {
        const int c = n < 4096 ? 4096 : n;
        s.cores = (struct core*)std::realloc(s.cores, (size_t)c * sizeof(struct core));
        s.ccap = c;
    }
    return s.cores;
}

// Homopolymer-compressed copy of the current sequence + the map back to raw offsets.
// hs[i] is the i-th run's base; hmap[i] is where that run STARTS in the raw string, and
// hmap[hn] == len closes the last run, so a compressed span [a,b) is raw [hmap[a], hmap[b]).
thread_local std::string hs;
thread_local std::vector<uint32_t> hmap;

bool init_alphabet() { LCP_INIT(); return true; }

// Measured on this simulation's reads (scripts/homopolymer_errors.py, 20k reads vs the reference):
// 91.3% of insertion errors merely lengthen a homopolymer run, deletions are 4.6x enriched in runs
// of >= 10, and indels outnumber substitutions 38:1. Every one of those errors changes the bytes of
// the core it falls in, i.e. changes the core's SIGNATURE -- the same locus stops being the same node.
// Collapsing runs before parsing removes that whole error class: TTTT and TTTTT normalise to T.
static void hpc(const char* s, int len, std::string& out, std::vector<uint32_t>& map) {
    out.clear(); map.clear();
    out.reserve((size_t)len); map.reserve((size_t)len + 1);
    for (int i = 0; i < len; i++) {
        if (i == 0 || s[i] != s[i - 1]) { out.push_back(s[i]); map.push_back((uint32_t)i); }
    }
    map.push_back((uint32_t)len);            // sentinel: closes the final run
}

// ASCII complement, 'N' for everything else -- the same table graph.cpp's rc_into() applies to the
// raw read, applied here to the compressed one.
uint8_t RC_LUT[256];
bool init_rc_lut() {
    for (int i = 0; i < 256; i++) RC_LUT[i] = 'N';
    RC_LUT[(uint8_t)'A'] = 'T'; RC_LUT[(uint8_t)'T'] = 'A';
    RC_LUT[(uint8_t)'C'] = 'G'; RC_LUT[(uint8_t)'G'] = 'C';
    return true;
}
const bool RC_LUT_READY = init_rc_lut();

// hs/hmap in place: the compressed reverse complement of what they hold. Run i of the reverse is
// run (hn-1-i) of the forward complemented, and its raw span mirrors: map'[i] = len - map[hn-i].
// The ends agree with the forward convention by construction -- map'[0] = len - map[hn] = 0 and
// map'[hn] = len - map[0] = len.
void rc_in_place(std::string& h, std::vector<uint32_t>& map, int len) {
    const int hn = (int)h.size();
    for (int i = 0, j = hn - 1; i <= j; i++, j--) {
        const char a = (char)RC_LUT[(uint8_t)h[i]], b = (char)RC_LUT[(uint8_t)h[j]];
        h[i] = b; h[j] = a;
    }
    for (int i = 0, j = hn; i <= j; i++, j--) {
        const uint32_t a = map[i], b = map[j];
        map[i] = (uint32_t)len - b; map[j] = (uint32_t)len - a;
    }
}

// True when reversing in compressed space is exactly compressing the reversed raw read. It is not
// when two DIFFERENT bytes outside ACGT sit in adjacent runs: both complement to 'N', so the raw
// reverse would compress them into one run while the compressed reverse keeps two. Rare enough to
// be worth a scan and a fallback rather than a wrong answer.
bool rc_safe(const std::string& h) {
    const int hn = (int)h.size();
    bool prev_bad = false;
    for (int i = 0; i < hn; i++) {
        const char c = h[i];
        const bool bad = (c != 'A' && c != 'C' && c != 'G' && c != 'T');
        if (bad && prev_bad) return false;
        prev_bad = bad;
    }
    return true;
}

// parse1 + the deepening chain on a caller-owned buffer: lps_deepen1()'s body without the shrinking
// realloc, which only costs an allocator round trip on a buffer that is about to be reused anyway.
int parse_to_level(Scratch& s, const char* data, int n, int target) {
    struct core* cores = grow_cores(s, (int)(n / CONSTANT_FACTOR) + 8);
    int size = parse1(data, data + n, cores, 0);
    for (int lv = 1; lv < target && size; lv++) {
        struct lps p; p.level = lv; p.size = size; p.cores = cores;
        if (lcp_dct(&p) < 0) return 0;
        size = parse3(cores + DCT_ITERATION_COUNT, cores + size, cores);
    }
    return size;
}

// cores -> (raw start, raw end, signature over the compressed bytes). end is exclusive; the span
// hashed is [start, end+1), which is what the shipped extract() hashes.
void emit(Scratch& s, const struct core* cores, int n, const std::string& h,
          const std::vector<uint32_t>& map, uint64_t offset) {
    grow_out(s, n);
    for (int i = 0; i < n; i++) {
        const uint64_t a = cores[i].start, b = cores[i].end + 1;
        s.sig[i] = core_sig(h.data() + a, (int)(b - a));
        s.st[i] = offset + map[a];
        s.en[i] = offset + map[b];
    }
}

}  // namespace

extern "C" void extract(const char* str, int len, uint64_t offset, int target_level,
                        int* n_out, const uint64_t** starts_out,
                        const uint64_t** ends_out, const uint64_t** sigs_out) {
    static bool once = init_alphabet(); (void)once;
    *n_out = 0; *starts_out = nullptr; *ends_out = nullptr; *sigs_out = nullptr;
    if (target_level < 1 || len < 3) return;
    hpc(str, len, hs, hmap);
    const int hn = (int)hs.size();
    if (hn < 3) return;
    const int n = parse_to_level(sc, hs.data(), hn, target_level);
    if (n <= 0) return;
    emit(sc, sc.cores, n, hs, hmap, offset);
    *n_out = n; *starts_out = sc.st; *ends_out = sc.en; *sigs_out = sc.sig;
}

// Both strands of one read, one compression. Reverse-strand coordinates are positions in the
// reverse complement of the read -- exactly what a caller that reverse-complemented the read itself
// and called extract() on it would have got.
extern "C" void extract2(const char* str, int len, int target_level,
                         int* nf, const uint64_t** stf, const uint64_t** enf, const uint64_t** sgf,
                         int* nr, const uint64_t** str_, const uint64_t** enr, const uint64_t** sgr) {
    static bool once = init_alphabet(); (void)once;
    *nf = *nr = 0;
    *stf = *enf = *sgf = *str_ = *enr = *sgr = nullptr;
    if (target_level < 1 || len < 3) return;

    hpc(str, len, hs, hmap);
    const int hn = (int)hs.size();
    if (hn < 3) return;
    const int n1 = parse_to_level(sc, hs.data(), hn, target_level);
    if (n1 > 0) { emit(sc, sc.cores, n1, hs, hmap, 0); *nf = n1; *stf = sc.st; *enf = sc.en; *sgf = sc.sig; }

    if (rc_safe(hs)) {
        rc_in_place(hs, hmap, len);
    } else {                                   // adjacent non-ACGT runs: compress the raw reverse
        static thread_local std::string rcbuf;
        rcbuf.resize((size_t)len);
        for (int i = 0; i < len; i++) rcbuf[i] = (char)RC_LUT[(uint8_t)str[len - 1 - i]];
        hpc(rcbuf.data(), len, hs, hmap);
    }
    const int hn2 = (int)hs.size();
    if (hn2 < 3) return;
    const int n2 = parse_to_level(sc2, hs.data(), hn2, target_level);
    if (n2 > 0) { emit(sc2, sc2.cores, n2, hs, hmap, 0); *nr = n2; *str_ = sc2.st; *enr = sc2.en; *sgr = sc2.sig; }
}
