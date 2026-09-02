extern "C" {
#include "lps.h"
#include "encoding.h"

int parse1(const char* begin, const char* end, struct core* cores, uint64_t offset);
int parse3(struct core* begin, struct core* end, struct core* cores);
int lcp_dct(struct lps* lps_ptr);
}
#include "snape.h"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

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

struct Scratch {
    uint64_t *st = nullptr, *en = nullptr, *sig = nullptr; int cap = 0;
    struct core* cores = nullptr; int ccap = 0;
};
thread_local Scratch sc, sc2;

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

thread_local std::string hs;
thread_local std::vector<uint32_t> hmap;

bool init_alphabet() { LCP_INIT(); return true; }

static void hpc(const char* s, int len, std::string& out, std::vector<uint32_t>& map) {
    out.resize((size_t)len); map.resize((size_t)len + 1);
    size_t k = 0;
    char prev = 0;
    for (int i = 0; i < len; i++) {
        const char c = (char)(s[i] & 0xdf);
        out[k] = c; map[k] = (uint32_t)i;
        k += (c != prev);
        prev = c;
    }
    map[k] = (uint32_t)len;
    out.resize(k); map.resize(k + 1);
}

uint8_t RC_LUT[256];
bool init_rc_lut() {
    for (int i = 0; i < 256; i++) RC_LUT[i] = 'N';
    RC_LUT[(uint8_t)'A'] = 'T'; RC_LUT[(uint8_t)'T'] = 'A';
    RC_LUT[(uint8_t)'C'] = 'G'; RC_LUT[(uint8_t)'G'] = 'C';
    return true;
}
const bool RC_LUT_READY = init_rc_lut();

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

int parse_to_level(Scratch& s, const char* data, int n, int target) {
    struct core* cores = grow_cores(s, (int)(n / LCP_CONSTANT_FACTOR) + 8);
    int size = parse1(data, data + n, cores, 0);
    for (int lv = 1; lv < target && size; lv++) {
        struct lps p; p.level = lv; p.size = size; p.cores = cores;
        if (lcp_dct(&p) < 0) return 0;
        size = parse3(cores + LCP_DCT_ITERATION_COUNT, cores + size, cores);
    }
    return size;
}

void emit(Scratch& s, const struct core* cores, int n, const std::string& h,
          const std::vector<uint32_t>& map, uint64_t offset) {
    grow_out(s, n);
    for (int i = 0; i < n; i++) {
        const uint64_t a = cores[i].start, b = cores[i].end + 1;
        s.sig[i] = core_sig(h.data() + a, (int)(b - a));
        s.st[i] = offset + map[a];
        s.en[i] = offset + map[b < map.size() ? b : map.size() - 1];
    }
}

}

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
    } else {
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
