extern "C" {
#include "lps.h"
#include "encoding.h"
}
#include "snape.h"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

static uint64_t core_sig(const char* key, int len) {
    const uint64_t m = 0xc6a4a7935bd1e995ULL;
    const int r = 47;
    uint64_t h = (uint64_t)len * m;
    const int nblocks = len / 8;
    for (int i = 0; i < nblocks; i++) {
        uint64_t k;
        std::memcpy(&k, key + i * 8, 8);
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

struct Soa { uint64_t *rep = nullptr, *st = nullptr, *en = nullptr; uint32_t *rep32 = nullptr, *sz = nullptr; int cap = 0; };

struct Scratch {
    uint64_t *st = nullptr, *en = nullptr, *sig = nullptr; int cap = 0;
    Soa a, b; uint32_t *lab = nullptr, *lsz = nullptr; int lcap = 0;
    const Soa* fin = nullptr;
};

static void grow_soa(Soa& s, int n) {
    if (n <= s.cap) return;
    const int c = n < 4096 ? 4096 : n;
    s.rep   = (uint64_t*)std::realloc(s.rep,   (size_t)c * 8);
    s.st    = (uint64_t*)std::realloc(s.st,    (size_t)c * 8);
    s.en    = (uint64_t*)std::realloc(s.en,    (size_t)c * 8);
    s.rep32 = (uint32_t*)std::realloc(s.rep32, (size_t)c * 4);
    s.sz    = (uint32_t*)std::realloc(s.sz,    (size_t)c * 4);
    s.cap = c;
}

static void grow_lab(Scratch& s, int n) {
    if (n <= s.lcap) return;
    const int c = n < 4096 ? 4096 : n;
    s.lab = (uint32_t*)std::realloc(s.lab, (size_t)c * 4);
    s.lsz = (uint32_t*)std::realloc(s.lsz, (size_t)c * 4);
    s.lcap = c;
}

static void grow_out(Scratch& s, int n) {
    if (n <= s.cap) return;
    const int c = n < 1024 ? 1024 : n;
    s.st  = (uint64_t*)std::realloc(s.st,  (size_t)c * 8);
    s.en  = (uint64_t*)std::realloc(s.en,  (size_t)c * 8);
    s.sig = (uint64_t*)std::realloc(s.sig, (size_t)c * 8);
    s.cap = c;
}

static void hpc(const char* s, int len, std::string& out, std::vector<uint32_t>& map) {
    out.resize((size_t)len);
    map.resize((size_t)len + 1);
    size_t k = 0;
    char prev = 0;
    for (int i = 0; i < len; i++) {
        const char c = (char)(s[i] & 0xdf);
        out[k] = c;
        map[k] = (uint32_t)i;
        k += (c != prev);
        prev = c;
    }
    map[k] = (uint32_t)len;
    out.resize(k);
    map.resize(k + 1);
}

static uint8_t RC_LUT[256];

static void ensure_rc_lut(void) {
    static bool ready = false;
    if (ready) return;
    for (int i = 0; i < 256; i++) RC_LUT[i] = 'N';
    RC_LUT[(uint8_t)'A'] = 'T'; RC_LUT[(uint8_t)'T'] = 'A';
    RC_LUT[(uint8_t)'C'] = 'G'; RC_LUT[(uint8_t)'G'] = 'C';
    ready = true;
}

static void ensure_lcp(void) {
    static bool ready = false;
    if (ready) return;
    LCP_INIT();
    ready = true;
}

static void rc_in_place(std::string& h, std::vector<uint32_t>& map, int len) {
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

static bool rc_safe(const std::string& h) {
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

static int parse_to_level(Scratch& s, const char* data, int n, int target) {
    const int cap = (int)(n / LCP_CONSTANT_FACTOR) + 8;
    grow_soa(s.a, cap); grow_soa(s.b, cap); grow_lab(s, cap);
    Soa* a = &s.a; Soa* b = &s.b;
    struct lcp_soa ia = { a->rep, a->rep32, a->st, a->en, a->sz };
    int size = parse1_soa(data, data + n, &ia, 0);
    if (size <= 0) return 0;
    for (int lv = 1; lv < target && size; lv++) {
        if (size < 2) return 0;
        struct lcp_soa in = { a->rep, a->rep32, a->st, a->en, a->sz },
                     out = { b->rep, b->rep32, b->st, b->en, b->sz };
        lcp_dct_soa(&in, size, lv, s.lab, s.lsz);
        size = parse3_soa(s.lab, s.lsz, &in, size, &out);
        Soa* t = a; a = b; b = t;
    }
    s.fin = a;
    return size;
}

static void emit(Scratch& s, int n, const std::string& h, const std::vector<uint32_t>& map) {
    grow_out(s, n);
    for (int i = 0; i < n; i++) {
        const uint64_t a = s.fin->st[i];
        const uint64_t b = s.fin->en[i] + 1;
        s.sig[i] = core_sig(h.data() + a, (int)(b - a));
        s.st[i] = map[a];
        s.en[i] = map[b < map.size() ? b : map.size() - 1];
    }
}

static int parse_emit(Scratch& s, const std::string& h, const std::vector<uint32_t>& map, int target) {
    const int hn = (int)h.size();
    if (hn < 3) return 0;
    const int n = parse_to_level(s, h.data(), hn, target);
    if (n <= 0) return 0;
    emit(s, n, h, map);
    return n;
}

static thread_local Scratch sc, sc2;
static thread_local std::string hs;
static thread_local std::vector<uint32_t> hmap;

static int strand(const char* str, int len, int target, Scratch& s) {
    hpc(str, len, hs, hmap);
    return parse_emit(s, hs, hmap, target);
}

extern "C" void extract(const char* str, int len, int target_level,
                         int* n_out, const uint64_t** starts_out,
                         const uint64_t** ends_out, const uint64_t** sigs_out) {
    ensure_lcp();
    ensure_rc_lut();
    *n_out = 0; *starts_out = nullptr; *ends_out = nullptr; *sigs_out = nullptr;
    if (target_level < 1 || len < 3) return;
    const int n = strand(str, len, target_level, sc);
    if (n > 0) { *n_out = n; *starts_out = sc.st; *ends_out = sc.en; *sigs_out = sc.sig; }
}

extern "C" void extract2(const char* str, int len, int target_level,
                          int* nf, const uint64_t** stf, const uint64_t** enf, const uint64_t** sgf,
                          int* nr, const uint64_t** str_, const uint64_t** enr, const uint64_t** sgr) {
    ensure_lcp();
    ensure_rc_lut();
    *nf = *nr = 0;
    *stf = *enf = *sgf = *str_ = *enr = *sgr = nullptr;
    if (target_level < 1 || len < 3) return;

    const int n1 = strand(str, len, target_level, sc);
    if (n1 > 0) { *nf = n1; *stf = sc.st; *enf = sc.en; *sgf = sc.sig; }

    if (rc_safe(hs)) {
        rc_in_place(hs, hmap, len);
    } else {
        static thread_local std::string rcbuf;
        rcbuf.resize((size_t)len);
        for (int i = 0; i < len; i++) rcbuf[i] = (char)RC_LUT[(uint8_t)str[len - 1 - i]];
        hpc(rcbuf.data(), len, hs, hmap);
    }
    const int n2 = parse_emit(sc2, hs, hmap, target_level);
    if (n2 > 0) { *nr = n2; *str_ = sc2.st; *enr = sc2.en; *sgr = sc2.sig; }
}
