#include <cstring>
#include <immintrin.h>
#include "core.h"
#include "io.h"

const Bases BASES;

static uint8_t code(int c) {
    c &= 0xdf;
    return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3;
}

Bases::Bases() {
    for (int c = 0; c < 256; c++) {
        enc[c] = code(c);
        bad[c] = code(c) == 3 && (c & 0xdf) != 'T';
    }
}

bool Bases::plain(const uint8_t* s, uint32_t n) const {
    uint32_t odd = 0;
    for (uint32_t i = 0; i < n; i++) odd |= bad[s[i]];
    return !odd;
}

void Bases::pack(const uint8_t* s, uint32_t n, uint8_t* d) const {
    uint32_t i = 0, k = 0;
    for (; i + 4 <= n; i += 4) d[k++] = enc[s[i]] << 6 | enc[s[i + 1]] << 4 | enc[s[i + 2]] << 2 | enc[s[i + 3]];
    if (i < n) {
        uint8_t c = 0;
        for (int j = 0; i < n; i++, j++) c |= enc[s[i]] << (6 - 2 * j);
        d[k] = c;
    }
}

static uint64_t rev2(uint64_t x) {
    x = __builtin_bswap64(x);
    x = ((x >> 4) & 0x0f0f0f0f0f0f0f0full) | ((x & 0x0f0f0f0f0f0f0f0full) << 4);
    return ((x >> 2) & 0x3333333333333333ull) | ((x & 0x3333333333333333ull) << 2);
}

void Read::squeeze(const uint8_t* rec) {
    uint32_t tag;
    std::memcpy(&tag, rec, 4);
    uint32_t size = tag & 0x7fffffff, nw = (size + 31) / 32, runs = 0;
    const uint8_t* q = rec + 4;
    if (tag >> 31) {
        raw.assign(nw * 8 + 8, 0);
        BASES.pack(q, size, raw.data());
        q = raw.data();
    }
    if (mk.size() < nw + 2) { rs.resize(nw + 2); mk.resize(nw + 2); cum.resize(nw + 2); word.resize(nw + 2); }
    hs.assign(nw + 2, 0);
    uint64_t prev = ~(uint64_t)q[0] >> 6 & 3;
    for (uint32_t j = 0; j < nw; j++) {
        uint64_t x = be64(q + 8 * j), d = x ^ (x >> 2 | prev << 62), s = (d | d >> 1) & 0x5555555555555555ull;
        if (size - 32 * j < 32) s &= ~0ull << (64 - 2 * (size - 32 * j));
        prev = x & 3;
        mk[j] = s;
        cum[j] = runs;
        uint32_t c = __builtin_popcountll(s), o = 2 * (runs & 31);
        uint64_t v = c ? _pext_u64(x, s | s << 1) << (64 - 2 * c) : 0;
        hs[runs >> 5] |= v >> o;
        if (o) hs[(runs >> 5) + 1] |= v << (64 - o);
        runs += c;
    }
    mk[nw] = 0;
    cum[nw] = runs;
    uint32_t kw = (runs + 31) / 32;
    for (uint32_t j = 0; j < kw; j++) rs[j] = rev2(~hs[kw - 1 - j]);
    rs[kw] = 0;
    pad = 32 * kw - runs;
    len = size;
    k = runs;
}

void Read::letters(bool reverse, char* out) const {
    const uint64_t* w = reverse ? rs.data() : hs.data();
    for (uint32_t i = reverse ? pad : 0, e = i + k; i < e; i++) *out++ = "ACGT"[w[i >> 5] >> (62 - 2 * (i & 31)) & 3];
}

void Read::starts(uint32_t* map) const {
    uint32_t q = 0;
    for (uint32_t j = 0; j < (len + 31) / 32; j++)
        for (uint64_t m = mk[j]; m; m ^= 1ull << (63 - __builtin_clzll(m))) map[q++] = 32 * j + __builtin_clzll(m) / 2;
    map[k] = len;
}

std::string Read::bases(const uint8_t* rec, uint32_t lo, uint32_t hi, uint32_t rc) const {
    const uint8_t* q = rec[3] >> 7 ? raw.data() : rec + 4;
    std::string s(hi - lo, 0);
    for (uint32_t i = lo; i < hi; i++) {
        uint32_t p = rc ? len - 1 - i : i;
        s[i - lo] = "ACGT"[(q[p >> 2] >> (6 - 2 * (p & 3)) & 3) ^ rc];
    }
    return s;
}

__extension__ typedef unsigned __int128 u128;

static uint64_t mix(uint64_t a, uint64_t b) {
    u128 p = (u128)a * b;
    return (uint64_t)p ^ (uint64_t)(p >> 64);
}

uint64_t core_hash(const uint64_t* w, uint32_t k) {
    uint64_t h = mix(k ^ 0xa0761d6478bd642full, 0xe7037ed1a0b428dbull);
    for (uint32_t i = 0; i < (k + 31) / 32; i++) h = mix(h ^ w[i], 0x9e3779b97f4a7c15ull);
    return mix(h ^ k, 0xd6e8feb86659fd93ull);
}

void Read::load(const uint8_t* exc, int w) {
    uint32_t m = (1u << w) - 1, prev = 0;
    BitReader in(exc, 48);
    n = exc[0] << 8 | exc[1];
    nf = exc[2] << 8 | exc[3];
    if (a.size() < n) { a.resize(2 * n); b.resize(2 * n); }
    for (uint32_t i = 0; i < n; i++) {
        uint32_t d = in.get(w), l;
        if (d == m) { d = in.get(32); l = in.get(32); }
        else l = in.get(w);
        if (i == nf) prev = 0;
        a[i] = prev + d;
        b[i] = a[i] + l;
        prev = a[i];
    }
}

void Read::store(std::vector<char>& out, int w) const {
    BitWriter bw(&out);
    uint32_t m = (1u << w) - 1;
    size_t t0 = out.size();
    out.resize(t0 + 6);
    for (uint32_t i = 0, prev = 0; i < n; prev = a[i++]) {
        if (i == nf) prev = 0;
        uint32_t d = a[i] - prev, l = b[i] - a[i];
        bool e = d >= m || l > m;
        if (e) bw.put(m, w);
        bw.put(d, e ? 32 : w);
        bw.put(l, e ? 32 : w);
    }
    bw.flush();
    uint8_t* u = (uint8_t*)out.data() + t0;
    size_t ub = out.size() - t0;
    u[0] = n >> 8; u[1] = n; u[2] = nf >> 8; u[3] = nf; u[4] = ub >> 8; u[5] = ub;
}

const uint8_t* Read::next(const Block& bl, size_t* o, size_t* e, int w) {
    const uint8_t* rec = (const uint8_t*)bl.in[0] + *o;
    squeeze(rec);
    load((const uint8_t*)bl.in[1] + *e, w);
    *o += svb_size(bl.in[0] + *o, bl.inlen[0] - *o);
    *e += exc_size(bl.in[1] + *e, bl.inlen[1] - *e);
    return rec;
}

uint32_t Read::rank(uint32_t x) const {
    uint32_t i = x & 31;
    return cum[x >> 5] + __builtin_popcountll(mk[x >> 5] & (i ? ~0ull << (64 - 2 * i) : 0));
}

uint32_t Read::cut(uint32_t i, uint64_t* out) const {
    uint32_t st, en;
    const uint64_t* w;
    if (i < nf) { st = rank(a[i]); en = rank(b[i]); w = hs.data(); }
    else { st = k - rank(len - a[i]) + pad; en = k - rank(len - b[i]) + pad; w = rs.data(); }
    uint32_t l = en - st, c = (l + 31) / 32, q = st >> 5, s = 2 * (st & 31);
    for (uint32_t j = 0; j < c; j++) out[j] = s ? w[q + j] << s | w[q + j + 1] >> (64 - s) : w[q + j];
    if (l & 31) out[c - 1] &= ~0ull << (64 - 2 * (l & 31));
    return l;
}

uint64_t Read::hash(uint32_t i) {
    return core_hash(word.data(), cut(i, word.data()));
}
