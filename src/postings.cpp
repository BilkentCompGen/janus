#include "postings.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {
[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[index] FATAL: %s\n", m); std::exit(1); }
}

void Index::load(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[index] FATAL: open %s\n", path.c_str()); std::exit(1); }
    char mg[4]; uint32_t ver = 0, resv = 0;
    if (std::fread(mg, 1, 4, f) != 4 || std::memcmp(mg, "SNPI", 4)) die("not an index.bin");
    if (std::fread(&ver, 4, 1, f) != 1 || ver != 1) die("index.bin version is not 1");
    if (std::fread(&S, 4, 1, f) != 1) die("read S");
    if (std::fread(&pad, 4, 1, f) != 1) die("read pad");
    if (std::fread(&K, 8, 1, f) != 1) die("read K");
    if (std::fread(&N, 4, 1, f) != 1) die("read N");
    if (std::fread(&resv, 4, 1, f) != 1) die("read header");
    off.resize(K + 1);
    if (std::fread(off.data(), 4, K + 1, f) != K + 1) die("read off");
    const uint64_t E = off[K];
    w.resize(2 * E);
    if (E && std::fread(w.data(), 8, 2 * E, f) != 2 * E) die("read postings");
    std::fclose(f);
}

Post Index::get(uint32_t node, uint32_t j) const {
    const uint64_t* p = w.data() + (size_t)(off[node] + j) * 2;
    const uint64_t w0 = p[0], w1 = p[1];
    Post o;
    o.read = (uint32_t)(w0 & 0xfffffffull);
    o.samp = (uint8_t)((w0 >> 28) & 0x3f);
    o.ord  = (uint32_t)((w0 >> 34) & 0x3fff);
    o.clen = (uint32_t)((w0 >> 48) & 0xffff);
    o.st   = (uint32_t)(w1 & 0xffffffffull);
    o.dl   = (uint32_t)((w1 >> 32) & 0xffff);
    o.dr   = (uint32_t)((w1 >> 48) & 0xffff);
    return o;
}

void SvbFile::open(const std::string& svb) {
    f = std::fopen(svb.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[index] FATAL: open %s\n", svb.c_str()); std::exit(1); }
    char hdr[32];
    if (std::fread(hdr, 1, 32, f) != 32 || std::memcmp(hdr, "SVB1", 4)) die("not a .svb");

    const std::string ip = svb + ".idx";
    FILE* g = std::fopen(ip.c_str(), "rb");
    if (!g) { std::fprintf(stderr, "[index] FATAL: %s yok, `snape svbidx %s` calistir\n",
                           ip.c_str(), svb.c_str()); std::exit(1); }
    char ih[32];
    if (std::fread(ih, 1, 32, g) != 32 || std::memcmp(ih, "SVBI", 4)) die("not a .svb.idx");
    std::memcpy(&nreads, ih + 8, 8);
    std::memcpy(&maxlen, ih + 16, 8);
    off.resize(nreads + 1);
    if (std::fread(off.data(), 8, nreads + 1, g) != nreads + 1) die("read .svb.idx");
    std::fclose(g);
}

void SvbFile::close() { if (f) std::fclose(f); f = nullptr; }

bool SvbFile::read(uint64_t no, std::string& seq) {
    if (no >= nreads) return false;
    if (std::fseek(f, (long)off[no], SEEK_SET) != 0) return false;
    uint32_t tag;
    if (std::fread(&tag, 4, 1, f) != 1) return false;
    const uint32_t L = tag & 0x7fffffffu;
    const bool raw = (tag & 0x80000000u) != 0;
    seq.resize(L);
    if (raw) return std::fread(&seq[0], 1, L, f) == L;
    const uint32_t nb = (L + 3) >> 2;
    static thread_local std::vector<uint8_t> buf;
    buf.resize(nb);
    if (nb && std::fread(buf.data(), 1, nb, f) != nb) return false;
    static const char T[4] = { 'A', 'C', 'G', 'T' };
    for (uint32_t i = 0; i < L; i++) seq[i] = T[(buf[i >> 2] >> (6 - 2 * (i & 3))) & 3];
    return true;
}

void SigLookup::build(const std::vector<uint64_t>& sig) {
    v = sig.data();
    off.assign((size_t(1) << BITS) + 1, 0);
    for (uint64_t s : sig) off[(size_t)(s >> (64 - BITS)) + 1]++;
    for (size_t i = 1; i < off.size(); i++) off[i] += off[i - 1];
}

uint32_t SigLookup::find(uint64_t s) const {
    const size_t b = (size_t)(s >> (64 - BITS));
    for (uint32_t i = off[b], e = off[b + 1]; i < e; i++) if (v[i] == s) return i;
    return UINT32_MAX;
}

std::vector<uint64_t> graph_sigs(const std::string& gp, uint32_t& S, uint64_t& K) {
    FILE* f = std::fopen(gp.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[index] FATAL: open %s\n", gp.c_str()); std::exit(1); }
    char mg[4]; uint32_t ver = 0, resv = 0; uint64_t ovs = 0, ovc = 0;
    if (std::fread(mg, 1, 4, f) != 4 || std::memcmp(mg, "SNPG", 4)) die("not a graph.bin");
    if (std::fread(&ver, 4, 1, f) != 1 || ver != 6) die("graph.bin version is not 6 -- rebuild it");
    if (std::fread(&S, 4, 1, f) != 1) die("read S");
    if (std::fread(&resv, 4, 1, f) != 1) die("read header");
    if (std::fread(&ovs, 8, 1, f) != 1 || std::fread(&ovc, 8, 1, f) != 1) die("read overlap");
    if (std::fread(&K, 8, 1, f) != 1) die("read K");
    std::vector<uint64_t> sig(K);
    if (K && std::fread(sig.data(), 8, K, f) != K) die("read sigs");
    std::fclose(f);
    return sig;
}

std::vector<uint8_t> graph_normal_only(const std::string& gp, uint32_t& S, uint64_t& K,
                                       std::vector<uint64_t>& sig) {
    FILE* f = std::fopen(gp.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[index] FATAL: open %s\n", gp.c_str()); std::exit(1); }
    char mg[4]; uint32_t ver = 0, resv = 0; uint64_t ovs = 0, ovc = 0;
    if (std::fread(mg, 1, 4, f) != 4 || std::memcmp(mg, "SNPG", 4)) die("not a graph.bin");
    if (std::fread(&ver, 4, 1, f) != 1 || ver != 6) die("graph.bin version is not 6 -- rebuild it");
    if (std::fread(&S, 4, 1, f) != 1 || S < 2) die("read S");
    if (std::fread(&resv, 4, 1, f) != 1) die("read header");
    if (std::fread(&ovs, 8, 1, f) != 1 || std::fread(&ovc, 8, 1, f) != 1) die("read overlap");
    if (std::fread(&K, 8, 1, f) != 1) die("read K");
    sig.resize(K);
    if (K && std::fread(sig.data(), 8, K, f) != K) die("read sigs");
    if (std::fseek(f, (long)(K * sizeof(cov_t)), SEEK_CUR) != 0) die("seek past normal plane");
    std::vector<uint8_t> nrm(K, 1);
    constexpr size_t CH = 1u << 20;
    std::vector<cov_t> buf(CH);
    for (uint32_t j = 1; j < S; j++)
        for (uint64_t done = 0; done < K; ) {
            const size_t m = (size_t)std::min<uint64_t>(CH, K - done);
            if (std::fread(buf.data(), sizeof(cov_t), m, f) != m) die("read tumour plane");
            for (size_t i = 0; i < m; i++) if (buf[i]) nrm[done + i] = 0;
            done += m;
        }
    std::fclose(f);
    return nrm;
}

bool core_bytes(const Index& ix, std::vector<SvbFile>& svb, uint32_t node,
                std::string& out, std::string& scratch) {
    const uint32_t n = ix.count(node);
    for (uint32_t j = 0; j < n; j++) {
        const Post p = ix.get(node, j);
        if (p.samp >= svb.size()) continue;
        if (!svb[p.samp].read(p.read, scratch)) continue;
        if ((uint64_t)p.st + p.clen > scratch.size()) continue;
        out.clear();
        char prev = 0;
        for (uint32_t i = 0; i < p.clen; i++) {
            const char c = (char)(scratch[p.st + i] & 0xdf);
            if (c != prev) out.push_back(c);
            prev = c;
        }
        return !out.empty();
    }
    return false;
}
