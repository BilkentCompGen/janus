#include "postings.h"

#include <stdlib.h>
#include <cstring>

void Index::load(const std::string& path) {
    FILE* f = xopen(path.c_str(), "rb");
    if (std::fread(&pad, 4, 1, f) != 1 || std::fread(&N, 4, 1, f) != 1 || std::fread(&K, 8, 1, f) != 1)
        die("index", "bad index.bin %s", path.c_str());
    off.resize(K + 1);
    if (std::fread(off.data(), 4, K + 1, f) != K + 1) die("index", "read off from %s", path.c_str());
    const uint64_t E = off[K];
    w.assign(E * 2, 0);
    if (E && std::fread(w.data(), 8, E * 2, f) != E * 2) die("index", "read postings from %s", path.c_str());
    std::fclose(f);
}

uint32_t Index::count(uint32_t node) const {
    return off[node + 1] - off[node];
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
    f = xopen(svb.c_str(), "rb");
    if (std::fseek(f, 16, SEEK_SET) != 0) die("index", "skip .svb header %s", svb.c_str());

    const std::string ip = svb + ".idx";
    FILE* g = std::fopen(ip.c_str(), "rb");
    if (!g) die("index", "%s missing, run: snape svbidx %s", ip.c_str(), svb.c_str());
    char ih[16];
    if (std::fread(ih, 1, 16, g) != 16) die("index", "not a .svb.idx %s", ip.c_str());
    std::memcpy(&nreads, ih, 8);
    off.resize(nreads + 1);
    if (std::fread(off.data(), 8, nreads + 1, g) != nreads + 1) die("index", "read %s", ip.c_str());
    std::fclose(g);
}

void SvbFile::close() { if (f) std::fclose(f); f = nullptr; }

bool SvbFile::read(uint64_t no, std::string& seq) {
    if (no >= nreads) return false;
    if (std::fseek(f, (long)off[no], SEEK_SET) != 0) return false;
    return svb_read_seq(f, seq);
}

static const uint32_t* quad_table(void) {
    static uint32_t q[256];
    static bool ready = false;
    if (!ready) {
        static const char L[4] = { 'A', 'C', 'G', 'T' };
        for (int b = 0; b < 256; b++) {
            const char c[4] = { L[(b >> 6) & 3], L[(b >> 4) & 3], L[(b >> 2) & 3], L[b & 3] };
            std::memcpy(&q[b], c, 4);
        }
        ready = true;
    }
    return q;
}

bool svb_read_seq(FILE* f, std::string& seq) {
    uint32_t tag;
    if (std::fread(&tag, 4, 1, f) != 1) return false;
    const uint32_t L = tag & 0x7fffffffu;
    seq.resize(L);
    if (tag >> 31) {
        if (L && std::fread(&seq[0], 1, L, f) != L) die("svb", ".svb raw truncated");
        return true;
    }
    const uint32_t nb = (L + 3) >> 2;
    static thread_local std::vector<uint8_t> pb;
    pb.resize(nb);
    if (nb && std::fread(pb.data(), 1, nb, f) != nb) die("svb", ".svb packed truncated");
    static const uint32_t* QUAD = quad_table();
    static const char T[4] = { 'A', 'C', 'G', 'T' };
    const uint32_t full = L >> 2;
    char* out = &seq[0];
    for (uint32_t i = 0; i < full; i++) std::memcpy(out + 4 * i, &QUAD[pb[i]], 4);
    for (uint32_t i = full << 2; i < L; i++) seq[i] = T[(pb[i >> 2] >> (6 - 2 * (i & 3))) & 3];
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

std::vector<uint64_t> graph_sigs(const std::string& gp, uint64_t& K) {
    FILE* f = xopen(gp.c_str(), "rb");
    uint32_t S;
    if (std::fread(&S, 4, 1, f) != 1) die("index", "bad graph.bin %s", gp.c_str());
    if (std::fread(&K, 8, 1, f) != 1) die("index", "read K from %s", gp.c_str());
    std::vector<uint64_t> sig(K);
    if (K && std::fread(sig.data(), 8, K, f) != K) die("index", "read sigs from %s", gp.c_str());
    std::fclose(f);
    return sig;
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
