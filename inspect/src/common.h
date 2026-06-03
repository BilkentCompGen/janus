#pragma once
// Shared loaders for inspect/ binaries. Mirrors the binary file formats
// written by src/graph.cpp (lcp-cortex-merged HEAD b7e0663).
//
// counts.bin   u64 N + N * { u64 sig, u32 n, u32 t }                       16 B/rec
// edges.bin    u64 M + M * { u64 sig_a, u64 sig_b }                        16 B/rec  (sig-level)
// graph.bin    u64 K + u8 color[K] + u32 ncov[K] + u32 tcov[K] + u64 E + (u32,u32)[E]
// posmap.bin   u64 nc + nc * { u32 namelen, char name[], u64 chrom_len } + u64 ns + ns * { u64 sig, u32 np, np * Pos }
//                where Pos = { u16 chrom_idx, u16 _pad, u32 pos_start, u32 pos_end, u8 strand, u8 _, u16 _pad }   16 B/rec
// seqs.bin     TBD via inspect/src/seqs_probe.cpp (Stage 0)
// census.bin   TBD via inspect/src/bubble_census.cpp (Stage 1)

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

namespace insp {

[[noreturn]] inline void die(const char* msg) {
    std::fprintf(stderr, "[FATAL] %s\n", msg);
    std::abort();
}

inline FILE* xopen(const std::string& p, const char* m) {
    FILE* f = std::fopen(p.c_str(), m);
    if (!f) { std::fprintf(stderr, "[FATAL] fopen %s\n", p.c_str()); std::abort(); }
    return f;
}

inline void xread(FILE* f, void* p, size_t n) {
    if (n && std::fread(p, 1, n, f) != n) {
        std::fprintf(stderr, "[FATAL] truncated read (%zu B)\n", n);
        std::abort();
    }
}

// ---------- counts.bin ----------
struct SigCount {
    uint64_t sig;
    uint32_t n;
    uint32_t t;
};

inline void load_counts(const std::string& path, std::vector<SigCount>& out) {
    FILE* f = xopen(path, "rb");
    uint64_t N; xread(f, &N, 8);
    out.resize(N);
    xread(f, out.data(), N * sizeof(SigCount));
    std::fclose(f);
}

// ---------- graph.bin ----------
struct Graph {
    uint64_t K = 0;
    std::vector<uint8_t> color;      // K
    std::vector<uint32_t> ncov;      // K
    std::vector<uint32_t> tcov;      // K
    std::vector<std::pair<uint32_t, uint32_t>> edges;  // E
};

inline void load_graph(const std::string& path, Graph& g) {
    FILE* f = xopen(path, "rb");
    xread(f, &g.K, 8);
    g.color.resize(g.K); xread(f, g.color.data(), g.K);
    g.ncov.resize(g.K);  xread(f, g.ncov.data(),  g.K * 4);
    g.tcov.resize(g.K);  xread(f, g.tcov.data(),  g.K * 4);
    uint64_t E; xread(f, &E, 8);
    g.edges.resize(E);
    xread(f, g.edges.data(), E * sizeof(g.edges[0]));
    std::fclose(f);
}

// ---------- posmap.bin ----------
struct Pos {
    uint16_t chrom_idx;
    uint16_t _pad1;
    uint32_t pos_start;
    uint32_t pos_end;
    uint8_t  strand;
    uint8_t  _pad2;
    uint16_t _pad3;
};
static_assert(sizeof(Pos) == 16, "Pos must be 16 B");

struct Posmap {
    std::vector<std::string> chrom_names;
    std::vector<uint64_t> chrom_lens;
    std::unordered_map<uint64_t, std::vector<Pos>> sig_pos;
};

inline void load_posmap(const std::string& path, Posmap& pm) {
    FILE* f = xopen(path, "rb");
    uint64_t nc; xread(f, &nc, 8);
    pm.chrom_names.resize(nc);
    pm.chrom_lens.resize(nc);
    for (uint64_t i = 0; i < nc; i++) {
        uint32_t L; xread(f, &L, 4);
        pm.chrom_names[i].resize(L);
        xread(f, pm.chrom_names[i].data(), L);
        xread(f, &pm.chrom_lens[i], 8);
    }
    uint64_t ns; xread(f, &ns, 8);
    pm.sig_pos.reserve(ns);
    for (uint64_t i = 0; i < ns; i++) {
        uint64_t sig; xread(f, &sig, 8);
        uint32_t np; xread(f, &np, 4);
        std::vector<Pos> v(np);
        if (np) xread(f, v.data(), np * sizeof(Pos));
        pm.sig_pos.emplace(sig, std::move(v));
    }
    std::fclose(f);
}

// ---------- seqs.bin ----------
// Schema (confirmed by inspect/src/seqs_probe.cpp on 2026-06-03):
//   u64 K                            // node count (== graph.bin K)
//   u64 total_payload                // = sum(seqlen[0..K))
//   u16 seqlen[K]                    // per-node window length
//   u8  payload[total_payload]       // concatenated ACGT sequences, no separators
//
// Random access: prefix-sum seqlens into offset[K+1] at load.
struct Seqs {
    uint64_t K = 0;
    uint64_t total_payload = 0;
    std::vector<uint16_t> seqlen;     // K entries
    std::vector<uint64_t> offset;     // K+1 entries (prefix sum)
    std::vector<uint8_t>  payload;    // total_payload bytes

    // Random-access getter. Returns pointer + length.
    const uint8_t* seq(uint64_t i, uint16_t& len_out) const {
        len_out = seqlen[i];
        return payload.data() + offset[i];
    }
};

inline void load_seqs(const std::string& path, Seqs& s) {
    FILE* f = xopen(path, "rb");
    xread(f, &s.K, 8);
    xread(f, &s.total_payload, 8);
    s.seqlen.resize(s.K);
    xread(f, s.seqlen.data(), s.K * 2);
    s.payload.resize(s.total_payload);
    xread(f, s.payload.data(), s.total_payload);
    std::fclose(f);
    s.offset.resize(s.K + 1);
    s.offset[0] = 0;
    for (uint64_t i = 0; i < s.K; i++) s.offset[i + 1] = s.offset[i] + s.seqlen[i];
    if (s.offset[s.K] != s.total_payload) die("seqs.bin: prefix-sum != declared total");
}

}  // namespace insp
