// sig_rc_pairing_v2.cpp — proper RC pairing test.
//
// For each sig X in counts.bin: take its stored window bytes, reverse-complement them,
// then re-run extract_fast on the RC window to obtain the LCP core (with its proper
// label) that the SAME window would yield if processed as an RC sequence. Hash that
// core's window with the LCP-derived label, then look up in counts.bin.
//
// This bypasses the chiral asymmetry of the byte hashing because we re-derive the
// label from the RC content directly.
//
// usage: sig_rc_pairing_v2 <counts.bin> <provenance.bin> [out_pairing.bin] [--level 5] [--flank 10]

#include "svp.h"
#include "MurmurHash3.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern "C" {
#include "encoding.h"
}

extern "C" void extract_fast(const char *str, int len, uint64_t offset, int target_level,
                             int *n_out, const uint64_t **starts_out,
                             const uint64_t **ends_out, const uint32_t **labels_out);

static char comp(char c) {
    switch (c & 0xDF) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
    }
    return 'N';
}

struct Prov { uint32_t n_fwd, n_rc, t_fwd, t_rc; };

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: sig_rc_pairing_v2 <counts.bin> <provenance.bin> [out.bin] [--level L] [--flank F]\n");
        return 2;
    }
    const char* cp = argv[1]; const char* pp = argv[2];
    const char* op = nullptr;
    int level = 5, flank = 10;
    int ai = 3;
    if (argc > ai && argv[ai][0] != '-') op = argv[ai++];
    for (; ai + 1 < argc; ai += 2) {
        std::string k = argv[ai]; const char* v = argv[ai + 1];
        if      (k == "--level") level = std::atoi(v);
        else if (k == "--flank") flank = std::atoi(v);
        else { std::fprintf(stderr, "[FATAL] arg %s\n", argv[ai]); std::abort(); }
    }
    LCP_INIT();

    FILE* fc = std::fopen(cp, "rb");
    if (!fc) { std::fprintf(stderr, "[FATAL] open %s\n", cp); std::abort(); }
    uint64_t K; std::fread(&K, 8, 1, fc);
    std::vector<SigCount> sigs(K);
    std::fread(sigs.data(), sizeof(SigCount), K, fc);
    uint64_t total; std::fread(&total, 8, 1, fc);
    std::vector<uint16_t> seqlen(K); std::fread(seqlen.data(), 2, K, fc);
    std::vector<uint8_t> payload(total); std::fread(payload.data(), 1, total, fc);
    std::fclose(fc);
    std::fprintf(stderr, "[v2] counts.bin: K=%lu total=%lu\n", (unsigned long)K, (unsigned long)total);

    FILE* fp = std::fopen(pp, "rb");
    if (!fp) { std::fprintf(stderr, "[FATAL] open %s\n", pp); std::abort(); }
    uint64_t Kp; std::fread(&Kp, 8, 1, fp);
    if (Kp != K) { std::fprintf(stderr, "[FATAL] K mismatch\n"); std::abort(); }
    std::vector<Prov> prov(K); std::fread(prov.data(), sizeof(Prov), K, fp);
    std::fclose(fp);

    // sig set for lookup
    std::unordered_set<uint64_t> sigset; sigset.reserve(K * 2);
    for (uint64_t i = 0; i < K; i++) sigset.insert(sigs[i].sig);
    std::fprintf(stderr, "[v2] sigset size = %zu\n", sigset.size());

    std::vector<uint64_t> seq_off(K + 1, 0);
    for (uint64_t i = 0; i < K; i++) seq_off[i + 1] = seq_off[i] + seqlen[i];

    // For each sig: RC its bytes; run extract_fast on a 'sentinel-padded' version so the
    // LCP has enough context to anchor cores; pick the central core; hash its window with
    // its label; look up in sigset.
    std::vector<uint8_t> rcbuf; rcbuf.reserve(4096);
    std::vector<uint8_t> paired(K, 0);  // 0=no, 1=yes, 2=self
    uint64_t no_core = 0, n_paired = 0, n_self = 0, n_unpaired = 0;

    const uint64_t mu = (uint64_t)flank;
    for (uint64_t i = 0; i < K; i++) {
        const uint16_t l = seqlen[i];
        rcbuf.assign(l, 0);
        const uint8_t* b = payload.data() + seq_off[i];
        for (uint16_t j = 0; j < l; j++) rcbuf[j] = (uint8_t)comp((char)b[l - 1 - j]);

        int n; const uint64_t *st, *en; const uint32_t *lb;
        extract_fast((const char*)rcbuf.data(), (int)l, 0, level, &n, &st, &en, &lb);
        if (n == 0) { no_core++; paired[i] = 0; n_unpaired++; continue; }

        // Pick the "central" core — the one whose flanking window most fills the RC bytes.
        // Heuristic: core whose center is closest to len/2.
        int best = 0;
        uint64_t best_dist = UINT64_MAX;
        const uint64_t mid = (uint64_t)l / 2;
        for (int k = 0; k < n; k++) {
            const uint64_t c = (st[k] + en[k]) / 2;
            const uint64_t d = (c > mid) ? c - mid : mid - c;
            if (d < best_dist) { best_dist = d; best = k; }
        }
        const uint64_t L = st[best] > mu ? st[best] - mu : 0;
        const uint64_t R = en[best] + mu < (uint64_t)l ? en[best] + mu : (uint64_t)l;
        const uint64_t g = murmur64((const char*)rcbuf.data() + L, (int)(R - L), lb[best]);
        if (sigset.count(g)) {
            if (g == sigs[i].sig) { paired[i] = 2; n_self++; }
            else                  { paired[i] = 1; n_paired++; }
        } else { paired[i] = 0; n_unpaired++; }

        if ((i & ((1u<<22) - 1)) == 0) {
            std::fprintf(stderr, "[v2] %lu/%lu  paired=%lu self=%lu unpaired=%lu no_core=%lu\n",
                         (unsigned long)i, (unsigned long)K,
                         (unsigned long)n_paired, (unsigned long)n_self,
                         (unsigned long)n_unpaired, (unsigned long)no_core);
        }
    }
    std::fprintf(stderr, "[v2] done: paired=%lu  palindromic-self=%lu  unpaired=%lu  no_core=%lu\n",
                 (unsigned long)n_paired, (unsigned long)n_self,
                 (unsigned long)n_unpaired, (unsigned long)no_core);

    if (op) {
        FILE* fo = std::fopen(op, "wb");
        std::fwrite(&K, 8, 1, fo); std::fwrite(paired.data(), 1, K, fo); std::fclose(fo);
        std::fprintf(stderr, "[v2] wrote %s\n", op);
    }

    auto cls = [&](uint64_t i) -> int {
        const Prov& p = prov[i];
        const bool n = (p.n_fwd + p.n_rc) > 0;
        const bool t = (p.t_fwd + p.t_rc) > 0;
        if (n && t) return 0;
        if (n)      return 1;
        if (t)      return 2;
        return 3;
    };
    const char* CN[] = {"BOTH", "N_only", "T_only", "TRUE_PHANTOM"};
    uint64_t total_by[4] = {0,0,0,0};
    uint64_t paired_by[4] = {0,0,0,0};
    uint64_t self_by[4]  = {0,0,0,0};
    for (uint64_t i = 0; i < K; i++) {
        const int c = cls(i);
        total_by[c]++;
        if      (paired[i] == 1) paired_by[c]++;
        else if (paired[i] == 2) self_by[c]++;
    }
    std::printf("\n== RC-pairing v2 by provenance class ==\n");
    std::printf("  %-14s  %-12s  %-22s  %-22s\n",
                "class", "n_sigs", "RC partner in counts", "palindromic-self");
    for (int c = 0; c < 4; c++) {
        uint64_t t = total_by[c]; if (!t) continue;
        std::printf("  %-14s  %-12lu  %lu  (%.3f%%)            %lu  (%.3f%%)\n",
                    CN[c], (unsigned long)t,
                    (unsigned long)paired_by[c], 100.0 * paired_by[c] / t,
                    (unsigned long)self_by[c],   100.0 * self_by[c]   / t);
    }
    return 0;
}
