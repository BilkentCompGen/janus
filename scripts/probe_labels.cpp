#include "../src/extract.cpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
static uint64_t HB[8][65], HS[8][65], NL[8], PAIRS = 0, UB = 0, MAXM = 0;
int main(int argc, char** argv) {
    const char* path = argv[1]; const uint64_t max_bases = std::strtoull(argv[2], nullptr, 10);
    LCP_INIT();
    FILE* f = std::fopen(path, "rb"); char hdr[32]; if (std::fread(hdr, 1, 32, f) != 32) return 1;
    static const char LUT[4] = {'A','C','G','T'};
    std::string seq; std::vector<uint8_t> pb; uint64_t bases = 0;
    while (bases < max_bases) {
        uint32_t tag; if (std::fread(&tag, 4, 1, f) != 1) break;
        const uint32_t L = tag & 0x7fffffffu; seq.resize(L);
        if (tag >> 31) { if (std::fread(&seq[0], 1, L, f) != L) return 1; }
        else { const uint32_t nb = (L + 3) >> 2; pb.resize(nb); if (std::fread(pb.data(), 1, nb, f) != nb) return 1;
               for (uint32_t i = 0; i < L; i++) seq[i] = LUT[(pb[i >> 2] >> (6 - 2 * (i & 3))) & 3]; }
        bases += L; if (L < 3) continue;
        hpc(seq.data(), (int)L, hs, hmap);
        if (argc > 3) for (size_t k = 1; k < hs.size(); k++) { if (k % 397 == 0) hs[k] = hs[k - 1]; if (k % 1013 == 0) hs[k] = 'N'; }
        const int n = (int)hs.size(); if (n < 3) continue;
        struct core* cores = grow_cores(sc, (int)(n / LCP_CONSTANT_FACTOR) + 8);
        int size = parse1(hs.data(), hs.data() + n, cores, 0);
        for (int i = 1; i < size; i++) {
            const uint64_t L = cores[i - 1].bit_rep & LCP_LEVEL1_MASK, R = cores[i].bit_rep & LCP_LEVEL1_MASK;
            const uint64_t x = (L ^ R) & 15u, L2 = (L >> 2) & 3, L1 = (L >> 4) & 3, R2 = (R >> 2) & 3, R1 = (R >> 4) & 3, Lm = L >> 6, Rm = R >> 6;
            PAIRS++;
            if (x == 0 && Lm != Rm) { const uint64_t a2 = Lm < Rm ? L1 : L2, b2 = Lm < Rm ? R2 : R1; UB += (a2 == b2); }
            if (Lm > MAXM) MAXM = Lm;
        }
        for (int lv = 1; lv < 5 && size; lv++) {
            struct lps p; p.level = lv; p.size = size; p.cores = cores;
            if (lcp_dct(&p) < 0) break;
            for (int i = 1; i < size; i++) { HB[lv][bitlen_u64(cores[i].bit_rep)]++; NL[lv]++; }
            size = parse3(cores + 1, cores + size, cores);
            for (int i = 0; i < size; i++) HS[lv + 1][cores[i].bit_size]++;
        }
    }
    std::printf("bases %.3f Gbp (forward strand, hpc)  level-1 pairs %llu  ctz(0) case %llu  max mid count %llu\n", bases / 1e9, (unsigned long long)PAIRS, (unsigned long long)UB, (unsigned long long)MAXM);
    for (int lv = 1; lv < 5; lv++) {
        uint64_t n = NL[lv], le8 = 0, le16 = 0, le12 = 0; int mx = 0;
        for (int b = 0; b <= 64; b++) { if (HB[lv][b]) mx = b; if (b <= 8) le8 += HB[lv][b]; if (b <= 12) le12 += HB[lv][b]; if (b <= 16) le16 += HB[lv][b]; }
        std::printf("DCT out at level %d: labels %llu  max bits %d  <=8b %.6f  <=12b %.6f  <=16b %.6f  >16b %llu\n", lv, (unsigned long long)n, mx,
                    (double)le8 / n, (double)le12 / n, (double)le16 / n, (unsigned long long)(n - le16));
        uint64_t m = 0, cap = HS[lv + 1][63]; for (int b = 0; b <= 64; b++) m += HS[lv + 1][b];
        std::printf("  level %d cores %llu  bit_size==63 (capped) %llu (%.4f%%)\n", lv + 1, (unsigned long long)m, (unsigned long long)cap, m ? 100.0 * cap / m : 0.0);
    }
    return 0;
}
