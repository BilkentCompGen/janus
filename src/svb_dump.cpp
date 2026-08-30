#include "svb.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

int cmd_dump(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: snape dump in.svb > reads.fastq\n"); return 1; }
    FILE* f = std::fopen(argv[1], "rb");
    if (!f) { std::perror("fopen"); return 1; }
    char hdr[32];
    if (std::fread(hdr, 1, 32, f) != 32 || std::memcmp(hdr, "SVB1", 4) != 0) {
        std::fprintf(stderr, "[FATAL] bad .svb header %s\n", argv[1]); return 1; }

    static const char LUT[4] = { 'A', 'C', 'G', 'T' };
    std::string seq, qual; std::vector<uint8_t> pb;
    uint32_t tag; uint64_t idx = 0;
    while (std::fread(&tag, 4, 1, f) == 1) {
        const uint32_t L = tag & 0x7fffffffu;
        seq.resize(L);
        if (tag >> 31) {
            if (L && std::fread(&seq[0], 1, L, f) != L) { std::fprintf(stderr, "[FATAL] raw truncated\n"); return 1; }
        } else {
            const uint32_t nb = (L + 3) >> 2;
            pb.resize(nb);
            if (nb && std::fread(pb.data(), 1, nb, f) != nb) { std::fprintf(stderr, "[FATAL] packed truncated\n"); return 1; }
            for (uint32_t i = 0; i < L; i++) seq[i] = LUT[(pb[i >> 2] >> (6 - 2 * (i & 3))) & 3];
        }
        qual.assign(L, 'I');
        std::printf("@r%llu\n%s\n+\n%s\n", (unsigned long long)idx++, seq.c_str(), qual.c_str());
    }
    std::fclose(f);
    return 0;
}
