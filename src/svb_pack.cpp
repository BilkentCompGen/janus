#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "svb.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static inline int code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
    }
    return -1;
}

int cmd_pack(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: snape pack out.svb < reads.fastq\n"); return 1; }
    FILE* out = std::fopen(argv[1], "wb");
    if (!out) { std::perror("fopen out"); return 1; }

    char hdr[16] = {0};
    std::fwrite(hdr, 1, 16, out);

    uint64_t nreads = 0, total = 0;
    char* line = nullptr; size_t cap = 0; ssize_t got; int ln = 0;
    std::vector<uint8_t> pb;
    while ((got = getline(&line, &cap, stdin)) != -1) {
        if (ln == 1) {
            size_t L = (size_t)got;
            while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) L--;
            bool raw = false;
            for (size_t i = 0; i < L; i++) if (code(line[i]) < 0) { raw = true; break; }
            uint32_t tag = (uint32_t)L | (raw ? 0x80000000u : 0u);
            std::fwrite(&tag, 4, 1, out);
            if (raw) {
                std::fwrite(line, 1, L, out);
            } else {
                uint32_t nb = (uint32_t)((L + 3) >> 2);
                pb.assign(nb, 0);
                for (size_t i = 0; i < L; i++) pb[i >> 2] |= (uint8_t)(code(line[i]) << (6 - 2 * (i & 3)));
                std::fwrite(pb.data(), 1, nb, out);
            }
            nreads++; total += L;
        }
        ln = (ln + 1) & 3;
    }
    std::free(line);

    std::memcpy(hdr, &nreads, 8);
    std::memcpy(hdr + 8, &total, 8);
    std::fseek(out, 0, SEEK_SET);
    std::fwrite(hdr, 1, 16, out);
    std::fclose(out);
    std::fprintf(stderr, "[svb_pack] %s: %llu reads, %llu bases\n",
                 argv[1], (unsigned long long)nreads, (unsigned long long)total);
    return 0;
}
