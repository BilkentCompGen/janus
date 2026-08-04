// svp pack: FASTQ (stdin) -> one .svb file (argv[1]). 2-bit packs pure-ACGT reads
// (A=00 C=01 G=10 T=11, base i in byte i/4, high 2 bits first); RAW-escapes any read
// containing a non-ACGT base (e.g. N) so packing is lossless. svp reads .svb natively.
// Streaming, O(1) memory; only the sequence line of each FASTQ record is kept.
// Format matches src/graph.cpp open_svb/svb_next.
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1          // getline
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
    return -1;                 // non-ACGT -> RAW escape
}

int cmd_pack(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: svp pack out.svb < reads.fastq\n"); return 1; }
    FILE* out = std::fopen(argv[1], "wb");
    if (!out) { std::perror("fopen out"); return 1; }

    char hdr[32] = {0};
    std::memcpy(hdr, "SVB1", 4);
    uint32_t ver = 1; std::memcpy(hdr + 4, &ver, 4);
    std::fwrite(hdr, 1, 32, out);                     // placeholder counts; patched at EOF

    uint64_t nreads = 0, total = 0;
    char* line = nullptr; size_t cap = 0; ssize_t got; int ln = 0;
    std::vector<uint8_t> pb;
    while ((got = getline(&line, &cap, stdin)) != -1) {
        if (ln == 1) {                                // FASTQ sequence line
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

    std::memcpy(hdr + 8, &nreads, 8);
    std::memcpy(hdr + 16, &total, 8);
    std::fseek(out, 0, SEEK_SET);
    std::fwrite(hdr, 1, 32, out);
    std::fclose(out);
    std::fprintf(stderr, "[svb_pack] %s: %llu reads, %llu bases\n",
                 argv[1], (unsigned long long)nreads, (unsigned long long)total);
    return 0;
}
