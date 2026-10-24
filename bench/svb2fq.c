// Packed reads (SVB1) back to FASTQ with the PacBio adapters that Badread adds trimmed off; used for the simulated normal.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    FILE* f = fopen(argv[1], "rb");
    const char* name = argc > 2 ? argv[2] : "r";
    char hdr[32];
    if (!f || fread(hdr, 1, 32, f) != 32 || memcmp(hdr, "SVB1", 4)) { fprintf(stderr, "not an SVB1 file\n"); return 1; }
    static char out[1 << 20];
    setvbuf(stdout, out, _IOFBF, sizeof out);
    const char B[4] = {'A', 'C', 'G', 'T'};
    size_t cap = 1 << 20;
    uint8_t* p = malloc(cap);
    char* s = malloc(cap * 4 + 1);
    char* q = malloc(cap * 4 + 1);
    uint64_t n = 0;
    uint32_t tag;
    while (fread(&tag, 4, 1, f) == 1) {
        uint32_t L = tag & 0x7fffffffu, raw = tag >> 31, nb = raw ? L : (L + 3) >> 2;
        if (nb > cap || L > cap * 4) {
            cap = (nb > L ? nb : L) * 2;
            p = realloc(p, cap); s = realloc(s, cap * 4 + 1); q = realloc(q, cap * 4 + 1);
        }
        if (fread(p, 1, nb, f) != nb) { fprintf(stderr, "truncated at read %lu\n", (unsigned long)n); return 1; }
        if (raw) memcpy(s, p, L);
        else for (uint32_t i = 0; i < L; i++) s[i] = B[(p[i >> 2] >> (6 - 2 * (i & 3))) & 3];
        char* r = s; uint32_t n = L;
        const char* A = "AATGTACTTCGTTCAGTTACGTATTGCT";
        const char* E = "GCAATACGTAACTGAACGAAGTACATT";
        uint32_t cut = 0;
        for (uint32_t j = 0; j + 21 <= n && j < 100; j++) if (!memcmp(r + j, A + 7, 21)) { cut = j + 21; break; }
        if (!cut) for (uint32_t k = 20; k >= 8; k--) if (k <= n && !memcmp(r, A + 28 - k, k)) { cut = k; break; }
        r += cut; n -= cut;
        uint32_t end = n;
        for (uint32_t j = n > 100 ? n - 100 : 0; j + 21 <= n; j++) if (!memcmp(r + j, E, 21)) { end = j; break; }
        if (end == n) for (uint32_t k = 20; k >= 8; k--) if (k <= n && !memcmp(r + n - k, E, k)) { end = n - k; break; }
        n = end;
        if (n < 100) continue;
        memmove(s, r, n); L = n;
        memset(q, 'I', L);
        printf("@%s_%lu\n", name, (unsigned long)n++);
        fwrite(s, 1, L, stdout); fputs("\n+\n", stdout);
        fwrite(q, 1, L, stdout); fputc('\n', stdout);
    }
    fprintf(stderr, "%lu reads\n", (unsigned long)n);
    return 0;
}
