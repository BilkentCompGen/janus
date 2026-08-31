#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "svb.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

int cmd_svbidx(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: snape svbidx <in.svb>\n"); return 1; }
    const std::string in = argv[1];
    FILE* f = std::fopen(in.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "[svbidx] FATAL: open %s\n", in.c_str()); return 1; }

    char hdr[32];
    if (std::fread(hdr, 1, 32, f) != 32 || std::memcmp(hdr, "SVB1", 4)) {
        std::fprintf(stderr, "[svbidx] FATAL: %s is not a .svb\n", in.c_str()); return 1;
    }
    uint32_t ver; std::memcpy(&ver, hdr + 4, 4);
    if (ver != 1) { std::fprintf(stderr, "[svbidx] FATAL: .svb version %u != 1\n", ver); return 1; }
    uint64_t nreads, total;
    std::memcpy(&nreads, hdr + 8, 8);
    std::memcpy(&total, hdr + 16, 8);

    std::vector<uint64_t> off;
    off.reserve(nreads + 1);

    constexpr size_t BLOCK = 1u << 22;
    std::vector<uint8_t> buf(BLOCK);
    uint64_t pos = 32, seen = 0, bases = 0;
    size_t have = 0, at = 0;
    uint64_t maxlen = 0;

    for (;;) {
        if (have - at < 4) {
            std::memmove(buf.data(), buf.data() + at, have - at);
            have -= at; at = 0;
            const size_t got = std::fread(buf.data() + have, 1, BLOCK - have, f);
            have += got;
            if (have - at < 4) break;
        }
        uint32_t tag; std::memcpy(&tag, buf.data() + at, 4);
        const uint64_t L = tag & 0x7fffffffu;
        const bool raw = (tag & 0x80000000u) != 0;
        const uint64_t payload = raw ? L : ((L + 3) >> 2);
        off.push_back(pos);
        seen++; bases += L;
        if (L > maxlen) maxlen = L;
        pos += 4 + payload;
        const uint64_t skip = 4 + payload;
        if (have - at >= skip) {
            at += (size_t)skip;
        } else {
            const uint64_t rest = skip - (have - at);
            at = have = 0;
            if (std::fseek(f, (long)rest, SEEK_CUR) != 0) {
                std::fprintf(stderr, "[svbidx] FATAL: seek past record %llu\n", (unsigned long long)seen);
                return 1;
            }
        }
    }
    std::fclose(f);
    off.push_back(pos);

    if (seen != nreads) {
        std::fprintf(stderr, "[svbidx] FATAL: header says %llu reads, scan found %llu\n",
                     (unsigned long long)nreads, (unsigned long long)seen);
        return 1;
    }
    if (bases != total) {
        std::fprintf(stderr, "[svbidx] FATAL: header says %llu bases, scan found %llu\n",
                     (unsigned long long)total, (unsigned long long)bases);
        return 1;
    }

    const std::string out = in + ".idx", tmp = out + ".tmp";
    FILE* o = std::fopen(tmp.c_str(), "wb");
    if (!o) { std::fprintf(stderr, "[svbidx] FATAL: open %s\n", tmp.c_str()); return 1; }
    char oh[32] = {0};
    std::memcpy(oh, "SVBI", 4);
    const uint32_t ov = 1; std::memcpy(oh + 4, &ov, 4);
    std::memcpy(oh + 8, &nreads, 8);
    std::memcpy(oh + 16, &maxlen, 8);
    if (std::fwrite(oh, 1, 32, o) != 32 ||
        std::fwrite(off.data(), 8, off.size(), o) != off.size()) {
        std::fprintf(stderr, "[svbidx] FATAL: write %s\n", tmp.c_str()); return 1;
    }
    std::fclose(o);
    if (std::rename(tmp.c_str(), out.c_str())) {
        std::fprintf(stderr, "[svbidx] FATAL: rename %s\n", out.c_str()); return 1;
    }
    std::fprintf(stderr, "[svbidx] %s: %llu reads, longest %llu bp -> %s (%.1f MB)\n",
                 in.c_str(), (unsigned long long)nreads, (unsigned long long)maxlen,
                 out.c_str(), off.size() * 8 / 1e6);
    return 0;
}
