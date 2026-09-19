#include "cmd.h"
#include "postings.h"
#include "log.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
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

int cmd_pack(int argc, const char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: snape pack out.svb < reads.fastq\n"); return 1; }
    FILE* out = xopen(argv[1], "wb");

    char hdr[16] = {0};
    wr(out, hdr, 16);

    uint64_t nreads = 0, total = 0;
    char* line = nullptr;
    size_t cap = 0;
    ssize_t got;
    int ln = 0;
    std::vector<uint8_t> pb;
    while ((got = getline(&line, &cap, stdin)) != -1) {
        if (ln == 1) {
            size_t L = (size_t)got;
            while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) L--;
            bool raw = false;
            for (size_t i = 0; i < L; i++) if (code(line[i]) < 0) { raw = true; break; }
            uint32_t tag = (uint32_t)L | (raw ? 0x80000000u : 0u);
            wr(out, &tag, 4);
            if (raw) {
                wr(out, line, L);
            } else {
                uint32_t nb = (uint32_t)((L + 3) >> 2);
                pb.assign(nb, 0);
                for (size_t i = 0; i < L; i++) pb[i >> 2] |= (uint8_t)(code(line[i]) << (6 - 2 * (i & 3)));
                wr(out, pb.data(), nb);
            }
            nreads++;
            total += L;
        }
        ln = (ln + 1) & 3;
    }
    std::free(line);

    std::memcpy(hdr, &nreads, 8);
    std::memcpy(hdr + 8, &total, 8);
    std::fseek(out, 0, SEEK_SET);
    wr(out, hdr, 16);
    std::fclose(out);
    say("svb_pack", "%s: %llu reads, %llu bases", argv[1], (unsigned long long)nreads, (unsigned long long)total);
    return 0;
}

int cmd_dump(int argc, const char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: snape dump in.svb > reads.fastq\n"); return 1; }
    FILE* f = xopen(argv[1], "rb");
    if (std::fseek(f, 16, SEEK_SET) != 0) die("svb", "bad .svb %s", argv[1]);

    std::string seq, qual;
    uint64_t idx = 0;
    while (svb_read_seq(f, seq)) {
        qual.assign(seq.size(), 'I');
        std::printf("@r%llu\n%s\n+\n%s\n", (unsigned long long)idx++, seq.c_str(), qual.c_str());
    }
    std::fclose(f);
    return 0;
}

int cmd_svbidx(int argc, const char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: snape svbidx <in.svb>\n"); return 1; }
    const std::string in = argv[1];
    FILE* f = xopen(in.c_str(), "rb");

    char hdr[16];
    if (std::fread(hdr, 1, 16, f) != 16) die("svbidx", "%s is not a .svb", in.c_str());
    uint64_t nreads, total;
    std::memcpy(&nreads, hdr, 8);
    std::memcpy(&total, hdr + 8, 8);

    std::vector<uint64_t> off;
    off.reserve(nreads + 1);

    constexpr size_t BLOCK = 1u << 22;
    std::vector<uint8_t> buf(BLOCK);
    uint64_t pos = 16, seen = 0, bases = 0;
    size_t have = 0, at = 0;
    uint64_t maxlen = 0;

    for (;;) {
        if (have - at < 4) {
            std::memmove(buf.data(), buf.data() + at, have - at);
            have -= at;
            at = 0;
            const size_t got = std::fread(buf.data() + have, 1, BLOCK - have, f);
            have += got;
            if (have - at < 4) break;
        }
        uint32_t tag;
        std::memcpy(&tag, buf.data() + at, 4);
        const uint64_t L = tag & 0x7fffffffu;
        const bool raw = (tag & 0x80000000u) != 0;
        const uint64_t payload = raw ? L : ((L + 3) >> 2);
        off.push_back(pos);
        seen++;
        bases += L;
        if (L > maxlen) maxlen = L;
        pos += 4 + payload;
        const uint64_t skip = 4 + payload;
        if (have - at >= skip) {
            at += (size_t)skip;
        } else {
            const uint64_t rest = skip - (have - at);
            at = have = 0;
            if (std::fseek(f, (long)rest, SEEK_CUR) != 0) die("svbidx", "seek past record %llu", (unsigned long long)seen);
        }
    }
    std::fclose(f);
    off.push_back(pos);

    if (seen != nreads) die("svbidx", "header says %llu reads, scan found %llu", (unsigned long long)nreads, (unsigned long long)seen);
    if (bases != total) die("svbidx", "header says %llu bases, scan found %llu", (unsigned long long)total, (unsigned long long)bases);

    const std::string out = in + ".idx", tmp = out + ".tmp";
    FILE* o = xopen(tmp.c_str(), "wb");
    char oh[16] = {0};
    std::memcpy(oh, &nreads, 8);
    std::memcpy(oh + 8, &maxlen, 8);
    wr(o, oh, 16);
    wr(o, off.data(), off.size() * 8);
    atomic_close(o, tmp.c_str(), out.c_str());
    say("svbidx", "%s: %llu reads, longest %llu bp -> %s (%.1f MB)",
        in.c_str(), (unsigned long long)nreads, (unsigned long long)maxlen,
        out.c_str(), off.size() * 8 / 1e6);
    return 0;
}
