#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include "hpc.h"
#include "io.h"

static size_t hpc(const uint8_t* p, uint32_t tag, int len, char* out) {
    size_t k = 0;
    int prev = -1;
    if (tag & 0x80000000) {
        for (int i = 0; i < len; i++) {
            int v = CODE.t[(uint8_t)(p[i] & 0xdf)];
            if (v == prev) continue;
            prev = v;
            if ((k & 3) == 0) out[k >> 2] = 0;
            out[k >> 2] |= v << (6 - 2 * (k & 3));
            k++;
        }
    } else {
        for (int i = 0; i < len; i++) {
            int v = (p[i >> 2] >> (6 - 2 * (i & 3))) & 3;
            if (v == prev) continue;
            prev = v;
            if ((k & 3) == 0) out[k >> 2] = 0;
            out[k >> 2] |= v << (6 - 2 * (k & 3));
            k++;
        }
    }
    return k;
}

static void hpc_work(block* b) {
    for (size_t off = 0; off < b->inlen;) {
        const uint8_t* rec = (const uint8_t*)b->in + off;
        uint32_t tag;
        std::memcpy(&tag, rec, 4);
        uint32_t len = tag & 0x7fffffff;
        bgrow(b, (len + 3) / 4);
        size_t k = hpc(rec + 4, tag, len, b->buf);
        uint32_t klen = (uint32_t)k;
        bput(b, &klen, 4);
        bput(b, b->buf, (k + 3) / 4);
        off += svb_size(b->in + off, b->inlen - off);
    }
}

int cmd_hpc(const char* in, const char* out, int t) {
    std::string path;
    if (!out) { path = std::string(in) + ".hpc"; out = path.c_str(); }
    FILE* fi = std::fopen(in, "rb");
    bin_header h;
    std::fread(&h, 1, 32, fi);
    FILE* fo = open_tmp(out);
    stream_run(fi, nullptr, fo, 32, t, svb_size, nullptr, hpc_work, nullptr);
    std::fclose(fi);
    bin_header o = {{'S', 'V', 'B', 'H'}, 2, h.nreads, h.total, {0}};
    std::fseek(fo, 0, SEEK_SET);
    std::fwrite(&o, 1, 32, fo);
    commit_file(fo, out);
    return 0;
}
