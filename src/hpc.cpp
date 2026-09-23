#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "core.h"
#include "hpc.h"
#include "io.h"

static void hpc_work(block* b) {
    rdcore* r = &((rdcore*)b->ctx)[b->id];
    for (size_t off = 0; off < b->inlen; off += svb_size(b->in + off, b->inlen - off)) {
        collapse(r, (const uint8_t*)b->in + off);
        uint32_t nw = (r->k + 31) / 32;
        for (uint32_t j = 0; j < nw; j++) r->hs[j] = __builtin_bswap64(r->hs[j]);
        bput(b, &r->k, 4);
        bput(b, r->hs, (r->k + 3) / 4);
    }
}

int cmd_hpc(const char* in, const char* out, int t) {
    std::string path;
    if (!out) { path = std::string(in) + ".hpc"; out = path.c_str(); }
    FILE* fi = std::fopen(in, "rb");
    bin_header h;
    std::fread(&h, 1, 32, fi);
    FILE* fo = open_tmp(out);
    rdcore* r = (rdcore*)std::calloc(t, sizeof(rdcore));
    stream_run(fi, nullptr, fo, 32, t, svb_size, nullptr, hpc_work, r);
    for (int i = 0; i < t; i++) free_core(&r[i]);
    std::free(r);
    std::fclose(fi);
    bin_header o = {{'S', 'V', 'B', 'H'}, 2, h.nreads, h.total, {0}};
    std::fseek(fo, 0, SEEK_SET);
    std::fwrite(&o, 1, 32, fo);
    commit_file(fo, out);
    return 0;
}
