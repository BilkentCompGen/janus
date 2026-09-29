#include <zlib.h>
#include <string>
#include <vector>
#include "io.h"
#include "cmd.h"

static bool line(gzFile g, std::string& s, char* b) {
    s.clear();
    while (gzgets(g, b, 1 << 16)) {
        s += b;
        if (s.back() == '\n') break;
    }
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return !s.empty() || !gzeof(g);
}

int cmd_pack(const char* in, const char* out) {
    gzFile g = gzopen(in, "rb");
    gzbuffer(g, 1 << 22);
    std::string x = std::string(out) + ".idx", s, q;
    FILE *fo = open_tmp(out), *fx = open_tmp(x.c_str());
    std::fseek(fo, sizeof(bin_header), SEEK_SET);
    std::fseek(fx, sizeof(bin_header), SEEK_SET);
    std::vector<char> b(1 << 16);
    std::vector<uint8_t> p;
    uint64_t n = 0, tot = 0, pos = sizeof(bin_header), mx = 0;
    while (line(g, q, b.data()) && line(g, s, b.data()) && line(g, q, b.data()) && line(g, q, b.data())) {
        uint32_t L = s.size(), raw = 0;
        for (char c : s) raw |= (c & 0xdf) != 'A' && (c & 0xdf) != 'C' && (c & 0xdf) != 'G' && (c & 0xdf) != 'T';
        uint32_t tag = L | raw << 31;
        p.assign(raw ? L : (L + 3) / 4, 0);
        for (uint32_t i = 0; i < L; i++) if (raw) p[i] = s[i]; else p[i >> 2] |= CODE.t[s[i] & 0xdf] << (6 - 2 * (i & 3));
        std::fwrite(&pos, 8, 1, fx);
        std::fwrite(&tag, 4, 1, fo);
        std::fwrite(p.data(), 1, p.size(), fo);
        pos += 4 + p.size();
        tot += L;
        mx = std::max<uint64_t>(mx, L);
        n++;
    }
    std::fwrite(&pos, 8, 1, fx);
    gzclose(g);
    bin_header ho = {{'S', 'V', 'B', '1'}, 1, n, tot, {}}, hx = {{'S', 'V', 'B', 'I'}, 1, n, mx, {}};
    finish(fo, out, &ho);
    finish(fx, x.c_str(), &hx);
    std::fprintf(stderr, "[pack] %lu reads, %lu bases\n", n, tot);
    return 0;
}
