#include <algorithm>
#include <cstdlib>
#include "cmd.h"
#include "core.h"
#include "io.h"

Parser::~Parser() {
    std::free(mem);
}

void Parser::grow(Read& r) {
    uint32_t k = r.k;
    if (r.a.size() < 2 * k + 2) { r.a.resize(4 * k + 4); r.b.resize(4 * k + 4); }
    if (map.size() > k) return;
    k *= 2;
    uint64_t c = (uint64_t)(k / LCP_CONSTANT_FACTOR) + 8, *p = mem = (uint64_t*)std::realloc(mem, c * 72);
    for (int i = 0; i < 2; i++, p += 4 * c) o[i] = {p, (uint32_t*)(p + 3 * c), p + c, p + 2 * c, (uint32_t*)(p + 3 * c) + c};
    lab = (uint32_t*)p;
    lsz = lab + c;
    map.resize(k + 1);
    hs.resize(k);
}

uint32_t Parser::parse(Read& r, uint32_t at) {
    int n = parse1_soa(hs.data(), hs.data() + r.k, &o[0], 0), c = 0;
    for (int v = 1; v < LEV && n > 0; v++, c ^= 1) {
        lcp_dct_soa(&o[c], n, v, lab, lsz);
        n = parse3_soa(lab, lsz, &o[c], n, &o[c ^ 1]);
    }
    for (int i = 0; i < n; i++) { r.a[at + i] = map[o[c].st[i]]; r.b[at + i] = map[o[c].en[i]]; }
    return n > 0 ? n : 0;
}

void Parser::cores(Read& r, const uint8_t* rec) {
    r.squeeze(rec);
    grow(r);
    r.starts(map.data());
    r.n = r.nf = 0;
    if (r.k < 3) return;
    r.letters(false, hs.data());
    r.nf = parse(r, 0);
    r.letters(true, hs.data());
    std::reverse(map.begin(), map.begin() + r.k + 1);
    for (uint32_t i = 0; i <= r.k; i++) map[i] = r.len - map[i];
    r.n = r.nf + parse(r, r.nf);
}

class Extract : public Stream {
public:
    explicit Extract(const Args& a);
    int run();

private:
    const Args& a;
    std::vector<Parser> parser;
    std::vector<Read> read;
    std::vector<uint64_t> cores;
    void block(Block& b) override;
};

static const int W = LEV + 4;

Extract::Extract(const Args& a) : a(a), parser(a.threads), read(a.threads), cores(a.threads, 0) {}

void Extract::block(Block& b) {
    Read& r = read[b.id];
    for (size_t off = 0; off < b.inlen[0]; off += svb_size(b.in[0] + off, b.inlen[0] - off)) {
        parser[b.id].cores(r, (const uint8_t*)b.in[0] + off);
        r.store(b.out, W);
        cores[b.id] += r.n;
    }
}

int Extract::run() {
    Output out(a.out);
    Input in(a.in);
    Header h('X', 3, pass(in.f, nullptr, out.f, a.threads), 0);
    h.resv[0] = LEV;
    h.resv[1] = W;
    for (uint64_t c : cores) h.total += c;
    out.finish(h);
    return 0;
}

int cmd_extract(const Args& a) {
    return Extract(a).run();
}
