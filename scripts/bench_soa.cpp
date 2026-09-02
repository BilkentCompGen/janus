#include "../src/extract.cpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>
#include <utility>

static double cpu_s() { timespec t; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static double T_P1 = 0, T_DCT[8] = {0}, T_P3[8] = {0};
static uint64_t N_LV[8] = {0};
static bool CHECK = false;
static uint64_t H = 0x243F6A8885A308D3ull;
static inline void fold(uint64_t v) { H ^= v; H *= 0x9E3779B97F4A7C15ull; H ^= H >> 29; }

struct Soa { uint64_t *rep = nullptr, *st = nullptr, *en = nullptr; uint32_t *rep32 = nullptr, *sz = nullptr; int cap = 0; };
static Soa SA, SB; static uint32_t* LAB = nullptr; static uint32_t* LSZ = nullptr; static int LCAP = 0;
static void grow_soa(Soa& s, int n) {
    if (n <= s.cap) return; const int c = n < 4096 ? 4096 : n;
    s.rep = (uint64_t*)std::realloc(s.rep, (size_t)c * 8); s.st = (uint64_t*)std::realloc(s.st, (size_t)c * 8);
    s.en = (uint64_t*)std::realloc(s.en, (size_t)c * 8); s.sz = (uint32_t*)std::realloc(s.sz, (size_t)c * 4);
    s.rep32 = (uint32_t*)std::realloc(s.rep32, (size_t)c * 4); s.cap = c;
}
static void grow_lab(int n) {
    if (n <= LCAP) return; const int c = n < 4096 ? 4096 : n;
    LAB = (uint32_t*)std::realloc(LAB, (size_t)c * 4); LSZ = (uint32_t*)std::realloc(LSZ, (size_t)c * 4); LCAP = c;
}
static void fold_soa(int lv, const Soa& s, int n) {
    fold((uint64_t)lv << 32 | (uint32_t)n);
    if (lv == 1) for (int i = 0; i < n; i++) { fold(s.st[i]); fold(s.en[i]); fold(LCP_LEVEL1_TAG | (uint64_t)s.rep32[i]); fold(s.sz[i]); }
    else for (int i = 0; i < n; i++) { fold(s.st[i]); fold(s.en[i]); fold(s.rep[i]); fold(s.sz[i]); }
}
static void fold_cores(int lv, const struct core* c, int n) {
    fold((uint64_t)lv << 32 | (uint32_t)n);
    for (int i = 0; i < n; i++) { fold(c[i].start); fold(c[i].end); fold(c[i].bit_rep); fold(c[i].bit_size); }
}
static int parse_timed_struct(Scratch& s, const char* data, int n, int target) {
    struct core* cores = grow_cores(s, (int)(n / LCP_CONSTANT_FACTOR) + 8);
    double t0 = cpu_s();
    int size = parse1(data, data + n, cores, 0);
    double t1 = cpu_s(); T_P1 += t1 - t0; N_LV[1] += size;
    if (CHECK) fold_cores(1, cores, size);
    for (int lv = 1; lv < target && size; lv++) {
        struct lps p; p.level = lv; p.size = size; p.cores = cores;
        if (lcp_dct(&p) < 0) return 0;
        double t2 = cpu_s(); T_DCT[lv] += t2 - t1;
        size = parse3(cores + LCP_DCT_ITERATION_COUNT, cores + size, cores);
        t1 = cpu_s(); T_P3[lv] += t1 - t2; N_LV[lv + 1] += size;
        if (CHECK) fold_cores(lv + 1, cores, size);
    }
    return size;
}
static uint64_t FALLBACK = 0;
static int parse_timed_soa(Soa*& cur, const char* data, int n, int target) {
    const int cap = (int)(n / LCP_CONSTANT_FACTOR) + 8;
    grow_soa(SA, cap); grow_soa(SB, cap); grow_lab(cap);
    Soa* a = &SA; Soa* b = &SB;
    struct lcp_soa ia = { a->rep, a->rep32, a->st, a->en, a->sz };
    double t0 = cpu_s();
    int size = parse1_soa(data, data + n, &ia, 0);
    if (size < 0) return -1;
    double t1 = cpu_s(); T_P1 += t1 - t0; N_LV[1] += size;
    if (CHECK) fold_soa(1, *a, size);
    for (int lv = 1; lv < target && size; lv++) {
        if (size < 2) return 0;
        struct lcp_soa in = { a->rep, a->rep32, a->st, a->en, a->sz }, out = { b->rep, b->rep32, b->st, b->en, b->sz };
        lcp_dct_soa(&in, size, lv, LAB, LSZ);
        double t2 = cpu_s(); T_DCT[lv] += t2 - t1;
        size = parse3_soa(LAB, LSZ, &in, size, &out);
        t1 = cpu_s(); T_P3[lv] += t1 - t2; N_LV[lv + 1] += size;
        std::swap(a, b);
        if (CHECK) fold_soa(lv + 1, *a, size);
    }
    cur = a; return size;
}
static void emit_soa(Scratch& s, const Soa& c, int n, const std::string& h, const std::vector<uint32_t>& map, uint64_t offset) {
    grow_out(s, n);
    for (int i = 0; i < n; i++) {
        const uint64_t a = c.st[i], b = c.en[i] + 1;
        s.sig[i] = core_sig(h.data() + a, (int)(b - a));
        s.st[i] = offset + map[a];
        s.en[i] = offset + map[b < map.size() ? b : map.size() - 1];
    }
}

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: bench_soa <in.svb> <max_bases> [level=5] [check|noisy|synthetic]\n"); return 2; }
    const char* path = argv[1];
    const uint64_t max_bases = std::strtoull(argv[2], nullptr, 10);
    const int level = argc > 3 ? std::atoi(argv[3]) : 5;
    const bool synthetic = argc > 4 && std::strcmp(argv[4], "synthetic") == 0;
    CHECK = argc > 4 && (std::strcmp(argv[4], "check") == 0 || std::strcmp(argv[4], "noisy") == 0);
    const bool noisy = argc > 4 && std::strcmp(argv[4], "noisy") == 0;
    uint64_t raw_reads = 0;
    LCP_INIT();
    FILE* f = std::fopen(path, "rb");
    if (!f) { std::perror(path); return 1; }
    char hdr[32];
    if (std::fread(hdr, 1, 32, f) != 32 || std::memcmp(hdr, "SVB1", 4)) { std::fprintf(stderr, "not a .svb\n"); return 1; }
    static const char LUT[4] = { 'A', 'C', 'G', 'T' };
    static uint32_t QUAD[256];
    for (int b = 0; b < 256; b++) { const char c[4] = { LUT[(b >> 6) & 3], LUT[(b >> 4) & 3], LUT[(b >> 2) & 3], LUT[b & 3] }; std::memcpy(&QUAD[b], c, 4); }
    std::string seq, rcbuf;
    std::vector<uint8_t> pb;
    uint64_t reads = 0, bases = 0, cores_f = 0, cores_r = 0, hpc_bases = 0;
    double t_read = 0, t_hpc = 0, t_rc = 0, t_parse = 0, t_hash = 0; uint64_t rc_fallback = 0;
    while (bases < max_bases) {
        double t0 = cpu_s();
        uint32_t tag;
        if (std::fread(&tag, 4, 1, f) != 1) break;
        const uint32_t L = tag & 0x7fffffffu;
        seq.resize(L);
        if (tag >> 31) { raw_reads++; if (L && std::fread(&seq[0], 1, L, f) != L) { std::fprintf(stderr, "truncated\n"); return 1; } }
        else {
            const uint32_t nb = (L + 3) >> 2; pb.resize(nb);
            if (nb && std::fread(pb.data(), 1, nb, f) != nb) { std::fprintf(stderr, "truncated\n"); return 1; }
            const uint32_t full = L >> 2; char* out = &seq[0];
            for (uint32_t i = 0; i < full; i++) std::memcpy(out + 4 * i, &QUAD[pb[i]], 4);
            for (uint32_t i = full << 2; i < L; i++) seq[i] = LUT[(pb[i >> 2] >> (6 - 2 * (i & 3))) & 3];
        }
        if (synthetic) for (uint32_t i = 0; i < L; i++) seq[i] = LUT[i & 3];
        double t1 = cpu_s(); t_read += t1 - t0;
        reads++; bases += L;
        if (L < 3) continue;
        hpc(seq.data(), (int)L, hs, hmap);
        if (noisy) for (size_t k = 1; k < hs.size(); k++) { if (k % 397 == 0) hs[k] = hs[k - 1]; if (k % 1013 == 0) hs[k] = 'N'; }
        double t2 = cpu_s(); t_hpc += t2 - t1;
        hpc_bases += hs.size();
        int n1 = 0; Soa* c1 = nullptr;
        if ((int)hs.size() >= 3) n1 = parse_timed_soa(c1, hs.data(), (int)hs.size(), level);
        if (n1 < 0) { FALLBACK++; n1 = parse_timed_struct(sc, hs.data(), (int)hs.size(), level); }
        double t3 = cpu_s(); t_parse += t3 - t2;
        if (n1 > 0) { if (c1) emit_soa(sc, *c1, n1, hs, hmap, 0); else emit(sc, sc.cores, n1, hs, hmap, 0); cores_f += n1; if (CHECK) for (int i = 0; i < n1; i++) { fold(sc.sig[i]); fold(sc.st[i]); fold(sc.en[i]); } }
        double t4 = cpu_s(); t_hash += t4 - t3;
        if (rc_safe(hs)) { rc_in_place(hs, hmap, (int)L); }
        else { rc_fallback++; rcbuf.resize(L); for (uint32_t i = 0; i < L; i++) rcbuf[i] = (char)RC_LUT[(uint8_t)seq[L - 1 - i]]; hpc(rcbuf.data(), (int)L, hs, hmap); }
        double t5 = cpu_s(); t_rc += t5 - t4;
        int n2 = 0; Soa* c2 = nullptr;
        if ((int)hs.size() >= 3) n2 = parse_timed_soa(c2, hs.data(), (int)hs.size(), level);
        if (n2 < 0) { FALLBACK++; n2 = parse_timed_struct(sc2, hs.data(), (int)hs.size(), level); }
        double t6 = cpu_s(); t_parse += t6 - t5;
        if (n2 > 0) { if (c2) emit_soa(sc2, *c2, n2, hs, hmap, 0); else emit(sc2, sc2.cores, n2, hs, hmap, 0); cores_r += n2; if (CHECK) for (int i = 0; i < n2; i++) { fold(sc2.sig[i]); fold(sc2.st[i]); fold(sc2.en[i]); } }
        double t7 = cpu_s(); t_hash += t7 - t6;
    }
    std::fclose(f);
    const double tot = t_read + t_hpc + t_rc + t_parse + t_hash, gb = bases / 1e9;
    std::printf("%s level=%d%s  [SoA]\n", path, level, CHECK ? "  CHECK" : "");
    if (CHECK) std::printf("checksum %016llx  (raw reads %llu)\n", (unsigned long long)H, (unsigned long long)raw_reads);
    std::printf("struct fallback strands %llu\n", (unsigned long long)FALLBACK);
    std::printf("reads %llu  bases %.3f Gbp  hpc bases %.3f Gbp  rc via second hpc %llu  cores fwd %llu rc %llu (%.1f per kbp per strand)\n",
                (unsigned long long)reads, gb, hpc_bases / 1e9, (unsigned long long)rc_fallback, (unsigned long long)cores_f, (unsigned long long)cores_r, cores_f / (bases / 1e3));
    std::printf("%-24s %9s %10s %7s\n", "phase", "CPU-s", "CPU-s/Gbp", "share");
    std::printf("%-24s %9.2f %10.3f %6.1f%%\n", "svb decode", t_read, t_read / gb, 100 * t_read / tot);
    std::printf("%-24s %9.2f %10.3f %6.1f%%\n", "homopolymer compress", t_hpc, t_hpc / gb, 100 * t_hpc / tot);
    std::printf("%-24s %9.2f %10.3f %6.1f%%\n", "rc derive (in place)", t_rc, t_rc / gb, 100 * t_rc / tot);
    std::printf("%-24s %9.2f %10.3f %6.1f%%\n", "LCP parse (2 strands)", t_parse, t_parse / gb, 100 * t_parse / tot);
    std::printf("%-24s %9.2f %10.3f %6.1f%%\n", "MurmurHash64A per core", t_hash, t_hash / gb, 100 * t_hash / tot);
    std::printf("%-24s %9.2f %10.3f\n", "total", tot, tot / gb);
    std::printf("\nLCP by level (both strands)   CPU-s   cores out (M)   per kbp\n");
    std::printf("  parse1 -> level 1        %8.2f   %12.1f   %7.1f   %.2f ns per hpc char (both strands: %.1f G chars)\n", T_P1, N_LV[1] / 1e6, N_LV[1] / (2 * bases / 1e3), T_P1 / (2.0 * hpc_bases) * 1e9, 2.0 * hpc_bases / 1e9);
    for (int lv = 1; lv < level; lv++)
        std::printf("  dct+parse3 -> level %d    %8.2f + %8.2f   %12.1f   %7.1f\n", lv + 1, T_DCT[lv], T_P3[lv], N_LV[lv + 1] / 1e6, N_LV[lv + 1] / (2 * bases / 1e3));
    return 0;
}
