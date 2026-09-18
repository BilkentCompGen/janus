#include "junctions.h"
#include "samples.h"
#include "postings.h"
#include <sys/stat.h>

#include "snape.h"
#include "unordered_dense.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {
int LEVEL = 5;
constexpr int MARGIN = 10;
constexpr size_t BATCH = 2000;

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[jcontigs] FATAL: %s\n", m); std::exit(1); }

struct Chunks {
    std::vector<uint64_t> carry;
    int lines = 0;

    ankerl::unordered_dense::map<uint64_t, uint32_t> of;
    size_t n = 0;
};

uint32_t findp(std::vector<uint32_t>& p, uint32_t x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

ankerl::unordered_dense::set<uint64_t> load_drop(const std::string& graph) {
    ankerl::unordered_dense::set<uint64_t> d;
    const size_t cut = graph.find_last_of('/');
    const std::string p = (cut == std::string::npos ? std::string(".") : graph.substr(0, cut)) + "/drop.sigs";
    FILE* f = std::fopen(p.c_str(), "r");
    if (!f) { std::fprintf(stderr, "[jcontigs] %s missing - no bubble filter applied\n", p.c_str()); return d; }
    unsigned long long v;
    while (std::fscanf(f, "%llu", &v) == 1) d.insert((uint64_t)v);
    std::fclose(f);
    std::fprintf(stderr, "[jcontigs] %s: %zu cores to skip\n", p.c_str(), d.size());
    return d;
}

Chunks build_chunks(const std::string& path, uint64_t& K, size_t& nsom, double min_vaf) {
    FILE* f = std::fopen(path.c_str(), "rb"); if (!f) die("open graph.bin");
    uint32_t S = 0;
    if (std::fread(&S, 4, 1, f) != 1 || S < 2 || S > 64) die("graph.bin sample count out of range");
    const uint32_t L = S - 1;
    if (std::fread(&K, 8, 1, f) != 1) die("read K");
    std::vector<uint64_t> sig(K);
    if (std::fread(sig.data(), 8, K, f) != K) die("read sigs");
    std::vector<cov_t> nc(K);
    if (std::fread(nc.data(), sizeof(cov_t), K, f) != K) die("read normal plane");

    std::vector<cov_t> tc((size_t)K * L);
    if (L && std::fread(tc.data(), sizeof(cov_t), (size_t)K * L, f) != (size_t)K * L) die("read tumour planes");
    uint64_t E = 0;
    if (std::fread(&E, 8, 1, f) != 1) die("read E");
    std::vector<uint32_t> es(E), ed(E);
    if (std::fread(es.data(), 4, E, f) != E) die("read esrc");
    if (std::fread(ed.data(), 4, E, f) != E) die("read edst");
    std::fclose(f);

    const auto drop = load_drop(path);
    std::vector<uint8_t> som(K, 0);
    uint64_t ndrop = 0;
    for (uint64_t i = 0; i < K; i++) {

        som[i] = (nc[i] == 0);
        if (som[i] && drop.contains(sig[i])) { som[i] = 0; ndrop++; }
    }
    if (ndrop) std::fprintf(stderr, "[jcontigs] %llu somatic cores dropped by drop.sigs\n", (unsigned long long)ndrop);
    std::vector<uint32_t> par(K);
    for (uint64_t i = 0; i < K; i++) par[i] = (uint32_t)i;
    for (uint64_t e = 0; e < E; e++) {
        const uint32_t a = es[e], b = ed[e];
        if (!som[a] || !som[b]) continue;
        const uint32_t ra = findp(par, a), rb = findp(par, b);
        if (ra != rb) par[ra] = rb;
    }
    std::unordered_map<uint32_t, uint32_t> id;
    std::vector<uint32_t> cid(K, UINT32_MAX);
    for (uint64_t i = 0; i < K; i++) {
        if (!som[i]) continue;
        const uint32_t r = findp(par, (uint32_t)i);
        auto it = id.find(r);
        if (it == id.end()) it = id.emplace(r, (uint32_t)id.size()).first;
        cid[i] = it->second;
    }
    const size_t NC = id.size();

    std::vector<uint32_t> num((size_t)NC * L, 0), den((size_t)NC * L, 0);
    std::vector<uint64_t> carry(NC, 0);
    for (uint64_t i = 0; i < K; i++) {
        if (!som[i]) continue;
        const size_t c = cid[i];
        for (uint32_t l = 0; l < L; l++) {
            const cov_t v = tc[(size_t)l * K + i];
            if (v) carry[c] |= (1ull << l);
            if (v > num[c * L + l]) num[c * L + l] = v;
        }
    }
    for (uint64_t e = 0; e < E; e++) {
        const uint32_t a = es[e], b = ed[e];
        if (som[a] == som[b]) continue;
        const uint32_t sm = som[a] ? a : b, o = som[a] ? b : a;
        const size_t c = cid[sm];
        for (uint32_t l = 0; l < L; l++) {
            const cov_t v = tc[(size_t)l * K + o];
            if (v > den[c * L + l]) den[c * L + l] = v;
        }
    }

    auto vaf = [&](size_t c) {
        double best = 0;
        for (uint32_t l = 0; l < L; l++) {
            const uint32_t d = den[c * L + l];
            best = std::max(best, d ? (double)num[c * L + l] / d : 1.0);
        }
        return best;
    };

    size_t noanch = 0;
    for (size_t c = 0; c < NC; c++) {
        bool any = false;
        for (uint32_t l = 0; l < L && !any; l++) any = den[c * L + l] != 0;
        noanch += !any;
    }
    std::fprintf(stderr, "[jcontigs] chunk VAF (chunk cov / anchor cov, best line), %zu chunks, %zu with no germline anchor\n",
                 NC, noanch);
    for (double t : {0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50}) {
        size_t k = 0;
        for (size_t c = 0; c < NC; c++) k += (vaf(c) >= t);
        std::fprintf(stderr, "[jcontigs]   VAF >= %.2f : %9zu chunks (%5.2f%%)\n", t, k, 100.0 * k / NC);
    }

    std::vector<uint32_t> keep(NC, UINT32_MAX);
    uint32_t nk = 0;
    for (size_t c = 0; c < NC; c++) if (vaf(c) >= min_vaf) keep[c] = nk++;
    Chunks c;
    c.carry.assign(nk, 0);
    for (uint64_t i = 0; i < K; i++) {
        if (!som[i] || keep[cid[i]] == UINT32_MAX) continue;
        c.of.emplace(sig[i], keep[cid[i]]);
        c.carry[keep[cid[i]]] = carry[cid[i]];
    }
    nsom = c.of.size();
    c.n = nk;
    c.lines = (int)L;
    return c;
}

int32_t chunk_of(const Chunks& c, uint64_t s) {
    const auto it = c.of.find(s);
    return it == c.of.end() ? -1 : (int32_t)it->second;
}

struct Read { std::string name, seq; };
}

int cmd_junctions(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: snape junctions <graph.bin> <sheet.tsv> <outdir> "
                                 "[threads=16] [pad_cores=50] [per_chunk=5] [min_flank=300] [min_vaf=0]\n"); return 2; }
    const std::string gpath = argv[1], spath = argv[2], odir = argv[3];
    const int NT = argc > 4 ? std::atoi(argv[4]) : 16;

    const int PAD = argc > 5 ? std::atoi(argv[5]) : 50;
    const uint32_t PER = argc > 6 ? (uint32_t)std::atoi(argv[6]) : 5;
    const int64_t MINF = argc > 7 ? std::atoll(argv[7]) : 300;
    const double MINVAF = argc > 8 ? std::atof(argv[8]) : 0.0;
    if (const char* lv = std::getenv("SNAPE_LEVEL")) LEVEL = std::atoi(lv);

    const std::vector<Sample> smp = samples_read(spath);
    uint64_t K = 0; size_t nsom = 0;

    const Chunks C = build_chunks(gpath, K, nsom, MINVAF);
    const size_t L = (size_t)C.lines;
    if ((size_t)smp.size() != L + 1) die("sheet and graph.bin disagree on the sample count");
    std::fprintf(stderr, "[jcontigs] %s: %llu nodes, %zu somatic -> %zu chunks over %zu lines; pad=%d cores, "
                 "<=%u reads/chunk/line, min flank %lld bp, min VAF %.2f, level=%d, %d threads\n",
                 gpath.c_str(), (unsigned long long)K, nsom, C.n, L, PAD, PER, (long long)MINF, MINVAF, LEVEL, NT);
    if (!C.n) die("no somatic chunks");

    mkdir(odir.c_str(), 0755);
    std::vector<std::atomic<uint32_t>> taken(C.n * L);
    for (auto& t : taken) t.store(0);
    std::vector<std::atomic<uint8_t>> touched(C.n * L);
    for (auto& t : touched) t.store(0);

    Index ix;
    {
        const size_t cut = gpath.find_last_of('/');
        ix.load((cut == std::string::npos ? std::string(".") : gpath.substr(0, cut)) + "/index.bin");
    }
    if (ix.K != K) die("index.bin and graph.bin disagree on the node count");
    if ((int)ix.pad != PAD) die("index.bin was built with a different pad_cores");
    std::fprintf(stderr, "[jcontigs] index.bin: %llu records, N=%u, pad=%u\n",
                 (unsigned long long)ix.off[ix.K], ix.N, ix.pad);

    struct Hit { uint32_t read, ord, chunk, st, clen, dl, dr; };
    std::vector<std::vector<Hit>> byline(L);
    {
        std::vector<uint64_t> sig(K);
        FILE* f = std::fopen(gpath.c_str(), "rb");
        if (!f) die("open graph.bin");
        std::fseek(f, 12, SEEK_SET);
        if (K && std::fread(sig.data(), 8, K, f) != K) die("read sigs");
        std::fclose(f);
        for (uint64_t i = 0; i < K; i++) {
            const int32_t ch = chunk_of(C, sig[i]);
            if (ch < 0) continue;
            const uint32_t n = ix.count((uint32_t)i);
            for (uint32_t j = 0; j < n; j++) {
                const Post p = ix.get((uint32_t)i, j);
                if (p.samp == 0 || (size_t)p.samp > L) continue;
                byline[p.samp - 1].push_back({p.read, p.ord, (uint32_t)ch, p.st, p.clen, p.dl, p.dr});
            }
        }
    }

    uint64_t sf_total = 0;
    for (size_t l = 0; l < L; l++) {
    const Sample& sm = smp[l + 1];
    const std::string opath = odir + "/ctg." + sm.name + ".fa";
    SvbFile rd; rd.open(sm.path);
    FILE* out = std::fopen((opath + ".tmp").c_str(), "w"); if (!out) die("open out");
    uint64_t emitted = 0, emit_bp = 0, short_flank = 0, nreads = 0;

    std::vector<Hit>& h = byline[l];
    std::sort(h.begin(), h.end(), [](const Hit& a, const Hit& b) {
        return a.read != b.read ? a.read < b.read : a.ord < b.ord;
    });

    std::string seq, buf;
    for (size_t a0 = 0; a0 < h.size();) {
        size_t e0 = a0;
        while (e0 < h.size() && h[e0].read == h[a0].read) e0++;
        nreads++;
        int w = 0; bool have_seq = false;
        for (size_t a = a0; a < e0;) {
            size_t b = a;
            for (;;) {
                while (b + 1 < e0 && h[b + 1].ord == h[b].ord + 1) b++;
                size_t nxt = b + 1;
                if (nxt >= e0 || h[nxt].ord - h[b].ord > (uint32_t)(2 * PAD)) break;
                b = nxt;
            }
            const uint32_t ch = h[a].chunk;
            touched[(size_t)ch * L + l].store(1, std::memory_order_relaxed);
            const int64_t ws = (int64_t)h[a].st - (int64_t)h[a].dl;
            const int64_t we = (int64_t)h[b].st + (int64_t)h[b].clen + (int64_t)h[b].dr;
            const int64_t lf = (int64_t)h[a].dl, rf = (int64_t)h[b].dr;
            a = b + 1;
            if (lf < MINF || rf < MINF) { short_flank++; continue; }
            if (taken[(size_t)ch * L + l].fetch_add(1) >= PER) continue;
            if (!have_seq) { if (!rd.read(h[a0].read, seq)) break; have_seq = true; }
            if (ws < 0 || (uint64_t)we > seq.size()) continue;
            char hd[256];
            std::snprintf(hd, sizeof hd, ">%s|r%u_%d chunk=%u len=%lld flank=%lld,%lld\n",
                          sm.name.c_str(), h[a0].read, w++, ch,
                          (long long)(we - ws), (long long)lf, (long long)rf);
            buf.assign(hd);
            buf.append(seq, (size_t)ws, (size_t)(we - ws));
            buf += "\n";
            std::fwrite(buf.data(), 1, buf.size(), out);
            emitted++; emit_bp += (uint64_t)(we - ws);
        }
        a0 = e0;
    }
    rd.close(); std::fclose(out);
    if (std::rename((opath + ".tmp").c_str(), opath.c_str())) die("rename");
    std::fprintf(stderr, "[jcontigs] %-8s %llu reads taken -> %llu contigs (%.2f Gbp) -> %s\n",
                 sm.name.c_str(), (unsigned long long)nreads,
                 (unsigned long long)emitted, emit_bp / 1e9, opath.c_str());
    sf_total += short_flank;
    }

    size_t got = 0, seen_ch = 0;
    for (size_t c = 0; c < C.n; c++) {
        bool g = false, t = false;
        for (size_t l = 0; l < L; l++) { g |= taken[c * L + l].load() > 0; t |= touched[c * L + l].load() != 0; }
        got += g; seen_ch += t;
    }
    std::fprintf(stderr, "[jcontigs] chunks: %zu total, %zu carried by some read, %zu with a contig "
                 "(%zu carried but NO read had %lld bp of anchor both sides -- these are the ones a "
                 "cross-read merge would be for), %llu windows rejected for short flanks\n",
                 C.n, seen_ch, got, seen_ch - got, (long long)MINF, (unsigned long long)sf_total);
    return 0;
}
