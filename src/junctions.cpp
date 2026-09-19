#include "cmd.h"
#include "samples.h"
#include "postings.h"
#include "log.h"
#include <sys/stat.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

static uint32_t findp(std::vector<uint32_t>& p, uint32_t x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

static ankerl::unordered_dense::set<uint64_t> load_drop(const std::string& graph) {
    ankerl::unordered_dense::set<uint64_t> d;
    const std::string p = dir_of(graph) + "/drop.sigs";
    FILE* f = std::fopen(p.c_str(), "r");
    if (!f) { say("jcontigs", "%s missing - no bubble filter applied", p.c_str()); return d; }
    unsigned long long v;
    while (std::fscanf(f, "%llu", &v) == 1) d.insert((uint64_t)v);
    std::fclose(f);
    say("jcontigs", "%s: %zu cores to skip", p.c_str(), d.size());
    return d;
}

static double chunk_vaf(const std::vector<uint32_t>& num, const std::vector<uint32_t>& den, size_t c, uint32_t L) {
    double best = 0;
    for (uint32_t l = 0; l < L; l++) {
        const uint32_t d = den[c * L + l];
        best = std::max(best, d ? (double)num[c * L + l] / d : 1.0);
    }
    return best;
}

Chunks build_chunks(const std::string& path, uint64_t& K, size_t& nsom, double min_vaf) {
    FILE* f = xopen(path.c_str(), "rb");
    uint32_t S = 0;
    if (std::fread(&S, 4, 1, f) != 1 || S < 2 || S > 64) die("jcontigs", "graph.bin sample count out of range");
    const uint32_t L = S - 1;
    if (std::fread(&K, 8, 1, f) != 1) die("jcontigs", "read K");
    std::vector<uint64_t> sig(K);
    if (std::fread(sig.data(), 8, K, f) != K) die("jcontigs", "read sigs");
    std::vector<cov_t> nc(K);
    if (std::fread(nc.data(), sizeof(cov_t), K, f) != K) die("jcontigs", "read normal plane");

    std::vector<cov_t> tc((size_t)K * L);
    if (L && std::fread(tc.data(), sizeof(cov_t), (size_t)K * L, f) != (size_t)K * L) die("jcontigs", "read tumour planes");
    uint64_t E = 0;
    if (std::fread(&E, 8, 1, f) != 1) die("jcontigs", "read E");
    std::vector<uint32_t> es(E), ed(E);
    if (std::fread(es.data(), 4, E, f) != E) die("jcontigs", "read esrc");
    if (std::fread(ed.data(), 4, E, f) != E) die("jcontigs", "read edst");
    std::fclose(f);

    const auto drop = load_drop(path);
    std::vector<uint8_t> som(K, 0);
    uint64_t ndrop = 0;
    for (uint64_t i = 0; i < K; i++) {
        som[i] = (nc[i] == 0);
        if (som[i] && drop.contains(sig[i])) { som[i] = 0; ndrop++; }
    }
    if (ndrop) say("jcontigs", "%llu somatic cores dropped by drop.sigs", (unsigned long long)ndrop);
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

    size_t noanch = 0;
    for (size_t c = 0; c < NC; c++) {
        bool any = false;
        for (uint32_t l = 0; l < L && !any; l++) any = den[c * L + l] != 0;
        noanch += !any;
    }
    say("jcontigs", "chunk VAF (chunk cov / anchor cov, best line), %zu chunks, %zu with no germline anchor",
        NC, noanch);
    for (double t : {0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50}) {
        size_t k = 0;
        for (size_t c = 0; c < NC; c++) k += (chunk_vaf(num, den, c, L) >= t);
        say("jcontigs", "  VAF >= %.2f : %9zu chunks (%5.2f%%)", t, k, 100.0 * k / NC);
    }

    std::vector<uint32_t> keep(NC, UINT32_MAX);
    uint32_t nk = 0;
    for (size_t c = 0; c < NC; c++) if (chunk_vaf(num, den, c, L) >= min_vaf) keep[c] = nk++;
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

static int32_t chunk_of(const Chunks& c, uint64_t s) {
    const auto it = c.of.find(s);
    return it == c.of.end() ? -1 : (int32_t)it->second;
}

static bool hit_cmp(const JHit& a, const JHit& b) {
    return a.read != b.read ? a.read < b.read : a.ord < b.ord;
}

int run_junctions(const std::string& gpath, const std::vector<Sample>& smp, const std::string& odir,
                  int NT, int PAD, uint32_t PER, int64_t MINF, double MINVAF) {
    const int LEVEL = level_env();

    uint64_t K = 0;
    size_t nsom = 0;
    const Chunks C = build_chunks(gpath, K, nsom, MINVAF);
    const size_t L = (size_t)C.lines;
    if (smp.size() != L + 1) die("jcontigs", "sheet and graph.bin disagree on the sample count");
    say("jcontigs", "%s: %llu nodes, %zu somatic -> %zu chunks over %zu lines; pad=%d cores, "
                   "<=%u reads/chunk/line, min flank %lld bp, min VAF %.2f, level=%d, %d threads",
        gpath.c_str(), (unsigned long long)K, nsom, C.n, L, PAD, PER, (long long)MINF, MINVAF, LEVEL, NT);
    if (!C.n) die("jcontigs", "no somatic chunks");

    mkdir(odir.c_str(), 0755);
    std::vector<std::atomic<uint32_t>> taken(C.n * L);
    for (auto& t : taken) t.store(0);
    std::vector<std::atomic<uint8_t>> touched(C.n * L);
    for (auto& t : touched) t.store(0);

    Index ix;
    ix.load(dir_of(gpath) + "/index.bin");
    if (ix.K != K) die("jcontigs", "index.bin and graph.bin disagree on the node count");
    if ((int)ix.pad != PAD) die("jcontigs", "index.bin was built with a different pad_cores");
    say("jcontigs", "index.bin: %llu records, N=%u, pad=%u",
        (unsigned long long)ix.off[ix.K], ix.N, ix.pad);

    std::vector<std::vector<JHit>> byline(L);
    {
        std::vector<uint64_t> sig(K);
        FILE* f = xopen(gpath.c_str(), "rb");
        std::fseek(f, 12, SEEK_SET);
        if (K && std::fread(sig.data(), 8, K, f) != K) die("jcontigs", "read sigs");
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
        SvbFile rd;
        rd.open(sm.path);
        FILE* out = xopen((opath + ".tmp").c_str(), "w");
        uint64_t emitted = 0, emit_bp = 0, short_flank = 0, nreads = 0;

        std::vector<JHit>& h = byline[l];
        std::sort(h.begin(), h.end(), hit_cmp);

        std::string seq, buf;
        for (size_t a0 = 0; a0 < h.size();) {
            size_t e0 = a0;
            while (e0 < h.size() && h[e0].read == h[a0].read) e0++;
            nreads++;
            int w = 0;
            bool have_seq = false;
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
                wr(out, buf.data(), buf.size());
                emitted++;
                emit_bp += (uint64_t)(we - ws);
            }
            a0 = e0;
        }
        atomic_close(out, (opath + ".tmp").c_str(), opath.c_str());
        rd.close();
        say("jcontigs", "%-8s %llu reads taken -> %llu contigs (%.2f Gbp) -> %s",
            sm.name.c_str(), (unsigned long long)nreads,
            (unsigned long long)emitted, emit_bp / 1e9, opath.c_str());
        sf_total += short_flank;
    }

    size_t got = 0, seen_ch = 0;
    for (size_t c = 0; c < C.n; c++) {
        bool g = false, t = false;
        for (size_t l = 0; l < L; l++) { g |= taken[c * L + l].load() > 0; t |= touched[c * L + l].load() != 0; }
        got += g;
        seen_ch += t;
    }
    say("jcontigs", "chunks: %zu total, %zu carried by some read, %zu with a contig "
                   "(%zu carried but NO read had %lld bp of anchor both sides -- these are the ones a "
                   "cross-read merge would be for), %llu windows rejected for short flanks",
        C.n, seen_ch, got, seen_ch - got, (long long)MINF, (unsigned long long)sf_total);
    return 0;
}

int cmd_junctions(int argc, const char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: snape junctions <graph.bin> <sheet.tsv> <outdir> "
                                 "[threads=16] [pad_cores=50] [per_chunk=5] [min_flank=300] [min_vaf=0]\n"); return 2; }
    return run_junctions(argv[1], samples_read(argv[2]), argv[3],
                          argi(argc, argv, 4, 16), argi(argc, argv, 5, 50),
                          (uint32_t)argi(argc, argv, 6, 5), (int64_t)argll(argc, argv, 7, 300),
                          argd(argc, argv, 8, 0.0));
}
