#include "cmd.h"
#include "samples.h"
#include "postings.h"
#include "log.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <thread>
#include <vector>

static void edge_chunk(FILE* f, long spos, long dpos, uint64_t done, size_t m, uint32_t* bs, uint32_t* bd) {
    if (std::fseek(f, spos + (long)(done * 4), SEEK_SET) != 0) die("bubbles", "seek esrc");
    if (std::fread(bs, 4, m, f) != m) die("bubbles", "read esrc");
    if (std::fseek(f, dpos + (long)(done * 4), SEEK_SET) != 0) die("bubbles", "seek edst");
    if (std::fread(bd, 4, m, f) != m) die("bubbles", "read edst");
}

static Graph load_graph(const std::string& path) {
    FILE* f = xopen(path.c_str(), "rb");
    Graph g;
    uint32_t S = 0;
    if (std::fread(&S, 4, 1, f) != 1 || S < 2 || S > 64) die("bubbles", "graph.bin sample count out of range");
    g.S = (int)S;
    if (std::fread(&g.K, 8, 1, f) != 1) die("bubbles", "read K");
    g.sig.resize(g.K);
    if (std::fread(g.sig.data(), 8, g.K, f) != g.K) die("bubbles", "read sigs");

    g.som.assign(g.K, 0);
    constexpr size_t CH = 1u << 20;
    std::vector<cov_t> buf(CH);
    for (uint64_t done = 0; done < g.K;) {
        const size_t m = (size_t)std::min<uint64_t>(CH, g.K - done);
        if (std::fread(buf.data(), sizeof(cov_t), m, f) != m) die("bubbles", "read normal plane");
        for (size_t i = 0; i < m; i++) g.som[done + i] = (buf[i] == 0);
        done += m;
    }
    if (std::fseek(f, (long)((size_t)(S - 1) * g.K * sizeof(cov_t)), SEEK_CUR) != 0)
        die("bubbles", "seek past tumour planes");

    uint64_t E = 0;
    if (std::fread(&E, 8, 1, f) != 1) die("bubbles", "read E");
    if (E > 0xffffffffull) die("bubbles", "graph.bin has more than 2^32 edges");
    const long spos = std::ftell(f), dpos = spos + (long)(E * 4);
    std::vector<uint32_t> od(g.K + 1, 0), idg(g.K + 1, 0);
    std::vector<uint32_t> bs(CH), bd(CH);
    for (uint64_t done = 0; done < E;) {
        const size_t m = (size_t)std::min<uint64_t>(CH, E - done);
        edge_chunk(f, spos, dpos, done, m, bs.data(), bd.data());
        for (size_t i = 0; i < m; i++) { od[bs[i]]++; idg[bd[i]]++; }
        done += m;
    }
    g.oo.assign(g.K + 1, 0);
    g.io.assign(g.K + 1, 0);
    for (uint64_t i = 0; i < g.K; i++) { g.oo[i + 1] = g.oo[i] + od[i]; g.io[i + 1] = g.io[i] + idg[i]; }
    std::vector<uint32_t>().swap(od);
    std::vector<uint32_t>().swap(idg);
    g.on.resize(g.oo[g.K]);
    g.in.resize(g.io[g.K]);
    std::vector<uint32_t> c1(g.oo.begin(), g.oo.end()), c2(g.io.begin(), g.io.end());
    for (uint64_t done = 0; done < E;) {
        const size_t m = (size_t)std::min<uint64_t>(CH, E - done);
        edge_chunk(f, spos, dpos, done, m, bs.data(), bd.data());
        for (size_t i = 0; i < m; i++) { g.on[c1[bs[i]]++] = bd[i]; g.in[c2[bd[i]]++] = bs[i]; }
        done += m;
    }
    std::fclose(f);
    return g;
}

static uint32_t trace_arm(const Graph& g, uint32_t start, int maxlen, std::vector<uint32_t>& arm) {
    arm.clear();
    uint32_t cur = start;
    for (int s = 0; s <= maxlen; s++) {
        arm.push_back(cur);
        if (g.oo[cur + 1] - g.oo[cur] != 1) return UINT32_MAX;
        const uint32_t nxt = g.on[g.oo[cur]];
        const uint64_t din = g.io[nxt + 1] - g.io[nxt];
        if (din >= 2) return nxt;
        if (din != 1) return UINT32_MAX;
        cur = nxt;
    }
    return UINT32_MAX;
}

static bool uniform(const Graph& g, const std::vector<uint32_t>& arm, uint8_t want) {
    for (uint32_t v : arm) if (g.som[v] != want) return false;
    return true;
}

struct BubJob {
    const Graph* g;
    int t, nt, maxlen;
    std::vector<std::vector<Bub>>* part;
    uint64_t* nsrc;
};

static void bub_worker(BubJob* j) {
    const Graph& g = *j->g;
    const uint64_t lo = g.K * (uint64_t)j->t / (uint64_t)j->nt, hi = g.K * (uint64_t)(j->t + 1) / (uint64_t)j->nt;
    std::vector<uint32_t> a1, a2;
    for (uint64_t s = lo; s < hi; s++) {
        if (g.oo[s + 1] - g.oo[s] < 2) continue;
        (*j->nsrc)++;
        for (uint64_t i = g.oo[s]; i < g.oo[s + 1]; i++) {
            const uint32_t T1 = trace_arm(g, g.on[i], j->maxlen, a1);
            if (T1 == UINT32_MAX) continue;
            for (uint64_t jj = i + 1; jj < g.oo[s + 1]; jj++) {
                if (trace_arm(g, g.on[jj], j->maxlen, a2) != T1) continue;
                if (g.som[a1[0]] == g.som[a2[0]]) continue;
                if (!uniform(g, a1, g.som[a1[0]]) || !uniform(g, a2, g.som[a2[0]])) continue;
                (*j->part)[j->t].push_back(g.som[a1[0]] ? Bub{a1, a2} : Bub{a2, a1});
            }
        }
    }
}

static std::vector<Bub> find_bubbles(const Graph& g, int maxlen, int NT) {
    const int nt = std::max(1, NT);
    std::vector<std::vector<Bub>> part(nt);
    std::vector<uint64_t> nsrc(nt, 0);
    std::vector<BubJob> jobs(nt);
    for (int t = 0; t < nt; t++)
        jobs[t] = BubJob{&g, t, nt, maxlen, &part, &nsrc[t]};
    std::vector<std::thread> th;
    for (int t = 0; t < nt; t++) th.emplace_back(bub_worker, &jobs[t]);
    for (auto& t : th) t.join();

    std::vector<Bub> out;
    uint64_t total = 0, srcs = 0;
    for (int t = 0; t < nt; t++) { total += part[t].size(); srcs += nsrc[t]; }
    out.reserve(total);
    for (int t = 0; t < nt; t++) {
        for (Bub& b : part[t]) out.push_back(std::move(b));
        std::vector<Bub>().swap(part[t]);
    }
    say("bubbles", "%llu nodes, %llu branching sources -> %zu clean bubbles (maxlen=%d)",
        (unsigned long long)g.K, (unsigned long long)srcs, out.size(), maxlen);
    return out;
}

static bool arm_seq(const Index& ix, std::vector<SvbFile>& svb, const std::vector<uint32_t>& arm,
                    std::string& out) {
    static thread_local std::string core, scratch;
    out.clear();
    for (size_t i = 0; i < arm.size(); i++) {
        if (!core_bytes(ix, svb, arm[i], core, scratch)) return false;
        if (core.empty()) return false;
        if (i == 0) { out = core; continue; }
        size_t k = std::min(out.size(), core.size());
        while (k && out.compare(out.size() - k, k, core, 0, k) != 0) k--;
        if (!k) return false;
        out.append(core, k, core.size() - k);
    }
    return !out.empty();
}

static void diff_span(const std::string& t, const std::string& n, size_t& a, size_t& dt, size_t& dn) {
    const size_t m = std::min(t.size(), n.size());
    a = 0;
    while (a < m && t[a] == n[a]) a++;
    size_t b = 0;
    while (b < m - a && t[t.size() - 1 - b] == n[n.size() - 1 - b]) b++;
    dt = t.size() - a - b;
    dn = n.size() - a - b;
}

int cmd_bubbles(int argc, const char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: snape bubbles <graph.bin> <sheet.tsv> [threads=16] [maxlen=6] "
                             "[min_diff=25]\n");
        return 2;
    }
    const std::string gpath = argv[1], spath = argv[2];
    const int NT = argc > 3 ? std::atoi(argv[3]) : 16;
    const int MAXLEN = argc > 4 ? std::atoi(argv[4]) : 6;
    const int MIN_DIFF = argc > 5 ? std::atoi(argv[5]) : 25;

    const std::vector<Sample> smp = samples_read(spath);
    const Graph g = load_graph(gpath);
    if ((size_t)g.S != smp.size()) die("bubbles", "sheet and graph.bin disagree on the sample count");

    const std::string dir = dir_of(gpath);
    Index ix;
    ix.load(dir + "/index.bin");
    if (ix.K != g.K) die("bubbles", "index.bin and graph.bin disagree on the node count");
    std::vector<SvbFile> svb(smp.size());
    for (size_t i = 0; i < smp.size(); i++) svb[i].open(smp[i].path);
    say("bubbles", "index.bin: %llu records, N=%u", (unsigned long long)ix.off[ix.K], ix.N);
    const std::vector<Bub> bubs = find_bubbles(g, MAXLEN, NT);
    if (bubs.empty()) say("bubbles", "no clean bubbles");

    ankerl::unordered_dense::set<uint64_t> drop;
    uint64_t nbub = 0, nomiss = 0;
    std::string ts, gs;
    for (const Bub& b : bubs) {
        const bool sized = arm_seq(ix, svb, b.som_arm, ts) && arm_seq(ix, svb, b.germ_arm, gs);
        size_t a = 0, dt = 0, dn = 0;
        if (sized) diff_span(ts, gs, a, dt, dn);
        nomiss += !sized;
        bool dropped = sized && dt < (size_t)MIN_DIFF && dn < (size_t)MIN_DIFF;

        if (dropped) {
            for (uint32_t v : b.som_arm) drop.insert(g.sig[v]);
            nbub++;
        }
    }

    for (auto& v : svb) v.close();
    const std::string out = dir + "/drop.sigs";
    FILE* f = xopen(out.c_str(), "w");
    for (uint64_t s : drop) std::fprintf(f, "%llu\n", (unsigned long long)s);
    std::fclose(f);
    say("bubbles", "both arms differ by <%d bp: dropped %llu/%zu bubbles, %zu cores "
                   "(%llu unsized bubbles kept) -> %s",
        MIN_DIFF, (unsigned long long)nbub, bubs.size(), drop.size(),
        (unsigned long long)nomiss, out.c_str());
    return 0;
}
