// snape bubbles <graph.bin> <normal> <tumour> [threads=16] [k=21] [w=22] [maxlen=40] [max_dbp=50] [max_unmatched=20]
//
// Most of the somatic material in the graph is not a junction. It is SNVs and sub-SV edits, and each
// one still becomes a chunk, a contig, and a minute of alignment. This stage names those cores before
// any contig exists, so the aligner only ever sees candidates that could be an SV.
//
// A BUBBLE is a source S with out-degree >= 2, two node-disjoint linear arms, and a sink T where they
// converge -- one arm somatic (ncov == 0), the other germline. Anything else is not a bubble and is
// left alone.
//
//         S --[ somatic arm: c1 c2 ... ]--> T
//          \--[ germline arm: g1 g2 ... ]--^
//
// An L5 core spans ~217 bp, so a core signature cannot see a 50 bp change: two cores match or they do
// not. Minimizers (k=21 w=22, ~13 per core) open the arms up, and the comparison is ORDERED -- a set
// intersection measures nothing here, because a junction core is built from S's tail and T's head and
// both of those are present in the germline arm anyway. What separates the two cases is colinearity:
//
//   SNV / small edit : same sequence with a local change. The minimizer sequences agree except for a
//                      short run in the middle -> the arms are the same length and the longest common
//                      subsequence covers nearly all of the somatic arm.
//   junction         : different sequence. Minimizers may still be shared, but they do not stay in
//                      order, and the LCS collapses.
//
// Two clauses, and they measure orthogonal things: HOW MUCH sequence changed, and whether what is
// there is the same sequence rearranged. Both are needed -- a balanced junction changes no length at
// all, and only the collapsing LCS separates it from an SNV.
//
// The length clause is in BASES, not in minimizers, and that is the whole of it. A minimizer is one
// sample per ~16.6 raw bases (measured, hcc1954: median arm 31 minimizers over 508 bases), so the old
// `max_dlen = 12` minimizers admitted everything up to ~200 bp into the "local edit, delete it" class
// -- four times the 50 bp the caller reports from. That is where the insertion and deletion recall
// went: over four CASTLE panels the gap to Severus is +33 recall at 50-100 bp and ~0 by 500 bp.
//
// An arm's length in bases is the sum of its cores' lengths minus the overlap consecutive cores of one
// read share, which is 122.5 bp at level 5 -- measured over 1.97e9 consecutive pairs on hcc1954 and
// 2.02e9 on h2009, identical to a decimal, because it is a property of the decomposition and not of
// the sample. So `max_dbp` is comparable with an SV length, and the value to give it is not a free
// parameter any more: it is the caller's own floor, `min_sv`. Below that a difference is not an event
// we would report even if we kept it.
//
// Measured against the truth loci (hcc1954, 188 INS/DEL loci carrying a somatic core), at equal
// recovery of the events the shipped pipeline misses, the base rule keeps far fewer cores than the
// minimizer rule -- 20/22 missed deletions for 32.9k cores against 39.6k, 21/22 for 42.6k against
// 57.0k, 22/22 for 53.2k against 78.2k. Same recall, a quarter to a third less to align.
#include "bubbles.h"
#include "samples.h"

#include "snape.h"
#include "unordered_dense.h"

extern "C" {
#include "encoding.h"
}

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

namespace {

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[bubbles] FATAL: %s\n", m); std::exit(1); }

int LEVEL = 5;                                 // runtime, from SNAPE_LEVEL
constexpr int MARGIN = 10;                     // MUST match the graph build
// Consecutive cores of one read overlap, so an arm of n cores spans sum(len) - overlap*(n-1) bases.
// MEASURED on this run's own reads rather than taken as a constant -- `collect` counts it over every
// consecutive pair it parses (1.97e9 pairs on hcc1954, 2.02e9 on h2009; both give 122.5 bp, because it
// is a property of the level-5 decomposition and not of the sample). The fallback is only for a run
// that parsed nothing.
double OVERLAP = 122.5;

// SNAPE_BUBBLE_TSV=<path>: dump one line per bubble with every statistic the drop rule could be
// built from, so a rule can be calibrated offline against the truth loci instead of by re-running
// this stage once per candidate threshold. Off by default and it changes nothing about drop.sigs.
const char* DUMP = nullptr;

inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return x;
}

inline int b2(char c) {
    switch (c) { case 'A': case 'a': return 0; case 'C': case 'c': return 1;
                 case 'G': case 'g': return 2; case 'T': case 't': return 3; }
    return -1;
}

// Canonical (strand-independent) minimizers, one per w-window, duplicates of the previous pick
// dropped so the result is the window's sequence of choices rather than one entry per position.
void minimizers(const char* s, int len, int k, int w, std::vector<uint64_t>& out) {
    out.clear();
    if (len < k) return;
    const uint64_t mask = (k >= 32) ? ~0ull : ((1ull << (2 * k)) - 1);
    uint64_t fw = 0, rv = 0; int have = 0;
    static thread_local std::vector<uint64_t> h;
    h.clear(); h.reserve((size_t)len);
    for (int i = 0; i < len; i++) {
        const int c = b2(s[i]);
        if (c < 0) { have = 0; fw = rv = 0; h.push_back(~0ull); continue; }
        fw = ((fw << 2) | (uint64_t)c) & mask;
        rv = (rv >> 2) | ((uint64_t)(3 - c) << (2 * (k - 1)));
        if (++have >= k) h.push_back(mix64(fw < rv ? fw : rv)); else h.push_back(~0ull);
    }
    uint64_t last = ~0ull;
    for (size_t i = 0; i + (size_t)w <= h.size(); i++) {
        uint64_t best = ~0ull;
        for (int j = 0; j < w; j++) best = std::min(best, h[i + j]);
        if (best != ~0ull && best != last) { out.push_back(best); last = best; }
    }
}

// Homopolymer-compressed, so a HiFi run-length wobble does not read as a sequence difference.
void hpc_local(const char* s, int len, std::string& out) {
    out.clear();
    for (int i = 0; i < len; i++) if (i == 0 || s[i] != s[i - 1]) out.push_back(s[i]);
}

// The cores this stage needs sequence for, each given a dense id so a hit costs one lookup and no
// allocation. Membership is read-only once built, so every worker can probe it at once.
struct Want {
    ankerl::unordered_dense::map<uint64_t, uint32_t> id;
    uint32_t add(uint64_t sig) {
        const uint32_t n = (uint32_t)id.size();
        return id.emplace(sig, n).first->second;
    }
    uint32_t find(uint64_t sig) const {
        const auto it = id.find(sig);
        return it == id.end() ? UINT32_MAX : it->second;
    }
    size_t size() const { return id.size(); }
};

size_t lcs(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    if (a.empty() || b.empty()) return 0;
    std::vector<uint16_t> prev(b.size() + 1, 0), cur(b.size() + 1, 0);
    for (size_t i = 1; i <= a.size(); i++) {
        for (size_t j = 1; j <= b.size(); j++)
            cur[j] = (a[i - 1] == b[j - 1]) ? (uint16_t)(prev[j - 1] + 1) : std::max(prev[j], cur[j - 1]);
        std::swap(prev, cur);
    }
    return prev[b.size()];
}

// ---- graph ----
struct Graph {
    uint64_t K = 0;
    std::vector<uint64_t> sig;
    std::vector<uint8_t> som;
    std::vector<uint64_t> oo, io;               // directed CSR: a bubble is a DIRECTED structure
    std::vector<uint32_t> on, in;
    int S = 0;
    std::vector<uint32_t> tc;                   // per-node tumour count, kept only under SNAPE_BUBBLE_TSV
    std::vector<uint32_t> slen;                 // per-node RAW core length, SEQ_RAW = stored verbatim
    std::vector<uint64_t> soff;                 // where that core's bytes start in `seq`
    std::vector<uint8_t> seq;                   // every surviving core's bytes, 2-bit packed
    uint64_t ov_s = 0, ov_c = 0;                // summed consecutive-core overlap, and its pair count
};

// A core's bytes, back from the 2-bit packing graph wrote.
void seq_get(const Graph& g, uint32_t v, std::string& out) {
    const uint32_t len = g.slen[v], n = len & ~SEQ_RAW;
    out.resize(n);
    if (!n) return;
    const uint8_t* p = g.seq.data() + g.soff[v];
    if (len & SEQ_RAW) { std::memcpy(&out[0], p, n); return; }
    static const char L[4] = { 'A', 'C', 'G', 'T' };
    for (uint32_t i = 0; i < n; i++) out[i] = L[(p[i >> 2] >> (6 - 2 * (i & 3))) & 3];
}

// graph.bin v2, streamed: nothing is materialised that the CSR replaces. The count planes are
// SoA now, so only plane 0 (the normal) is read at all -- `som` is exactly "the normal never saw
// this core", and survival already guarantees some sample did. The edge endpoints are two separate
// SoA columns, so each chunk costs two seeks instead of holding a 1.2 GB AoS array.
Graph load_graph(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb"); if (!f) die("open graph.bin");
    Graph g;
    char mg[4]; uint32_t ver = 0, S = 0, resv = 0;
    if (std::fread(mg, 1, 4, f) != 4 || std::memcmp(mg, "SNPG", 4)) die("not a graph.bin (bad magic)");
    if (std::fread(&ver, 4, 1, f) != 1 || ver != 4) die("graph.bin version is not 4 -- rebuild it");
    if (std::fread(&S, 4, 1, f) != 1 || S < 2 || S > 64) die("graph.bin sample count out of range");
    if (std::fread(&resv, 4, 1, f) != 1) die("read header");
    if (std::fread(&g.ov_s, 8, 1, f) != 1 || std::fread(&g.ov_c, 8, 1, f) != 1) die("read overlap");
    g.S = (int)S;
    if (std::fread(&g.K, 8, 1, f) != 1) die("read K");
    g.sig.resize(g.K);
    if (std::fread(g.sig.data(), 8, g.K, f) != g.K) die("read sigs");

    g.som.assign(g.K, 0);
    constexpr size_t CH = 1u << 20;
    std::vector<cov_t> buf(CH);
    for (uint64_t done = 0; done < g.K; ) {                       // plane 0 = normal: somatic iff 0
        const size_t m = (size_t)std::min<uint64_t>(CH, g.K - done);
        if (std::fread(buf.data(), sizeof(cov_t), m, f) != m) die("read normal plane");
        for (size_t i = 0; i < m; i++) g.som[done + i] = (buf[i] == 0);
        done += m;
    }
    if (DUMP) {                                                   // diagnostic only: max over lines
        g.tc.assign(g.K, 0);
        for (uint32_t j = 1; j < S; j++)
            for (uint64_t done = 0; done < g.K; ) {
                const size_t m = (size_t)std::min<uint64_t>(CH, g.K - done);
                if (std::fread(buf.data(), sizeof(cov_t), m, f) != m) die("read tumour plane");
                for (size_t i = 0; i < m; i++) g.tc[done + i] = std::max(g.tc[done + i], (uint32_t)buf[i]);
                done += m;
            }
    } else if (std::fseek(f, (long)((size_t)(S - 1) * g.K * sizeof(cov_t)), SEEK_CUR) != 0) {
        die("seek past tumour planes");
    }

    uint64_t sbytes = 0;
    if (std::fread(&sbytes, 8, 1, f) != 1) die("read core sequence size");
    g.slen.resize(g.K);
    if (g.K && std::fread(g.slen.data(), 4, g.K, f) != g.K) die("read core lengths");
    g.seq.resize(sbytes);
    if (sbytes && std::fread(g.seq.data(), 1, sbytes, f) != sbytes) die("read core sequence");
    g.soff.resize(g.K);
    { uint64_t at = 0;
      for (uint64_t i = 0; i < g.K; i++) { g.soff[i] = at; at += seq_bytes(g.slen[i]); }
      if (at != sbytes) die("core sequence section is the wrong size"); }

    uint64_t E = 0;
    if (std::fread(&E, 8, 1, f) != 1) die("read E");
    // The fill cursors below are 32-bit, which halves what they cost on a 75 M-edge graph. Node ids
    // are 32-bit anyway, so a graph this loader could not take is one nothing else could read
    // either -- but say so rather than truncate.
    if (E > 0xffffffffull) die("graph.bin has more than 2^32 edges");
    const long spos = std::ftell(f), dpos = spos + (long)(E * 4);
    std::vector<uint32_t> od(g.K + 1, 0), idg(g.K + 1, 0);
    std::vector<uint32_t> bs(CH), bd(CH);
    auto edge_chunk = [&](uint64_t done, size_t m) {
        if (std::fseek(f, spos + (long)(done * 4), SEEK_SET) != 0) die("seek esrc");
        if (std::fread(bs.data(), 4, m, f) != m) die("read esrc");
        if (std::fseek(f, dpos + (long)(done * 4), SEEK_SET) != 0) die("seek edst");
        if (std::fread(bd.data(), 4, m, f) != m) die("read edst");
    };
    for (uint64_t done = 0; done < E; ) {                         // pass 1: degrees
        const size_t m = (size_t)std::min<uint64_t>(CH, E - done);
        edge_chunk(done, m);
        for (size_t i = 0; i < m; i++) { od[bs[i]]++; idg[bd[i]]++; }
        done += m;
    }
    g.oo.assign(g.K + 1, 0); g.io.assign(g.K + 1, 0);
    for (uint64_t i = 0; i < g.K; i++) { g.oo[i + 1] = g.oo[i] + od[i]; g.io[i + 1] = g.io[i] + idg[i]; }
    std::vector<uint32_t>().swap(od); std::vector<uint32_t>().swap(idg);
    g.on.resize(g.oo[g.K]); g.in.resize(g.io[g.K]);
    std::vector<uint32_t> c1(g.oo.begin(), g.oo.end()), c2(g.io.begin(), g.io.end());
    for (uint64_t done = 0; done < E; ) {                         // pass 2: fill
        const size_t m = (size_t)std::min<uint64_t>(CH, E - done);
        edge_chunk(done, m);
        for (size_t i = 0; i < m; i++) { g.on[c1[bs[i]]++] = bd[i]; g.in[c2[bd[i]]++] = bs[i]; }
        done += m;
    }
    std::fclose(f);
    return g;
}

// Follow one linear arm until it converges on a node with in-degree >= 2 (the sink) or stops being
// linear. Returns the sink, or UINT32_MAX when this is not an arm.
uint32_t trace_arm(const Graph& g, uint32_t start, int maxlen, std::vector<uint32_t>& arm) {
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

bool uniform(const Graph& g, const std::vector<uint32_t>& arm, uint8_t want) {
    for (uint32_t v : arm) if (g.som[v] != want) return false;
    return true;
}

struct Bub { std::vector<uint32_t> som_arm, germ_arm; };

// Every branching node is independent of every other one, so the sweep is split into contiguous
// ranges and the ranges are concatenated in range order -- the same bubbles in the same order as the
// single-threaded sweep produced, which is what keeps drop.sigs byte-for-byte the same.
std::vector<Bub> find_bubbles(const Graph& g, int maxlen, int NT) {
    const int nt = std::max(1, NT);
    std::vector<std::vector<Bub>> part(nt);
    std::vector<uint64_t> nsrc(nt, 0);
    std::vector<std::thread> th;
    for (int t = 0; t < nt; t++) th.emplace_back([&, t]{
        const uint64_t lo = g.K * t / nt, hi = g.K * (t + 1) / nt;
        std::vector<uint32_t> a1, a2;
        for (uint64_t s = lo; s < hi; s++) {
            if (g.oo[s + 1] - g.oo[s] < 2) continue;
            nsrc[t]++;
            for (uint64_t i = g.oo[s]; i < g.oo[s + 1]; i++) {
                const uint32_t T1 = trace_arm(g, g.on[i], maxlen, a1);
                if (T1 == UINT32_MAX) continue;
                for (uint64_t j = i + 1; j < g.oo[s + 1]; j++) {
                    if (trace_arm(g, g.on[j], maxlen, a2) != T1) continue;
                    if (g.som[a1[0]] == g.som[a2[0]]) continue;      // need exactly one somatic arm
                    if (!uniform(g, a1, g.som[a1[0]]) || !uniform(g, a2, g.som[a2[0]])) continue;
                    part[t].push_back(g.som[a1[0]] ? Bub{a1, a2} : Bub{a2, a1});
                }
            }
        }
    });
    for (auto& t : th) t.join();
    std::vector<Bub> out;
    uint64_t total = 0, srcs = 0;
    for (int t = 0; t < nt; t++) { total += part[t].size(); srcs += nsrc[t]; }
    out.reserve(total);
    for (int t = 0; t < nt; t++) {
        for (Bub& b : part[t]) out.push_back(std::move(b));
        std::vector<Bub>().swap(part[t]);
    }
    std::fprintf(stderr, "[bubbles] %llu dugum, %llu dallanma -> %zu temiz bubble (maxlen=%d)\n",
                 (unsigned long long)g.K, (unsigned long long)srcs, out.size(), maxlen);
    return out;
}

// Unpack each wanted core and take its homopolymer-compressed length and its minimizers. ~225k
// cores on the 21-sample graph, so one thread is ample -- this replaces a 2.4 TB read.
void derive(const Graph& g, const std::vector<uint32_t>& nod, int k, int w,
            std::vector<std::vector<uint64_t>>& mins_of, std::vector<uint32_t>& raw_of,
            std::vector<uint32_t>& hpc_of) {
    const size_t W = nod.size();
    mins_of.assign(W, {}); raw_of.assign(W, 0); hpc_of.assign(W, 0);
    std::string core, hs; std::vector<uint64_t> mins;
    for (size_t i = 0; i < W; i++) {
        seq_get(g, nod[i], core);
        if (core.empty()) continue;
        hpc_local(core.data(), (int)core.size(), hs);
        minimizers(hs.data(), (int)hs.size(), k, w, mins);
        mins_of[i] = mins;
        raw_of[i] = (uint32_t)core.size();
        hpc_of[i] = (uint32_t)hs.size();
    }
}

void arm_mins(const Graph& g, const std::vector<uint32_t>& arm, const Want& want,
              const std::vector<std::vector<uint64_t>>& mins, std::vector<uint64_t>& out) {
    out.clear();
    for (uint32_t v : arm) {
        const uint32_t id = want.find(g.sig[v]);
        if (id != UINT32_MAX) out.insert(out.end(), mins[id].begin(), mins[id].end());
    }
}

}  // namespace

int cmd_bubbles(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: snape bubbles <graph.bin> <sheet.tsv> [threads=16] [k=21] "
                             "[w=22] [maxlen=40] [max_dbp=50] [max_unmatched=20]\n");
        return 2;
    }
    const std::string gpath = argv[1], spath = argv[2];
    const int NT = argc > 3 ? std::atoi(argv[3]) : 16;
    const int KMER = argc > 4 ? std::atoi(argv[4]) : 21;
    const int W = argc > 5 ? std::atoi(argv[5]) : 22;
    const int MAXLEN = argc > 6 ? std::atoi(argv[6]) : 40;
    // 50 and not a swept number: it is `snape call`'s own min_sv. Two arms that differ by
    // less than the smallest event we would report differ by nothing we would report.
    const int MAX_DBP = argc > 7 ? std::atoi(argv[7]) : 50;
    const int MAX_UNM = argc > 8 ? std::atoi(argv[8]) : 20;
    if (const char* lv = std::getenv("SNAPE_LEVEL")) LEVEL = std::atoi(lv);
    DUMP = std::getenv("SNAPE_BUBBLE_TSV");

    const std::vector<Sample> smp = samples_read(spath);
    const Graph g = load_graph(gpath);
    if ((size_t)g.S != smp.size()) die("sheet and graph.bin disagree on the sample count");
    const std::vector<Bub> bubs = find_bubbles(g, MAXLEN, NT);
    if (bubs.empty()) std::fprintf(stderr, "[bubbles] temiz bubble yok\n");

    // The arms' bytes come from graph.bin, which carries every surviving core's raw sequence. This
    // used to be a scan of every read set to recover them: measured on the 21-sample simulation,
    // 2.4 TB of reading to obtain 49 MB of sequence. The signature IS the hash of a core's
    // homopolymer-compressed bytes, so whichever read first gave the graph this core, its
    // compression and its minimizers are the ones any other read would have produced.
    Want want_s, want_g;
    std::vector<uint32_t> nod_s, nod_g;
    auto addn = [](Want& w, std::vector<uint32_t>& nod, uint64_t sig, uint32_t v) {
        const uint32_t i = w.add(sig);
        if (i >= nod.size()) nod.resize(i + 1);
        nod[i] = v;
    };
    for (const Bub& b : bubs) {
        for (uint32_t v : b.som_arm) addn(want_s, nod_s, g.sig[v], v);
        for (uint32_t v : b.germ_arm) addn(want_g, nod_g, g.sig[v], v);
    }
    std::vector<std::vector<uint64_t>> mm_s, mm_g;
    std::vector<uint32_t> raw_s, hpc_s, raw_g, hpc_g;
    derive(g, nod_s, KMER, W, mm_s, raw_s, hpc_s);
    derive(g, nod_g, KMER, W, mm_g, raw_g, hpc_g);
    // graph measured the overlap over every read of every sample as it went; the mean is what makes
    // two arms of different core count comparable in bases.
    OVERLAP = g.ov_c ? (double)g.ov_s / (double)g.ov_c : OVERLAP;
    std::fprintf(stderr, "[bubbles] %zu + %zu cekirdek grafdan cozuldu, ortusme %.1f bp (%llu cift)\n",
                 nod_s.size(), nod_g.size(), OVERLAP, (unsigned long long)g.ov_c);

    // Both arms the same length and the somatic arm almost entirely colinear with the germline one:
    // a local edit, not a novel adjacency. Its cores never become a contig.
    FILE* dmp = nullptr;
    if (DUMP) {
        dmp = std::fopen(DUMP, "w");
        if (!dmp) die("open SNAPE_BUBBLE_TSV");
        std::fprintf(dmp, "ns\tng\tms\tmg\tlcs\tsraw\tgraw\tshpc\tghpc\tsmiss\tgmiss\ttmin\ttmax\tdrop\tsigs\n");
    }
    // Sum an arm's core lengths, and count the cores no read gave us bytes for.
    auto arm_len = [&](const std::vector<uint32_t>& arm, const Want& want,
                       const std::vector<uint32_t>& raw, const std::vector<uint32_t>& hpc,
                       uint64_t& sraw, uint64_t& shpc, uint64_t& miss) {
        sraw = shpc = miss = 0;
        for (uint32_t v : arm) {
            const uint32_t id = want.find(g.sig[v]);
            if (id == UINT32_MAX || raw[id] == 0) { miss++; continue; }
            sraw += raw[id]; shpc += hpc[id];
        }
    };

    ankerl::unordered_dense::set<uint64_t> drop;
    uint64_t nbub = 0, nomiss = 0;
    std::vector<uint64_t> ms, mg;
    for (const Bub& b : bubs) {
        arm_mins(g, b.som_arm, want_s, mm_s, ms);
        arm_mins(g, b.germ_arm, want_g, mm_g, mg);
        uint64_t sraw, shpc, smiss, graw, ghpc, gmiss;
        arm_len(b.som_arm, want_s, raw_s, hpc_s, sraw, shpc, smiss);
        arm_len(b.germ_arm, want_g, raw_g, hpc_g, graw, ghpc, gmiss);
        // A core no read gave us bytes for has no length, so this arm's span is not measurable and
        // the bubble is kept rather than judged on a number that is short by a whole core.
        const bool sized = !ms.empty() && !mg.empty() && ms.size() <= 400 && mg.size() <= 400 &&
                           smiss == 0 && gmiss == 0;
        nomiss += (smiss || gmiss);
        const double dspan = std::fabs(((double)sraw - (double)graw) -
                                       OVERLAP * ((double)b.som_arm.size() - (double)b.germ_arm.size()));
        const bool dropped = sized && dspan <= (double)MAX_DBP &&
            (long long)(ms.size() - lcs(ms, mg)) <= MAX_UNM;
        if (dropped) {
            for (uint32_t v : b.som_arm) drop.insert(g.sig[v]);
            nbub++;
        }
        if (!dmp) continue;
        uint32_t tmin = UINT32_MAX, tmax = 0;
        for (uint32_t v : b.som_arm) { tmin = std::min(tmin, g.tc[v]); tmax = std::max(tmax, g.tc[v]); }
        std::fprintf(dmp, "%zu\t%zu\t%zu\t%zu\t%zu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%u\t%u\t%d\t",
                     b.som_arm.size(), b.germ_arm.size(), ms.size(), mg.size(),
                     sized ? lcs(ms, mg) : (size_t)0, (unsigned long long)sraw, (unsigned long long)graw,
                     (unsigned long long)shpc, (unsigned long long)ghpc, (unsigned long long)smiss,
                     (unsigned long long)gmiss, tmin == UINT32_MAX ? 0 : tmin, tmax, (int)dropped);
        for (size_t i = 0; i < b.som_arm.size(); i++)
            std::fprintf(dmp, "%s%llu", i ? "," : "", (unsigned long long)g.sig[b.som_arm[i]]);
        std::fputc('\n', dmp);
    }
    if (dmp) { std::fclose(dmp); std::fprintf(stderr, "[bubbles] bubble istatistikleri -> %s\n", DUMP); }

    const size_t cut = gpath.find_last_of('/');
    const std::string out = (cut == std::string::npos ? std::string(".") : gpath.substr(0, cut)) + "/drop.sigs";
    FILE* f = std::fopen(out.c_str(), "w"); if (!f) die("open drop.sigs");
    for (uint64_t s : drop) std::fprintf(f, "%llu\n", (unsigned long long)s);
    std::fclose(f);
    std::fprintf(stderr, "[bubbles] |dspan|<=%dbp & eslesmeyen<=%d -> %llu/%zu bubble, %zu cekirdek "
                 "(%llu bubble olculemedi, tutuldu) -> %s\n",
                 MAX_DBP, MAX_UNM, (unsigned long long)nbub, bubs.size(), drop.size(),
                 (unsigned long long)nomiss, out.c_str());
    return 0;
}
