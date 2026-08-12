// snape junctions <graph.bin> <reads> <out.fa> [threads=16] [pad_cores=50] [per_chunk=5] [min_flank=300] [min_vaf=0]
//
// Junction contigs instead of whole reads: only the sequence around each somatic run, and only from
// as many reads as the breakpoint caller actually needs.
//
// A CHUNK is a connected component of the graph's somatic subgraph (ncov == 0, tcov > 0) -- one
// junction candidate, found without a reference. The reads are then scanned once and, for each chunk,
// the first `per_chunk` reads that carry it are cut down to
//
//     [ start of the (first somatic core - pad_cores) , end of the (last somatic core + pad_cores) ]
//
// which is the novel sequence plus about a kilobase of normal anchor either side. No assembly and no
// overlap merging is needed to get there: consecutive LCP cores of ONE read are already adjacent and
// overlapping in that read's own bytes, so the span between two core boundaries IS the sequence. The
// only case that would need merging across reads is a chunk no single read spans with anchor to
// spare, and `min_flank` counts exactly those instead of assuming there are none.
//
// A record is named <read>_<n> so that the caller downstream still counts one supporting READ per
// record, while a read that carries two different junctions cannot inflate either one's support.
//
// `drop.sigs` beside the graph -- written by `snape bubbles` -- lists somatic cores to exclude before
// chunks are built: a core whose bubble has two arms of equal length that agree minimizer-for-
// minimizer is an SNV or a sub-SV edit, not a junction, and carving it out here means no contig is
// ever emitted for it. Absent file, no filter.
//
// `min_vaf` drops a chunk whose tumour support is too thin a fraction of the local depth. The
// denominator cannot come from the somatic subgraph -- a somatic core has ncov == 0 by definition, so
// tcov/(ncov+tcov) is identically 1 -- it comes from the germline cores the chunk hangs off, whose
// tcov counts every tumour read through that locus, variant-carrying or not. That ratio is the same
// quantity a caller's VAF measures, and unlike a flat tcov floor it is blind neither to a
// low-coverage region nor to a repeat whose counts are inflated.
#include "junctions.h"

#include "snape.h"
#include "unordered_dense.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {
int LEVEL = 5;                             // runtime, from SNAPE_LEVEL; echoed in the banner
constexpr int MARGIN = 10;                 // MUST match the graph build
// 256 reads a batch, 32 batches in flight: the whole point is that a worker never waits on the
// input, and 2000-read batches (34 MB each, one per worker) were both coarse and heavy.
constexpr size_t BATCH = 256, Q_CAP = 32;

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[jcontigs] FATAL: %s\n", m); std::exit(1); }

// ---- somatic subgraph -> chunks ----
struct Chunks {
    std::vector<uint64_t> sig;             // sorted somatic signatures
    std::vector<uint32_t> chunk;           // parallel: which component each belongs to
    size_t n = 0;                          // number of components
};

uint32_t findp(std::vector<uint32_t>& p, uint32_t x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

// `snape bubbles` output, if that stage has run: the somatic cores it identified as sub-SV edits
ankerl::unordered_dense::set<uint64_t> load_drop(const std::string& graph) {
    ankerl::unordered_dense::set<uint64_t> d;
    const size_t cut = graph.find_last_of('/');
    const std::string p = (cut == std::string::npos ? std::string(".") : graph.substr(0, cut)) + "/drop.sigs";
    FILE* f = std::fopen(p.c_str(), "r");
    if (!f) { std::fprintf(stderr, "[jcontigs] %s yok - bubble filtresi uygulanmiyor\n", p.c_str()); return d; }
    unsigned long long v;
    while (std::fscanf(f, "%llu", &v) == 1) d.insert((uint64_t)v);
    std::fclose(f);
    std::fprintf(stderr, "[jcontigs] %s: %zu cekirdek atlanacak\n", p.c_str(), d.size());
    return d;
}

Chunks build_chunks(const std::string& path, uint64_t& K, size_t& nsom, double min_vaf) {
    FILE* f = std::fopen(path.c_str(), "rb"); if (!f) die("open graph.bin");
    if (std::fread(&K, 8, 1, f) != 1) die("read K");
    std::vector<uint64_t> sig(K); std::vector<uint32_t> nc(K), tc(K);
    if (std::fread(sig.data(), 8, K, f) != K) die("read sigs");
    if (std::fread(nc.data(), 4, K, f) != K) die("read ncov");
    if (std::fread(tc.data(), 4, K, f) != K) die("read tcov");
    uint64_t E = 0;
    if (std::fread(&E, 8, 1, f) != 1) die("read E");
    std::vector<uint32_t> ed(E * 4);
    if (std::fread(ed.data(), 4, E * 4, f) != E * 4) die("read edges");
    std::fclose(f);

    const auto drop = load_drop(path);
    std::vector<uint8_t> som(K, 0);
    uint64_t ndrop = 0;
    for (uint64_t i = 0; i < K; i++) {
        som[i] = (nc[i] == 0 && tc[i] > 0);
        if (som[i] && drop.contains(sig[i])) { som[i] = 0; ndrop++; }
    }
    if (ndrop) std::fprintf(stderr, "[jcontigs] %llu somatik cekirdek drop.sigs ile elendi\n", (unsigned long long)ndrop);
    std::vector<uint32_t> par(K);
    for (uint64_t i = 0; i < K; i++) par[i] = (uint32_t)i;
    for (uint64_t e = 0; e < E; e++) {
        const uint32_t a = ed[e * 4], b = ed[e * 4 + 1];
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

    // Numerator: the best-supported somatic core in the chunk -- every read over the junction emits
    // most of the chunk's cores, so these counts sit close together and the max is the cleanest read
    // of "how many tumour reads saw this junction". Denominator: the highest-tcov germline neighbour,
    // which is the conservative choice (a bigger denominator only ever lowers the VAF).
    std::vector<uint32_t> num(NC, 0), den(NC, 0);
    for (uint64_t i = 0; i < K; i++) if (som[i]) num[cid[i]] = std::max(num[cid[i]], tc[i]);
    for (uint64_t e = 0; e < E; e++) {
        const uint32_t a = ed[e * 4], b = ed[e * 4 + 1];
        if (som[a] == som[b]) continue;
        const uint32_t s = som[a] ? a : b, o = som[a] ? b : a;
        den[cid[s]] = std::max(den[cid[s]], tc[o]);
    }
    // No germline neighbour at all means no denominator to judge by, so the chunk is kept rather than
    // guessed at; min_flank is what catches those anyway.
    auto vaf = [&](size_t c) { return den[c] ? (double)num[c] / den[c] : 1.0; };

    size_t noanch = 0;
    for (size_t c = 0; c < NC; c++) noanch += (den[c] == 0);
    std::fprintf(stderr, "[jcontigs] chunk VAF (chunk tcov / anchor tcov), %zu chunks, %zu with no germline anchor\n",
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
    for (uint64_t i = 0; i < K; i++) {
        if (!som[i] || keep[cid[i]] == UINT32_MAX) continue;
        c.sig.push_back(sig[i]); c.chunk.push_back(keep[cid[i]]);
    }
    nsom = c.sig.size();
    c.n = nk;
    // sig[] came out of graph.bin already sorted, so the parallel arrays are sorted with it
    return c;
}

int32_t chunk_of(const Chunks& c, uint64_t s) {
    const auto it = std::lower_bound(c.sig.begin(), c.sig.end(), s);
    if (it == c.sig.end() || *it != s) return -1;
    return (int32_t)c.chunk[(size_t)(it - c.sig.begin())];
}

// ---- reads ----
struct Read { std::string name, seq; };
}  // namespace

int cmd_junctions(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: snape junctions <graph.bin> <reads> <out.fa> "
                                 "[threads=16] [pad_cores=50] [per_chunk=5] [min_flank=300] [min_vaf=0]\n"); return 2; }
    const std::string gpath = argv[1], rpath = argv[2], opath = argv[3];
    const int NT = argc > 4 ? std::atoi(argv[4]) : 16;
    // 50 cores of anchor each side, not 10: half of what hcc1954 missed was a placement failure and
    // not a detection one, and a bigger pad cannot cost a contig -- the window is clamped to the
    // read's own cores and min_flank tests the flank ACHIEVED. Measured at 15/25/50/75; 50 dominates
    // 25 at every support level and 75 gives back what it gains.
    const int PAD = argc > 5 ? std::atoi(argv[5]) : 50;
    const uint32_t PER = argc > 6 ? (uint32_t)std::atoi(argv[6]) : 5;
    const int64_t MINF = argc > 7 ? std::atoll(argv[7]) : 300;
    const double MINVAF = argc > 8 ? std::atof(argv[8]) : 0.0;
    if (const char* lv = std::getenv("SNAPE_LEVEL")) LEVEL = std::atoi(lv);

    uint64_t K = 0; size_t nsom = 0;
    const Chunks C = build_chunks(gpath, K, nsom, MINVAF);
    std::fprintf(stderr, "[jcontigs] %s: %llu nodes, %zu somatic -> %zu chunks; pad=%d cores, "
                 "<=%u reads/chunk, min flank %lld bp, min VAF %.2f, level=%d, %d threads\n",
                 gpath.c_str(), (unsigned long long)K, nsom, C.n, PAD, PER, (long long)MINF, MINVAF, LEVEL, NT);
    if (!C.n) die("no somatic chunks");

    std::vector<std::atomic<uint32_t>> taken(C.n);
    for (auto& t : taken) t.store(0);
    std::vector<std::atomic<uint8_t>> touched(C.n); // a read carried it, whatever the flanks were
    for (auto& t : touched) t.store(0);

    Reads in = reads_open(rpath);
    FILE* out = std::fopen((opath + ".tmp").c_str(), "w"); if (!out) die("open out");
    std::mutex out_mu;
    std::atomic<uint64_t> seen{0}, emitted{0}, emit_bp{0}, seen_bp{0}, short_flank{0};

    // One reader fills batches, the workers take whole batches. The loop this replaces filled the
    // batch under the input mutex AND copied every read into it while holding that mutex, so with
    // 2000-read batches the lock was busy about as often as the workers were: measured 684% of the
    // 1600% it was given.
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<std::vector<Read>> q;
    bool eof = false;

    auto reader = [&]() {
        std::vector<Read> batch; batch.reserve(BATCH);
        Read r;
        auto push = [&]() {
            if (batch.empty()) return;
            std::unique_lock<std::mutex> lk(qmtx);
            cv_full.wait(lk, [&]{ return q.size() < Q_CAP; });
            q.push(std::move(batch)); cv_empty.notify_one();
            batch.clear(); batch.reserve(BATCH);
        };
        while (reads_next(in, r.seq, &r.name)) {
            batch.push_back(std::move(r));
            if (batch.size() >= BATCH) push();
        }
        push();
        { std::lock_guard<std::mutex> lk(qmtx); eof = true; }
        cv_empty.notify_all();
    };

    auto worker = [&]() {
        std::string chunk_out;
        std::vector<int> ci; std::vector<int32_t> cc;      // reused across reads
        for (;;) {
            std::vector<Read> batch;
            {
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || eof; });
                if (q.empty()) return;
                batch = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            chunk_out.clear();
            uint64_t sbp = 0, ebp = 0, en = 0, sf = 0;
            for (const Read& r : batch) {
                sbp += r.seq.size();
                int n; const uint64_t *st, *en_, *sg;
                extract(r.seq.data(), (int)r.seq.size(), 0, LEVEL, &n, &st, &en_, &sg);
                // the cores this read keeps, and which chunk each belongs to
                ci.clear(); cc.clear();
                for (int i = 0; i < n; i++) {
                    if (st[i] < (uint64_t)MARGIN || en_[i] + MARGIN > r.seq.size()) continue;
                    ci.push_back(i); cc.push_back(chunk_of(C, sg[i]));
                }
                int w = 0;
                for (size_t a = 0; a < ci.size();) {
                    if (cc[a] < 0) { a++; continue; }
                    size_t b = a;
                    // Maximal somatic run, joining runs closer than the padding would put them
                    // anyway. A mobile-element insertion is only somatic at its two JOINS -- its
                    // interior is an Alu or a L1 the reference already has -- so the two ends arrive
                    // as separate runs, and cutting them into separate windows means no contig ever
                    // contains the insertion. Merged, one window spans join to join.
                    for (;;) {
                        while (b + 1 < ci.size() && cc[b + 1] >= 0) b++;
                        size_t nxt = b + 1;
                        while (nxt < ci.size() && cc[nxt] < 0) nxt++;
                        if (nxt >= ci.size() || nxt - b > (size_t)(2 * PAD)) break;
                        b = nxt;
                    }
                    const uint32_t ch = (uint32_t)cc[a];
                    touched[ch].store(1, std::memory_order_relaxed);
                    const size_t la = a >= (size_t)PAD ? a - PAD : 0;
                    const size_t rb = std::min(ci.size() - 1, b + (size_t)PAD);
                    const int64_t ws = (int64_t)st[ci[la]], we = (int64_t)en_[ci[rb]];
                    const int64_t lf = (int64_t)st[ci[a]] - ws, rf = we - (int64_t)en_[ci[b]];
                    a = b + 1;
                    if (lf < MINF || rf < MINF) { sf++; continue; }
                    if (taken[ch].fetch_add(1) >= PER) continue;
                    char h[256];
                    std::snprintf(h, sizeof h, ">%s_%d chunk=%u len=%lld flank=%lld,%lld\n",
                                  r.name.c_str(), w++, ch, (long long)(we - ws), (long long)lf, (long long)rf);
                    chunk_out += h;
                    chunk_out.append(r.seq, (size_t)ws, (size_t)(we - ws));
                    chunk_out += "\n";
                    en++; ebp += (uint64_t)(we - ws);
                }
            }
            seen += batch.size(); seen_bp += sbp; emitted += en; emit_bp += ebp; short_flank += sf;
            {
                std::lock_guard<std::mutex> lk(out_mu);
                std::fwrite(chunk_out.data(), 1, chunk_out.size(), out);
            }
            const uint64_t s = seen.load();
            if (s / 1000000 != (s - batch.size()) / 1000000)
                std::fprintf(stderr, "[jcontigs]   %llu reads, %llu contigs\n",
                             (unsigned long long)s, (unsigned long long)emitted.load());
        }
    };
    std::vector<std::thread> th;
    th.emplace_back(reader);
    for (int i = 1; i < std::max(2, NT); i++) th.emplace_back(worker);
    for (auto& t : th) t.join();
    reads_close(in); std::fclose(out);
    if (std::rename((opath + ".tmp").c_str(), opath.c_str())) die("rename");

    size_t got = 0, seen_ch = 0;
    for (size_t i = 0; i < C.n; i++) { got += (taken[i].load() > 0); seen_ch += touched[i].load(); }
    std::fprintf(stderr, "[jcontigs] %llu reads (%.2f Gbp) -> %llu contigs (%.2f Gbp, %.2fx less sequence)\n",
                 (unsigned long long)seen.load(), seen_bp.load() / 1e9,
                 (unsigned long long)emitted.load(), emit_bp.load() / 1e9,
                 emit_bp.load() ? (double)seen_bp.load() / (double)emit_bp.load() : 0.0);
    std::fprintf(stderr, "[jcontigs] chunks: %zu total, %zu carried by some read, %zu with a contig "
                 "(%zu carried but NO read had %lld bp of anchor both sides -- these are the ones a "
                 "cross-read merge would be for), %llu windows rejected for short flanks\n",
                 C.n, seen_ch, got, seen_ch - got, (long long)MINF, (unsigned long long)short_flank.load());
    return 0;
}
