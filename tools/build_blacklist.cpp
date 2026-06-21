// Build a repeat-core blacklist from the reference assembly itself: LCP-parse the reference
// (forward + reverse-complement, owned-once per genomic position), hash each core with the SAME
// XXH3 core_sig the graph builder uses, count occurrences, and emit every sig seen > MINFREQ times
// (i.e. a core that occurs at >10 genomic loci -> a repeat). Output: u64 N + N sorted u64 sigs.
// The graph builder loads this via SVP_BLACKLIST and skips those cores (bypassing A->X->B to A->B).
//
// Built out-of-band (NOT part of `make`); links the same extract.o + lcptools as svp/trace so the
// LCP parse + signatures are byte-identical. Constants MUST match the graph build (LEVEL via
// SVP_BLACKLIST_LEVEL, default 5; MARGIN 10; FLANK 0).
#include "svp.h"                 // extract(), core_sig()
#include "unordered_dense.h"
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {
int LEVEL = 5;                                  // SVP_BLACKLIST_LEVEL; MUST match the graph --level
constexpr int MARGIN = 10;
constexpr uint64_t CHUNK = 4u << 20, PAD = 8192;

uint8_t COMP[256];
void init_comp() { for (int i = 0; i < 256; i++) COMP[i] = 'N';
    COMP[(uint8_t)'A']='T'; COMP[(uint8_t)'T']='A'; COMP[(uint8_t)'C']='G'; COMP[(uint8_t)'G']='C'; }

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[blacklist] FATAL: %s\n", m); std::exit(1); }
int env_int(const char* k, int d) { const char* v = std::getenv(k); return v ? std::atoi(v) : d; }

struct FaEntry { uint64_t len, offset; uint32_t linebases, linewidth; };
struct Faidx { int fd = -1; std::unordered_map<std::string, FaEntry> ix; std::vector<std::string> order; };

void load_faidx(const std::string& fa, Faidx& f) {
    f.fd = open(fa.c_str(), O_RDONLY); if (f.fd < 0) die("open fasta");
    FILE* fi = std::fopen((fa + ".fai").c_str(), "r"); if (!fi) die("open .fai (run: samtools faidx ref.fasta)");
    char name[512]; uint64_t len, off; uint32_t lb, lw;
    while (std::fscanf(fi, "%511s\t%lu\t%lu\t%u\t%u\n", name, &len, &off, &lb, &lw) == 5) {
        f.ix[name] = {len, off, lb, lw}; f.order.push_back(name); }
    std::fclose(fi);
}

std::string fa_fetch(const Faidx& f, const std::string& chrom, uint64_t start0, uint64_t want) {
    auto it = f.ix.find(chrom); if (it == f.ix.end()) return "";
    const FaEntry& e = it->second; if (start0 >= e.len) return "";
    if (start0 + want > e.len) want = e.len - start0;
    const uint64_t fstart = e.offset + (start0 / e.linebases) * e.linewidth + (start0 % e.linebases);
    const uint64_t nl = e.linewidth - e.linebases;
    const uint64_t fbytes = want + (want / e.linebases + 2) * nl;
    std::string raw(fbytes, '\0');
    const ssize_t got = pread(f.fd, &raw[0], fbytes, (off_t)fstart); if (got <= 0) return "";
    raw.resize((size_t)got);
    std::string out; out.reserve(want);
    for (char c : raw) { if (c=='\n'||c=='\r') continue; out.push_back(c>='a'&&c<='z'?(char)(c-32):c); if (out.size()==want) break; }
    return out;
}

std::string revcomp(const std::string& s) {
    std::string r(s.size(), 'N');
    for (size_t j = 0; j < s.size(); j++) r[j] = (char)COMP[(uint8_t)s[s.size()-1-j]];
    return r;
}

struct Job { std::string name; uint64_t clen, chunk_start; };

// count owned-once core sigs of one chunk (fwd + rc) into a thread-local map
void count_chunk(const Job& j, const Faidx& fa, ankerl::unordered_dense::map<uint64_t, uint32_t>& cnt) {
    const uint64_t own_lo = j.chunk_start, own_hi = std::min(j.clen, j.chunk_start + CHUNK);
    const uint64_t cs_lo = j.chunk_start > PAD ? j.chunk_start - PAD : 0;
    const uint64_t cs_hi = std::min(j.clen, own_hi + PAD);
    std::string s = fa_fetch(fa, j.name, cs_lo, cs_hi - cs_lo);
    if ((int)s.size() < 2*MARGIN+4) return;
    for (int pass = 0; pass < 2; pass++) {
        const std::string d = pass==0 ? s : revcomp(s);
        const int len = (int)d.size();
        int n; const uint64_t *st,*en,*sig; extract(d.data(), len, 0, LEVEL, &n, &st, &en, &sig);
        for (int i = 0; i < n; i++) {
            if (st[i] < (uint64_t)MARGIN || en[i] + MARGIN > (uint64_t)len) continue;
            const uint64_t gstart = pass==0 ? cs_lo+st[i] : cs_lo+(len-en[i]);
            if (gstart < own_lo || gstart >= own_hi) continue;              // each genomic position counted once
            cnt[sig[i]]++;
        }
    }
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: build_blacklist <ref.fa> <out.bin> [minfreq=10] [threads=16]\n"); return 1; }
    const std::string ref = argv[1], out = argv[2];
    const uint32_t minfreq = argc > 3 ? (uint32_t)std::atoi(argv[3]) : 10;   // blacklist sigs occurring > minfreq times
    const int NT = argc > 4 ? std::atoi(argv[4]) : 16;
    LEVEL = env_int("SVP_BLACKLIST_LEVEL", 5);
    init_comp();

    Faidx fa; load_faidx(ref, fa);
    std::vector<Job> jobs;
    for (auto& nm : fa.order) { const uint64_t L = fa.ix[nm].len;
        for (uint64_t cs = 0; cs < L; cs += CHUNK) jobs.push_back({nm, L, cs}); }
    std::fprintf(stderr, "[blacklist] ref=%s level=%d minfreq>%u : %zu chunks over %zu contigs on %d threads\n",
                 ref.c_str(), LEVEL, minfreq, jobs.size(), fa.order.size(), NT);

    std::vector<ankerl::unordered_dense::map<uint64_t, uint32_t>> local(NT);
    std::atomic<size_t> next{0};
    auto worker = [&](int wid) { size_t k;
        while ((k = next.fetch_add(1)) < jobs.size()) {
            count_chunk(jobs[k], fa, local[wid]);
            if ((k & 1023) == 0) std::fprintf(stderr, "\r[blacklist]   chunk %zu/%zu", k, jobs.size());
        } };
    std::vector<std::thread> th; for (int i = 0; i < NT; i++) th.emplace_back(worker, i); for (auto& t : th) t.join();
    std::fprintf(stderr, "\r[blacklist]   chunks done                         \n");

    // merge thread-local counts
    ankerl::unordered_dense::map<uint64_t, uint32_t> cnt; cnt.reserve(local[0].size() * 2);
    for (auto& m : local) { for (const auto& kv : m) cnt[kv.first] += kv.second; ankerl::unordered_dense::map<uint64_t, uint32_t>().swap(m); }

    // threshold + occurrence histogram
    std::vector<uint64_t> bl; uint64_t total_occ = 0; uint64_t hist[6] = {0};   // 1,2,3-5,6-10,11-100,>100
    for (const auto& kv : cnt) { total_occ += kv.second;
        const uint32_t c = kv.second;
        if (c==1) hist[0]++; else if (c==2) hist[1]++; else if (c<=5) hist[2]++; else if (c<=10) hist[3]++; else if (c<=100) hist[4]++; else hist[5]++;
        if (c > minfreq) bl.push_back(kv.first); }
    std::sort(bl.begin(), bl.end());

    FILE* f = std::fopen((out + ".tmp").c_str(), "wb"); if (!f) die("open out");
    uint64_t N = bl.size(); std::fwrite(&N, 8, 1, f); if (N) std::fwrite(bl.data(), 8, N, f); std::fclose(f);
    if (std::rename((out + ".tmp").c_str(), out.c_str())) die("rename out");

    std::fprintf(stderr, "[blacklist] distinct cores=%zu  total occurrences=%lu (incl. fwd+rc, owned-once)\n", cnt.size(), (unsigned long)total_occ);
    std::fprintf(stderr, "[blacklist]   occ histogram: =1:%lu  =2:%lu  3-5:%lu  6-10:%lu  11-100:%lu  >100:%lu\n",
                 (unsigned long)hist[0],(unsigned long)hist[1],(unsigned long)hist[2],(unsigned long)hist[3],(unsigned long)hist[4],(unsigned long)hist[5]);
    std::fprintf(stderr, "[blacklist] blacklisted (>%u): %lu sigs (%.2f%% of distinct cores) -> %s\n",
                 minfreq, (unsigned long)N, 100.0*(double)N/(double)cnt.size(), out.c_str());
    return 0;
}
