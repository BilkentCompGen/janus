// snape pick <ctg.fa> <calls.vcf> <normal reads> <out.fastq> [threads=16] [min_shared=5]
//
// Which normal reads are worth aligning.
//
// `snape call` asks the normal reads at the reference locus whether a call is germline, and that one
// question is the largest single source of precision in the pipeline: hcc1954's contigs alone produce
// 1601 calls at 0.587/0.849, and asking the normal at each of their loci removes 554 of them -- 552
// false, 2 true -- for 0.895/0.847. But answering it by aligning the WHOLE normal makes the aligner the
// dominant cost of a caller whose entire premise is that the graph replaces alignment -- two hours
// against about forty minutes for every other stage put together.
//
// It also aligns far more than the question needs. The question is only ever asked at the loci the
// candidate calls sit on: about 2,400 of them, a few kilobases each. The reads that answer it are the
// normal reads carrying that same sequence, and the graph already has a way to recognise those
// without a reference -- they share LCP cores with the contig the call came from.
//
// So: take the contig each candidate names (`CTG=` in the VCF), take its GERMLINE ANCHORS -- the
// flanks, which are sequence the normal has too, as opposed to the somatic middle, which it does not
// -- parse them into cores, and emit every normal read sharing at least `min_shared` of them. What
// comes out is 1.56% of the normal on hcc1954, and it answers the question as well as the whole
// sample does: 0.887/0.858 against 0.891/0.858, eleven and a half minutes of alignment against two
// hours.
#include "pick.h"

#include "snape.h"
#include "unordered_dense.h"

extern "C" {
#include "encoding.h"
}

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
#include <vector>

namespace {

int LEVEL = 5;                                 // runtime, from SNAPE_LEVEL
constexpr int MARGIN = 10;                     // MUST match the graph build

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[pick] FATAL: %s\n", m); std::exit(1); }

// The contigs the candidates name. One per call, which is all that is needed: any contig that
// reaches the locus carries the same anchors.
ankerl::unordered_dense::set<std::string> wanted_contigs(const char* vcf) {
    ankerl::unordered_dense::set<std::string> w;
    FILE* f = std::fopen(vcf, "r"); if (!f) die("open calls.vcf");
    char* buf = nullptr; size_t cap = 0;
    while (getline(&buf, &cap, f) > 0) {
        if (buf[0] == '#') continue;
        const char* p = std::strstr(buf, ";CTG=");
        if (!p) continue;
        p += 5;
        const char* e = p;
        while (*e && *e != ';' && *e != '\t' && *e != '\n') e++;
        w.emplace(p, (size_t)(e - p));
    }
    std::free(buf); std::fclose(f);
    return w;
}

// A junction contig is [germline anchor | somatic run | germline anchor], and the header says how
// long each anchor is. Only the anchors go in: the somatic middle is by construction absent from the
// normal, so cores from it would select nothing and cost a lookup per read.
size_t anchor_cores(const char* fa, const ankerl::unordered_dense::set<std::string>& want,
                    ankerl::unordered_dense::set<uint64_t>& sigs) {
    FILE* f = std::fopen(fa, "r"); if (!f) die("open ctg.fa");
    char* buf = nullptr; size_t cap = 0;
    bool take = false; int64_t fl = 0, fr = 0;
    size_t nctg = 0;
    std::string seq, name;
    auto flush = [&]() {
        if (!take || seq.empty()) return;
        nctg++;
        const int64_t L = (int64_t)seq.size();
        const int64_t lo = std::min(fl, L), hi = std::max<int64_t>(0, L - fr);
        for (int pass = 0; pass < 2; pass++) {
            const int64_t a = pass ? hi : 0, b = pass ? L : lo;
            if (b - a < 200) continue;
            int n; const uint64_t *st, *en, *sg;
            extract(seq.data() + a, (int)(b - a), 0, LEVEL, &n, &st, &en, &sg);
            for (int i = 0; i < n; i++) {
                if (st[i] < (uint64_t)MARGIN || en[i] + MARGIN > (uint64_t)(b - a)) continue;
                sigs.insert(sg[i]);
            }
        }
        seq.clear();
    };
    while (getline(&buf, &cap, f) > 0) {
        if (buf[0] == '>') {
            flush();
            char* sp = std::strchr(buf + 1, ' ');
            const std::string nm(buf + 1, sp ? (size_t)(sp - buf - 1) : std::strlen(buf + 1) - 1);
            take = want.contains(nm);
            name = nm;
            fl = fr = 0;
            if (take && sp) {
                const char* p = std::strstr(sp, "flank=");
                if (p) std::sscanf(p + 6, "%lld,%lld", (long long*)&fl, (long long*)&fr);
            }
        } else if (take) {
            size_t n = std::strlen(buf);
            while (n && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) n--;
            seq.append(buf, n);
        }
    }
    flush();
    std::free(buf); std::fclose(f);
    return nctg;
}

}  // namespace

int cmd_pick(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: snape pick <ctg.fa> <calls.vcf> <normal reads> <out.fastq> "
                             "[threads=16] [min_shared=5]\n");
        return 2;
    }
    const std::string fa = argv[1], vcf = argv[2], reads = argv[3], out = argv[4];
    const int NT = argc > 5 ? std::atoi(argv[5]) : 16;
    // 5 shared cores, not 1: a read that really covers the locus shares dozens, while a single
    // shared core is what a repetitive anchor collects from all over the genome. Measured on
    // hcc1954: 1 picks 405,604 reads (3.58%), 5 picks 176,912 (1.56%), and the calls that come out
    // are identical to the record -- TP 949, FP 121 either way.
    const int MINSH = argc > 6 ? std::atoi(argv[6]) : 5;
    if (const char* lv = std::getenv("SNAPE_LEVEL")) LEVEL = std::atoi(lv);

    const auto want = wanted_contigs(vcf.c_str());
    if (want.empty()) { std::fprintf(stderr, "[pick] adayda CTG= yok, secilecek bir sey yok\n"); return 1; }
    ankerl::unordered_dense::set<uint64_t> sigs;
    const size_t nctg = anchor_cores(fa.c_str(), want, sigs);
    std::fprintf(stderr, "[pick] %zu aday contig -> %zu contig bulundu -> %zu cipa cekirdegi\n",
                 want.size(), nctg, sigs.size());
    if (sigs.empty()) die("no anchor cores");

    Reads in = reads_open(reads);
    FILE* fo = std::fopen((out + ".tmp").c_str(), "w"); if (!fo) die("open out");
    std::mutex out_mu;
    std::atomic<uint64_t> seen{0}, kept{0};

    // One reader fills batches, the workers take whole batches. The loop this replaces took the
    // input mutex once per read and copied the read AND its name while holding it, so every thread
    // queued behind every other thread's memcpy -- the same shape bubbles' collect had, where it
    // cost more than half the machine.
    constexpr size_t Q_CAP = 128, BATCH = 32;
    struct Rd { std::string name, seq; };
    std::mutex qmtx; std::condition_variable cv_full, cv_empty;
    std::queue<std::vector<Rd>> q;
    bool eof = false;

    auto reader = [&]() {
        std::vector<Rd> batch; batch.reserve(BATCH);
        Rd r;
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
        std::string buf;
        for (;;) {
            std::vector<Rd> batch;
            {
                std::unique_lock<std::mutex> lk(qmtx);
                cv_empty.wait(lk, [&]{ return !q.empty() || eof; });
                if (q.empty()) return;
                batch = std::move(q.front()); q.pop(); cv_full.notify_one();
            }
            for (const Rd& r : batch) {
                const uint64_t s = ++seen;
                int n; const uint64_t *st, *en, *sg;
                extract(r.seq.data(), (int)r.seq.size(), 0, LEVEL, &n, &st, &en, &sg);
                int hit = 0;
                for (int i = 0; i < n && hit < MINSH; i++) {
                    if (st[i] < (uint64_t)MARGIN || en[i] + MARGIN > r.seq.size()) continue;
                    hit += sigs.contains(sg[i]);
                }
                if (hit >= MINSH) {
                    buf.assign("@"); buf += r.name; buf += "\n"; buf += r.seq;
                    buf += "\n+\n"; buf.append(r.seq.size(), 'I'); buf += "\n";
                    std::lock_guard<std::mutex> lk(out_mu);
                    std::fwrite(buf.data(), 1, buf.size(), fo);
                    kept++;
                }
                if (s % 2000000 == 0)
                    std::fprintf(stderr, "[pick]   %llu okuma, %llu secildi\n",
                                 (unsigned long long)s, (unsigned long long)kept.load());
            }
        }
    };
    std::vector<std::thread> th;
    th.emplace_back(reader);
    for (int i = 1; i < std::max(2, NT); i++) th.emplace_back(worker);
    for (auto& t : th) t.join();
    reads_close(in); std::fclose(fo);
    if (std::rename((out + ".tmp").c_str(), out.c_str())) die("rename");

    const uint64_t sn = seen.load(), kp = kept.load();
    std::fprintf(stderr, "[pick] %llu normal okumanin %llu tanesi secildi (%%%.2f) -> %s\n",
                 (unsigned long long)sn, (unsigned long long)kp, sn ? 100.0 * kp / sn : 0.0, out.c_str());
    return 0;
}
