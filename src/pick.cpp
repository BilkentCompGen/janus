#include "pick.h"
#include "postings.h"
#include <algorithm>

#include "snape.h"
#include "unordered_dense.h"

extern "C" {
#include "encoding.h"
}

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

namespace {

int LEVEL = 5;
constexpr int MARGIN = 10;

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[pick] FATAL: %s\n", m); std::exit(1); }

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

}

int cmd_pick(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: snape pick <ctg.fa[,ctg.fa...]> <calls.vcf> <normal reads> <out.fastq> "
                             "[min_shared=5] [trim_bp=0]\n");
        return 2;
    }
    const std::string fa = argv[1], vcf = argv[2], reads = argv[3], out = argv[4];
    const int MINSH = argc > 5 ? std::atoi(argv[5]) : 5;
    const int64_t TRIM = argc > 6 ? std::atoll(argv[6]) : 0;
    if (const char* lv = std::getenv("SNAPE_LEVEL")) LEVEL = std::atoi(lv);

    const auto want = wanted_contigs(vcf.c_str());
    if (want.empty()) { std::fprintf(stderr, "[pick] no CTG= in candidates, nothing to pick\n"); return 1; }

    ankerl::unordered_dense::set<uint64_t> sigs;
    size_t nctg = 0;
    for (size_t i = 0; i <= fa.size();) {
        const size_t c = std::min(fa.find(',', i), fa.size());
        nctg += anchor_cores(fa.substr(i, c - i).c_str(), want, sigs);
        i = c + 1;
    }
    std::fprintf(stderr, "[pick] %zu candidate contigs -> %zu found -> %zu anchor cores\n",
                 want.size(), nctg, sigs.size());
    if (sigs.empty()) die("no anchor cores");

    const size_t cut = fa.find_last_of('/');
    const std::string dir = cut == std::string::npos ? std::string(".") : fa.substr(0, cut);
    uint64_t K = 0;
    const std::vector<uint64_t> gsig = graph_sigs(dir + "/graph.bin", K);
    SigLookup sl; sl.build(gsig);
    Index ix; ix.load(dir + "/index.bin");
    if (ix.K != K) die("index.bin and graph.bin disagree on the node count");
    std::fprintf(stderr, "[pick] index.bin: %llu records, N=%u\n",
                 (unsigned long long)ix.off[ix.K], ix.N);

    ankerl::unordered_dense::map<uint32_t, uint32_t> hits;
    size_t found = 0;
    for (uint64_t sg : sigs) {
        const uint32_t node = sl.find(sg);
        if (node == UINT32_MAX) continue;
        found++;
        const uint32_t n = ix.count(node);
        for (uint32_t j = 0; j < n; j++) {
            const Post p = ix.get(node, j);
            if (p.samp != 0) continue;
            hits[p.read]++;
        }
    }
    std::vector<uint32_t> take;
    for (const auto& kv : hits) if (kv.second >= (uint32_t)MINSH) take.push_back(kv.first);
    std::sort(take.begin(), take.end());

    SvbFile rd; rd.open(reads);
    FILE* fo = std::fopen((out + ".tmp").c_str(), "w"); if (!fo) die("open out");
    std::string seq, buf;
    uint64_t kept = 0, kept_bp = 0;
    for (uint32_t r : take) {
        if (!rd.read(r, seq)) continue;
        uint64_t a = 0, b = seq.size();
        if (TRIM > 0) {
            a = 0; b = seq.size();
        }
        const size_t Ln = (size_t)(b - a);
        buf.assign("@r"); buf += std::to_string(r); buf += "\n";
        buf.append(seq, (size_t)a, Ln);
        buf += "\n+\n"; buf.append(Ln, 'I'); buf += "\n";
        std::fwrite(buf.data(), 1, buf.size(), fo);
        kept++; kept_bp += Ln;
    }
    rd.close(); std::fclose(fo);
    if (std::rename((out + ".tmp").c_str(), out.c_str())) die("rename");

    std::fprintf(stderr, "[pick] %zu of %zu anchor cores in the graft -> %zu candidate reads -> "
                 "%llu selected (>=%d shared cores), %.3f Gbp -> %s\n",
                 found, sigs.size(), hits.size(), (unsigned long long)kept, MINSH,
                 kept_bp / 1e9, out.c_str());
    return 0;
}
