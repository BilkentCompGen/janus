#include "cmd.h"
#include "postings.h"
#include "log.h"
#include <algorithm>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

static ankerl::unordered_dense::set<std::string> wanted_contigs(const char* vcf) {
    ankerl::unordered_dense::set<std::string> w;
    FILE* f = xopen(vcf, "r");
    char* buf = nullptr;
    size_t cap = 0;
    while (getline(&buf, &cap, f) > 0) {
        if (buf[0] == '#') continue;
        const char* p = std::strstr(buf, ";CTG=");
        if (!p) continue;
        p += 5;
        const char* e = p;
        while (*e && *e != ';' && *e != '\t' && *e != '\n') e++;
        w.emplace(p, (size_t)(e - p));
    }
    std::free(buf);
    std::fclose(f);
    return w;
}

static void flush_ctg(bool take, std::string& seq, int64_t fl, int64_t fr,
                       size_t& nctg, ankerl::unordered_dense::set<uint64_t>& sigs) {
    if (!take || seq.empty()) return;
    nctg++;
    const int64_t L = (int64_t)seq.size();
    const int64_t lo = std::min(fl, L), hi = std::max<int64_t>(0, L - fr);
    constexpr int MARGIN = 10;
    for (int pass = 0; pass < 2; pass++) {
        const int64_t a = pass ? hi : 0, b = pass ? L : lo;
        if (b - a < 200) continue;
        int n;
        const uint64_t *st, *en, *sg;
        extract(seq.data() + a, (int)(b - a), level_env(), &n, &st, &en, &sg);
        for (int i = 0; i < n; i++) {
            if (st[i] < (uint64_t)MARGIN || en[i] + MARGIN > (uint64_t)(b - a)) continue;
            sigs.insert(sg[i]);
        }
    }
    seq.clear();
}

static size_t anchor_cores(const char* fa, const ankerl::unordered_dense::set<std::string>& want,
                            ankerl::unordered_dense::set<uint64_t>& sigs) {
    FILE* f = xopen(fa, "r");
    char* buf = nullptr;
    size_t cap = 0;
    bool take = false;
    int64_t fl = 0, fr = 0;
    size_t nctg = 0;
    std::string seq, name;
    while (getline(&buf, &cap, f) > 0) {
        if (buf[0] == '>') {
            flush_ctg(take, seq, fl, fr, nctg, sigs);
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
    flush_ctg(take, seq, fl, fr, nctg, sigs);
    std::free(buf);
    std::fclose(f);
    return nctg;
}

int cmd_pick(int argc, const char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: snape pick <ctg.fa[,ctg.fa...]> <calls.vcf> <normal reads> <out.fastq> "
                             "[min_shared=5] [trim_bp=0]\n");
        return 2;
    }
    const std::string fa = argv[1], vcf = argv[2], reads = argv[3], out = argv[4];
    const int MINSH = argc > 5 ? std::atoi(argv[5]) : 5;

    const auto want = wanted_contigs(vcf.c_str());
    if (want.empty()) { say("pick", "no CTG= in candidates, nothing to pick"); return 1; }

    ankerl::unordered_dense::set<uint64_t> sigs;
    size_t nctg = 0;
    for (size_t i = 0; i <= fa.size();) {
        const size_t c = std::min(fa.find(',', i), fa.size());
        nctg += anchor_cores(fa.substr(i, c - i).c_str(), want, sigs);
        i = c + 1;
    }
    say("pick", "%zu candidate contigs -> %zu found -> %zu anchor cores", want.size(), nctg, sigs.size());
    if (sigs.empty()) die("pick", "no anchor cores");

    const std::string dir = dir_of(fa);
    uint64_t K = 0;
    const std::vector<uint64_t> gsig = graph_sigs(dir + "/graph.bin", K);
    SigLookup sl;
    sl.build(gsig);
    Index ix;
    ix.load(dir + "/index.bin");
    if (ix.K != K) die("pick", "index.bin and graph.bin disagree on the node count");
    say("pick", "index.bin: %llu records, N=%u", (unsigned long long)ix.off[ix.K], ix.N);

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

    SvbFile rd;
    rd.open(reads);
    FILE* fo = xopen((out + ".tmp").c_str(), "w");
    std::string seq, buf;
    uint64_t kept = 0, kept_bp = 0;
    for (uint32_t r : take) {
        if (!rd.read(r, seq)) continue;
        const size_t Ln = seq.size();
        buf.assign("@r");
        buf += std::to_string(r);
        buf += "\n";
        buf.append(seq, 0, Ln);
        buf += "\n+\n";
        buf.append(Ln, 'I');
        buf += "\n";
        wr(fo, buf.data(), buf.size());
        kept++;
        kept_bp += Ln;
    }
    rd.close();
    atomic_close(fo, (out + ".tmp").c_str(), out.c_str());

    say("pick", "%zu of %zu anchor cores in the graft -> %zu candidate reads -> "
               "%llu selected (>=%d shared cores), %.3f Gbp -> %s",
        found, sigs.size(), hits.size(), (unsigned long long)kept, MINSH,
        kept_bp / 1e9, out.c_str());
    return 0;
}
