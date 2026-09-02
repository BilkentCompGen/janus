#include "pick.h"
#include "postings.h"
#include <algorithm>

#include "snape.h"
#include "unordered_dense.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

[[noreturn]] void die(const char* m) { std::fprintf(stderr, "[pick] FATAL: %s\n", m); std::exit(1); }

ankerl::unordered_dense::set<uint64_t> load_drop(const std::string& path) {
    ankerl::unordered_dense::set<uint64_t> d;
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return d;
    char* buf = nullptr; size_t cap = 0;
    while (getline(&buf, &cap, f) > 0) d.insert(std::strtoull(buf, nullptr, 10));
    std::free(buf); std::fclose(f);
    return d;
}

}

int cmd_pick(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "usage: snape pick <ctg.fa[,ctg.fa...]> <calls.vcf> <normal reads> <out.fastq> "
                             "[threads=16] [min_shared=5]\n");
        return 2;
    }
    const std::string fa = argv[1], reads = argv[3], out = argv[4];
    const int MINSH = argc > 6 ? std::atoi(argv[6]) : 5;

    const size_t cut = fa.find_last_of('/');
    const std::string dir = cut == std::string::npos ? std::string(".") : fa.substr(0, cut);

    uint32_t S = 0; uint64_t K = 0;
    std::vector<uint64_t> gsig;
    const std::vector<uint8_t> nrm = graph_normal_only(dir + "/graph.bin", S, K, gsig);
    const auto drop = load_drop(dir + "/drop.n.sigs");

    Index ix; ix.load(dir + "/index.bin");
    if (ix.K != K) die("index.bin and graph.bin disagree on the node count");

    SvbFile rd; rd.open(reads);
    std::vector<uint32_t> hits(rd.nreads, 0);

    uint64_t nnrm = 0, ndrop = 0;
    for (uint64_t i = 0; i < K; i++) {
        if (!nrm[i]) continue;
        if (drop.contains(gsig[i])) { ndrop++; continue; }
        nnrm++;
        const uint32_t n = ix.count((uint32_t)i);
        for (uint32_t j = 0; j < n; j++) {
            const Post p = ix.get((uint32_t)i, j);
            if (p.samp != 0 || p.read >= rd.nreads) continue;
            hits[p.read]++;
        }
    }
    std::fprintf(stderr, "[pick] %llu normal-only cekirdek (%llu drop.n.sigs ile elendi), index.bin %llu kayit N=%u\n",
                 (unsigned long long)nnrm, (unsigned long long)ndrop,
                 (unsigned long long)ix.off[ix.K], ix.N);

    FILE* fo = std::fopen((out + ".tmp").c_str(), "w"); if (!fo) die("open out");
    std::string seq, buf;
    uint64_t cand = 0, kept = 0, kept_bp = 0;
    for (uint32_t r = 0; r < rd.nreads; r++) {
        if (!hits[r]) continue;
        cand++;
        if (hits[r] < (uint32_t)MINSH) continue;
        if (!rd.read(r, seq)) continue;
        buf.assign("@r"); buf += std::to_string(r); buf += "\n";
        buf += seq;
        buf += "\n+\n"; buf.append(seq.size(), 'I'); buf += "\n";
        std::fwrite(buf.data(), 1, buf.size(), fo);
        kept++; kept_bp += seq.size();
    }
    rd.close(); std::fclose(fo);
    if (std::rename((out + ".tmp").c_str(), out.c_str())) die("rename");

    std::fprintf(stderr, "[pick] %llu aday okuma -> %llu secildi (>=%d ortak cekirdek), %.3f Gbp -> %s\n",
                 (unsigned long long)cand, (unsigned long long)kept, MINSH,
                 kept_bp / 1e9, out.c_str());
    return 0;
}
