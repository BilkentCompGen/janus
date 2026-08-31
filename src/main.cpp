#include "bubbles.h"
#include "call.h"
#include "graph.h"
#include "junctions.h"
#include "pick.h"
#include "samples.h"
#include "snape.h"
#include "svb.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

extern "C" {
#include "encoding.h"
}

static bool exists(const std::string& p) { FILE* f = std::fopen(p.c_str(), "rb"); if (f) { std::fclose(f); return true; } return false; }

static void usage() {
    std::fprintf(stderr,
        "usage:\n"
        "  snape run       --samples <sheet.tsv> <ref.fa> <dir> [--threads 16] [--level 5] [--margin 10]\n"
        "                  [--min-mult 3] [--supp-discover 5] [--supp-assign 2]\n"
        "  snape graph     --samples <sheet.tsv> <dir> [--level 5] [--margin 10] [--min-mult 3] [--threads 16]\n"
        "                  [--index N] [--pad-cores 50]   N>0: <dir>/index.bin, dugum basina N kayit\n"
        "\n"
        "  the sheet is a three-column TSV, one row per read set:\n"
        "      name  role    path\n"
        "      N01   normal  hifi/normal.svb        exactly one normal\n"
        "      T01   tumor   hifi/t01.svb           one or more tumour lines\n"
        "  snape junctions <graph.bin> <sheet.tsv> <outdir> [threads=16] [pad_cores=50] [per_chunk=5] [min_flank=300] [min_vaf=0]\n"
        "  snape bubbles   <graph.bin> <sheet.tsv> [threads=16] [k=21] [w=22] [maxlen=40] [max_dbp=50] [max_unmatched=20]\n"
        "  snape call      <ctg.bam> <normal.bam|-> <out.vcf> [min_mapq=10] [min_aln=0] [min_sv=50] [cluster=50] [min_supp=5]\n"
        "  snape pick      <ctg.fa> <calls.vcf> <normal reads> <out.fastq> [threads=16] [min_shared=5]\n"
        "  snape pack      <out.svb> < reads.fastq\n"
        "  snape dump      <in.svb> > reads.fastq\n"
        "  snape svbidx    <in.svb>            per-read offsets -> in.svb.idx\n");
}

static int sh(const std::string& cmd) {
    std::fprintf(stderr, "[run] $ %s\n", cmd.c_str());
    const int rc = std::system(cmd.c_str());
    if (rc != 0) std::fprintf(stderr, "[run] FATAL: rc=%d\n", rc);
    return rc;
}

static std::string self_path() {
    char b[4096];
    const ssize_t n = readlink("/proc/self/exe", b, sizeof b - 1);
    if (n <= 0) { std::fprintf(stderr, "[run] FATAL: readlink /proc/self/exe\n"); std::exit(1); }
    b[n] = 0;
    return b;
}

static bool is_svb(const std::string& p) { return p.size() >= 4 && p.compare(p.size() - 4, 4, ".svb") == 0; }

static std::string as_fastq(const std::string& reads, const std::string& self) {
    return is_svb(reads) ? "'" + self + "' dump '" + reads + "' | " : std::string();
}

static std::string align(const std::string& reads, const std::string& ref, const std::string& out,
                         const std::string& tmp, int threads, const std::string& self) {
    const int st = threads > 4 ? 4 : 1;
    const std::string src = is_svb(reads) ? "-" : "'" + reads + "'";
    return as_fastq(reads, self) + "minimap2 -ax map-hifi -t " + std::to_string(threads > 2 ? threads - 2 : 1) +
           " -Y --secondary=no '" + ref + "' " + src +
           " | samtools sort -@ " + std::to_string(st) + " -m 3G -T '" + tmp + "' -o '" + out + "' -";
}

static int align_checked(const std::string& reads, const std::string& ref, const std::string& bam,
                         const std::string& tmp, int threads, const std::string& self, const std::string& T) {
    if (sh(align(reads, ref, bam + ".tmp", tmp, threads, self))) return 1;
    if (sh("samtools quickcheck '" + bam + ".tmp'")) return 1;
    if (std::rename((bam + ".tmp").c_str(), bam.c_str())) { std::fprintf(stderr, "[run] FATAL: rename %s\n", bam.c_str()); return 1; }
    return sh("samtools index -@ " + T + " '" + bam + "'");
}

static int run_pipeline(const std::vector<Sample>& smp, const std::string& ref,
                        const std::string& dir, int level, int margin, int min_mult, int threads,
                        int supp_discover, int supp_assign) {
    mkdir(dir.c_str(), 0755);
    samples_lock(smp, dir + "/samples.tsv");
    const std::string self = self_path();
    const std::string normal = smp[0].path;
    const std::string graph = dir + "/graph.bin", drop = dir + "/drop.sigs";
    const std::string ctg = dir + "/ctg." + smp.back().name + ".fa";

    std::string ctg_fa, ctg_bam;
    for (size_t i = 1; i < smp.size(); i++) {
        if (i > 1) { ctg_fa += ','; ctg_bam += ','; }
        ctg_fa  += dir + "/ctg." + smp[i].name + ".fa";
        ctg_bam += smp[i].name + "=" + dir + "/ctg." + smp[i].name + ".bam";
    }
    const std::string nbam = dir + "/normal.bam", vcf = dir + "/calls.vcf";
    const std::string cand = dir + "/candidates.vcf", nfq = dir + "/normal_at_candidates.fastq";
    const std::string T = std::to_string(threads);

    if (exists(graph)) std::fprintf(stderr, "[run] atlandi graph: %s var\n", graph.c_str());
    else run_graph_reads(smp, dir, level, margin, min_mult, threads);

    if (exists(drop)) std::fprintf(stderr, "[run] atlandi bubbles: %s var\n", drop.c_str());
    else {
        std::vector<std::string> a{"bubbles", graph, dir + "/samples.tsv", T};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_bubbles((int)v.size(), v.data())) return rc;
    }

    if (exists(ctg)) std::fprintf(stderr, "[run] atlandi junctions: %s var\n", ctg.c_str());
    else {
        std::vector<std::string> a{"junctions", graph, dir + "/samples.tsv", dir, T};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_junctions((int)v.size(), v.data())) return rc;
    }

    for (size_t i = 1; i < smp.size(); i++) {
        const std::string bam = dir + "/ctg." + smp[i].name + ".bam";
        if (exists(bam)) { std::fprintf(stderr, "[run] atlandi align %s: %s var\n", smp[i].name.c_str(), bam.c_str()); continue; }
        if (align_checked(dir + "/ctg." + smp[i].name + ".fa", ref, bam, dir + "/stc." + smp[i].name,
                          threads, self, T)) return 1;
    }

    if (exists(cand)) std::fprintf(stderr, "[run] atlandi aday call: %s var\n", cand.c_str());
    else {
        std::vector<std::string> a{"call", ctg_bam, "-", cand, "10", "0", "50", "50",
                                   std::to_string(supp_discover > 1 ? supp_discover - 1 : 1)};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_call((int)v.size(), v.data())) return rc;
    }

    if (exists(nfq)) std::fprintf(stderr, "[run] atlandi pick: %s var\n", nfq.c_str());
    else {
        std::vector<std::string> a{"pick", ctg_fa, cand, normal, nfq, T};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_pick((int)v.size(), v.data())) return rc;
    }

    if (exists(nbam)) std::fprintf(stderr, "[run] atlandi normal align: %s var\n", nbam.c_str());
    else {
        if (align_checked(nfq, ref, nbam, dir + "/stn", threads, self, T)) return 1;
    }

    if (exists(vcf)) std::fprintf(stderr, "[run] atlandi call: %s var\n", vcf.c_str());
    else {
        std::vector<std::string> a{"call", ctg_bam, nbam, vcf, "10", "0", "50", "50",
                                   std::to_string(supp_discover), "6", "300", "0.01", "0.40",
                                   "0.50", "0.25", std::to_string(supp_assign)};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_call((int)v.size(), v.data())) return rc;
    }
    std::fprintf(stderr, "[run] bitti -> %s\n", vcf.c_str());
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];

    if (cmd == "junctions") { LCP_INIT(); return cmd_junctions(argc - 1, argv + 1); }
    if (cmd == "bubbles")   { LCP_INIT(); return cmd_bubbles(argc - 1, argv + 1); }
    if (cmd == "call")      return cmd_call(argc - 1, argv + 1);
    if (cmd == "pick")      { LCP_INIT(); return cmd_pick(argc - 1, argv + 1); }
    if (cmd == "pack")      return cmd_pack(argc - 1, argv + 1);
    if (cmd == "dump")      return cmd_dump(argc - 1, argv + 1);
    if (cmd == "svbidx")    return cmd_svbidx(argc - 1, argv + 1);

    int level = 5, margin = 10, min_mult = 3, threads = 16, supp_discover = 5, supp_assign = 2;
    int index_n = 8, pad_cores = 50;
    std::string sheet;
    std::vector<std::string> pos;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        if (a.rfind("--", 0) != 0) { pos.push_back(a); continue; }
        if (i + 1 >= argc) { std::fprintf(stderr, "[FATAL] missing value for %s\n", a.c_str()); std::abort(); }
        const char* v = argv[++i];
        if (a == "--level") level = std::atoi(v);
        else if (a == "--margin") margin = std::atoi(v);
        else if (a == "--min-mult") min_mult = std::atoi(v);
        else if (a == "--threads") threads = std::atoi(v);
        else if (a == "--samples") sheet = v;
        else if (a == "--supp-discover") supp_discover = std::atoi(v);
        else if (a == "--supp-assign") supp_assign = std::atoi(v);
        else if (a == "--index") index_n = std::atoi(v);
        else if (a == "--pad-cores") pad_cores = std::atoi(v);
        else { std::fprintf(stderr, "[FATAL] arg: %s\n", a.c_str()); std::abort(); }
    }
    if (min_mult > (int)COV_MAX) { std::fprintf(stderr, "[FATAL] --min-mult %d exceeds the count width (%d)\n", min_mult, (int)COV_MAX); std::abort(); }
    LCP_INIT();

    if ((cmd == "run" || cmd == "graph") && sheet.empty()) {
        std::fprintf(stderr, "[FATAL] %s needs --samples <sheet.tsv>\n", cmd.c_str());
        usage(); return 2;
    }
    if (cmd == "run" && pos.size() == 2)
        return run_pipeline(samples_read(sheet), pos[0], pos[1], level, margin, min_mult, threads,
                            supp_discover, supp_assign);
    if (cmd == "graph" && pos.size() == 1) {
        const std::vector<Sample> smp = samples_read(sheet);
        mkdir(pos[0].c_str(), 0755);
        samples_lock(smp, pos[0] + "/samples.tsv");
        run_graph_reads(smp, pos[0], level, margin, min_mult, threads, index_n, pad_cores);
    }
    else { usage(); return 2; }
    return 0;
}
