#include "cmd.h"
#include "graph.h"
#include "samples.h"
#include "snape.h"
#include "log.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

extern "C" {
#include "encoding.h"
}

static bool exists(const std::string& p) {
    FILE* f = std::fopen(p.c_str(), "rb");
    if (f) { std::fclose(f); return true; }
    return false;
}

static void usage() {
    std::fprintf(stderr,
        "usage:\n"
        "  snape run       --samples <sheet.tsv> <ref.fa> <dir> [--threads 16] [--level 5] [--margin 10]\n"
        "                  [--min-mult 3] [--supp-discover 5] [--supp-assign 2] [--max-bases B]\n"
        "  snape graph     --samples <sheet.tsv> <dir> [--level 5] [--margin 10] [--min-mult 3] [--threads 16]\n"
        "                  [--index N] [--pad-cores 50] [--max-bases B]\n"
        "                  N>0: <dir>/index.bin, N records per node; B>0: first B bases of each sample are read\n"
        "\n"
        "  the sheet is a three-column TSV, one row per read set:\n"
        "      name  role    path\n"
        "      N01   normal  hifi/normal.svb        exactly one normal\n"
        "      T01   tumor   hifi/t01.svb           one or more tumour lines\n"
        "  snape junctions <graph.bin> <sheet.tsv> <outdir> [threads=16] [pad_cores=50] [per_chunk=5] [min_flank=300] [min_vaf=0]\n"
        "  snape bubbles   <graph.bin> <sheet.tsv> [threads=16] [maxlen=6] [min_diff=25]\n"
        "  snape call      <ctg.bam|name=bam,name=bam,...> <normal.bam|-> <out.vcf> [min_mapq=10] [min_aln=0] [min_sv=50] "
        "[cluster=50] [min_supp=5] [shred_n=6] [shred_win=300] [max_mm=0.01] [shred_rep=0.40] "
        "[anchor_rep=0.50] [normal_dirty=0.25] [supp_assign=2]\n"
        "  snape pick      <ctg.fa[,ctg.fa...]> <calls.vcf> <normal reads> <out.fastq> [min_shared=5] [trim_bp=0]\n"
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
    if (n <= 0) die("run", "readlink /proc/self/exe");
    b[n] = 0;
    return b;
}

static std::string as_fastq(const std::string& reads, const std::string& self) {
    return is_svb(reads.c_str()) ? "'" + self + "' dump '" + reads + "' | " : std::string();
}

static std::string align(const std::string& reads, const std::string& ref, const std::string& out,
                         const std::string& tmp, int threads, const std::string& self) {
    const int st = threads > 4 ? 4 : 1;
    const std::string src = is_svb(reads.c_str()) ? "-" : "'" + reads + "'";
    return as_fastq(reads, self) + "minimap2 -ax map-hifi -t " + std::to_string(threads > 2 ? threads - 2 : 1) +
           " -Y --secondary=no '" + ref + "' " + src +
           " | samtools sort -@ " + std::to_string(st) + " -m 3G -T '" + tmp + "' -o '" + out + "' -";
}

static int align_checked(const std::string& reads, const std::string& ref, const std::string& bam,
                         const std::string& tmp, int threads, const std::string& self, const std::string& T) {
    if (sh(align(reads, ref, bam + ".tmp", tmp, threads, self))) return 1;
    if (sh("samtools quickcheck '" + bam + ".tmp'")) return 1;
    if (std::rename((bam + ".tmp").c_str(), bam.c_str())) die("run", "rename %s", bam.c_str());
    return sh("samtools index -@ " + T + " '" + bam + "'");
}

static int run_pipeline(const std::vector<Sample>& smp, const std::string& ref,
                        const std::string& dir, int level, int margin, int min_mult, int threads,
                        int supp_discover, int supp_assign, uint64_t max_bases) {
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

    if (exists(graph)) say("run", "skipped graph: %s exists", graph.c_str());
    else run_graph_reads(smp, dir, level, margin, min_mult, threads, max_bases);

    if (exists(drop)) say("run", "skipped bubbles: %s exists", drop.c_str());
    else {
        std::string a1 = "bubbles", a2 = graph, a3 = dir + "/samples.tsv", a4 = T;
        const char* av[] = {a1.c_str(), a2.c_str(), a3.c_str(), a4.c_str()};
        if (int rc = cmd_bubbles(4, av)) return rc;
    }

    if (exists(ctg)) say("run", "skipped junctions: %s exists", ctg.c_str());
    else {
        std::string a1 = "junctions", a2 = graph, a3 = dir + "/samples.tsv", a4 = dir, a5 = T;
        const char* av[] = {a1.c_str(), a2.c_str(), a3.c_str(), a4.c_str(), a5.c_str()};
        if (int rc = cmd_junctions(5, av)) return rc;
    }

    for (size_t i = 1; i < smp.size(); i++) {
        const std::string bam = dir + "/ctg." + smp[i].name + ".bam";
        if (exists(bam)) { say("run", "skipped align %s: %s exists", smp[i].name.c_str(), bam.c_str()); continue; }
        if (align_checked(dir + "/ctg." + smp[i].name + ".fa", ref, bam, dir + "/stc." + smp[i].name,
                          threads, self, T)) return 1;
    }

    if (exists(cand)) say("run", "skipped candidate call: %s exists", cand.c_str());
    else {
        std::string a1 = "call", a2 = ctg_bam, a3 = "-", a4 = cand, a5 = "10", a6 = "0", a7 = "50", a8 = "50",
                   a9 = std::to_string(supp_discover > 1 ? supp_discover - 1 : 1);
        const char* av[] = {a1.c_str(), a2.c_str(), a3.c_str(), a4.c_str(), a5.c_str(), a6.c_str(),
                            a7.c_str(), a8.c_str(), a9.c_str()};
        if (int rc = cmd_call(9, av)) return rc;
    }

    if (exists(nfq)) say("run", "skipped pick: %s exists", nfq.c_str());
    else {
        std::string a1 = "pick", a2 = ctg_fa, a3 = cand, a4 = normal, a5 = nfq;
        const char* av[] = {a1.c_str(), a2.c_str(), a3.c_str(), a4.c_str(), a5.c_str()};
        if (int rc = cmd_pick(5, av)) return rc;
    }

    if (exists(nbam)) say("run", "skipped normal align: %s exists", nbam.c_str());
    else {
        if (align_checked(nfq, ref, nbam, dir + "/stn", threads, self, T)) return 1;
    }

    if (exists(vcf)) say("run", "skipped call: %s exists", vcf.c_str());
    else {
        std::string a1 = "call", a2 = ctg_bam, a3 = nbam, a4 = vcf, a5 = "10", a6 = "0", a7 = "50", a8 = "50",
                   a9 = std::to_string(supp_discover), a10 = "6", a11 = "300", a12 = "0.01", a13 = "0.40",
                   a14 = "0.50", a15 = "0.25", a16 = std::to_string(supp_assign);
        const char* av[] = {a1.c_str(), a2.c_str(), a3.c_str(), a4.c_str(), a5.c_str(), a6.c_str(),
                            a7.c_str(), a8.c_str(), a9.c_str(), a10.c_str(), a11.c_str(), a12.c_str(),
                            a13.c_str(), a14.c_str(), a15.c_str(), a16.c_str()};
        if (int rc = cmd_call(16, av)) return rc;
    }
    say("run", "done -> %s", vcf.c_str());
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];

    if (cmd == "junctions") return cmd_junctions(argc - 1, (const char**)(argv + 1));
    if (cmd == "bubbles")   return cmd_bubbles(argc - 1, (const char**)(argv + 1));
    if (cmd == "call")      return cmd_call(argc - 1, (const char**)(argv + 1));
    if (cmd == "pick")      { LCP_INIT(); return cmd_pick(argc - 1, (const char**)(argv + 1)); }
    if (cmd == "pack")      return cmd_pack(argc - 1, (const char**)(argv + 1));
    if (cmd == "dump")      return cmd_dump(argc - 1, (const char**)(argv + 1));
    if (cmd == "svbidx")    return cmd_svbidx(argc - 1, (const char**)(argv + 1));

    int level = 5, margin = 10, min_mult = 3, threads = 16, supp_discover = 5, supp_assign = 2;
    int index_n = 8, pad_cores = 50;
    uint64_t max_bases = 0;
    std::string sheet;
    std::vector<std::string> pos;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        if (a.rfind("--", 0) != 0) { pos.push_back(a); continue; }
        if (i + 1 >= argc) die("run", "missing value for %s", a.c_str());
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
        else if (a == "--max-bases") max_bases = std::strtoull(v, nullptr, 10);
        else die("run", "arg: %s", a.c_str());
    }
    if (min_mult > (int)COV_MAX) die("run", "--min-mult %d exceeds the count width (%d)", min_mult, (int)COV_MAX);
    LCP_INIT();

    if ((cmd == "run" || cmd == "graph") && sheet.empty()) {
        die("run", "%s needs --samples <sheet.tsv>", cmd.c_str());
    }
    if (cmd == "run" && pos.size() == 2)
        return run_pipeline(samples_read(sheet), pos[0], pos[1], level, margin, min_mult, threads,
                            supp_discover, supp_assign, max_bases);
    if (cmd == "graph" && pos.size() == 1) {
        const std::vector<Sample> smp = samples_read(sheet);
        mkdir(pos[0].c_str(), 0755);
        samples_lock(smp, pos[0] + "/samples.tsv");
        run_graph_reads(smp, pos[0], level, margin, min_mult, threads, max_bases, index_n, pad_cores);
    }
    else { usage(); return 2; }
    return 0;
}
