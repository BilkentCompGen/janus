// svp: 2-colour LCP-core graph built FROM READS (normal vs tumour reads), then the somatic material
// carved out of those reads as junction contigs. READS ONLY -- no assembly, and no reference until
// the contigs are handed to an aligner, outside this binary.
//
//   svp graph     <normal> <tumour> <dir>      reads           -> dir/graph.bin
//   svp bubbles   <graph.bin> <normal> <tumour> graph + reads  -> dir/drop.sigs (sub-SV cores)
//   svp junctions <graph.bin> <reads> <out.fa> graph + reads   -> junction contigs
//   svp call      <ctg.bam> <out.vcf>          aligned contigs -> breakpoints
//   svp pack      <out.svb>                    FASTQ (stdin)   -> 2-bit .svb
//   svp dump      <in.svb>                     .svb            -> FASTQ (stdout)
//   svp refine    <dir>                        graph.bin     -> refined.bin (topology EC)
//   svp run       <normal> <tumour> <dir>      graph -> refine, each skipped if done
//
// `graph` and `junctions` are the two halves of calling; `pack`/`dump` are the read format both of
// them consume, kept here rather than as separate binaries so there is exactly one thing to build.
#include "bubbles.h"
#include "call.h"
#include "graph.h"
#include "junctions.h"
#include "refinement.h"
#include "svb.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

extern "C" {
#include "encoding.h"
}

static bool exists(const std::string& p) { FILE* f = std::fopen(p.c_str(), "rb"); if (f) { std::fclose(f); return true; } return false; }

static void usage() {
    std::fprintf(stderr,
        "usage:\n"
        "  svp graph     <normal> <tumour> <dir> [--level 5] [--margin 10] [--min-mult 3] [--threads 16]\n"
        "  svp junctions <graph.bin> <reads.fastq> <out.fa> [threads] [pad_cores] [per_chunk] [min_flank] [min_vaf]\n"
        "  svp bubbles   <graph.bin> <normal> <tumour> [threads=16] [k=21] [w=22] [maxlen=40] [max_dlen=12] [max_unmatched=20]\n"
        "  svp call      <ctg.bam> <normal.bam|-> <out.vcf> [min_mapq=10] [min_aln=1000] [min_sv=50] [cluster=50] [min_supp=3]\n"
        "  svp pack      <out.svb> < reads.fastq\n"
        "  svp dump      <in.svb> > reads.fastq\n"
        "  svp refine    <dir>                                  (graph.bin -> refined.bin)\n"
        "  svp run       <normal> <tumour> <dir> [graph opts]   (graph -> refine; skips completed)\n");
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];

    // These take their arguments positionally, with no --flags to strip, so they are dispatched
    // before the graph option parser runs. junctions needs the LCP alphabet; pack/dump do not.
    if (cmd == "junctions") { LCP_INIT(); return cmd_junctions(argc - 1, argv + 1); }
    if (cmd == "bubbles")   { LCP_INIT(); return cmd_bubbles(argc - 1, argv + 1); }
    if (cmd == "call")      return cmd_call(argc - 1, argv + 1);
    if (cmd == "pack")      return cmd_pack(argc - 1, argv + 1);
    if (cmd == "dump")      return cmd_dump(argc - 1, argv + 1);

    int level = 5, margin = 10, min_mult = 3, threads = 16;
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
        else { std::fprintf(stderr, "[FATAL] arg: %s\n", a.c_str()); std::abort(); }
    }
    LCP_INIT();

    if (cmd == "graph" && pos.size() == 3) run_graph_reads(pos[0], pos[1], pos[2], level, margin, min_mult, threads);
    else if (cmd == "refine" && pos.size() == 1) run_refine(pos[0]);
    else if (cmd == "run" && pos.size() == 3) {
        const std::string& dir = pos[2];
        if (exists(dir + "/graph.bin")) std::fprintf(stderr, "[skip] graph: %s/graph.bin exists\n", dir.c_str());
        else run_graph_reads(pos[0], pos[1], dir, level, margin, min_mult, threads);
        if (exists(dir + "/refined.bin")) std::fprintf(stderr, "[skip] refine: %s/refined.bin exists\n", dir.c_str());
        else run_refine(dir);
    } else { usage(); return 2; }
    return 0;
}
