// svp (lcp-cortex): somatic SV caller over a 2-colour LCP-core graph built FROM READS
// (normal vs tumor HiFi reads). READS ONLY; no assembly in the pipeline. Runs as two
// resumable, checkpointed stages so the expensive LCP pass is paid once:
//   svp graph    <normal.fq.gz> <tumor.fq.gz> <dir>  -> dir/graph.bin     (2-pass bloom build:
//                                                        colour + node-id edges + inline payload)
//   svp simplify <dir>                               -> dir/simplify.bin  (in-RAM unitig
//                                                        compaction, then bubble-zipped graph)
//   svp run      <normal.fq.gz> <tumor.fq.gz> <dir>  -> both, skipping completed stages
#include "graph.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

extern "C" {
#include "encoding.h"
}

static bool exists(const std::string& p) { FILE* f = std::fopen(p.c_str(), "rb"); if (f) { std::fclose(f); return true; } return false; }

static void usage() {
    std::fprintf(stderr,
        "usage:\n"
        "  svp graph    <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5] [--flank 10] [--margin 10] [--min-mult 10] [--threads 16]\n"
        "  svp simplify <dir>\n"
        "  svp run      <normal.fq.gz> <tumor.fq.gz> <dir> [all opts]   (skips completed stages)\n");
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];
    int level = 5, flank = 10, margin = 10, min_mult = 10, threads = 16;
    std::vector<std::string> pos;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        if (a.rfind("--", 0) != 0) { pos.push_back(a); continue; }
        if (i + 1 >= argc) { std::fprintf(stderr, "[FATAL] missing value for %s\n", a.c_str()); std::abort(); }
        const char* v = argv[++i];
        if (a == "--level") level = std::atoi(v);
        else if (a == "--flank") flank = std::atoi(v);
        else if (a == "--margin") margin = std::atoi(v);
        else if (a == "--min-mult") min_mult = std::atoi(v);
        else if (a == "--threads") threads = std::atoi(v);
        else { std::fprintf(stderr, "[FATAL] arg: %s\n", a.c_str()); std::abort(); }
    }
    LCP_INIT();

    if (cmd == "graph" && pos.size() == 3) run_graph_reads(pos[0], pos[1], pos[2], level, flank, margin, min_mult, threads);
    else if (cmd == "simplify" && pos.size() == 1) run_simplify(pos[0]);
    else if (cmd == "run" && pos.size() == 3) {
        const std::string& dir = pos[2];
        if (exists(dir + "/graph.bin")) std::fprintf(stderr, "[skip] graph: %s/graph.bin exists\n", dir.c_str());
        else run_graph_reads(pos[0], pos[1], dir, level, flank, margin, min_mult, threads);
        if (exists(dir + "/simplify.bin")) std::fprintf(stderr, "[skip] simplify: %s/simplify.bin exists\n", dir.c_str());
        else run_simplify(dir);
    } else { usage(); return 2; }
    return 0;
}
