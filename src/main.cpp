// svp (lcp-cortex): somatic SV caller over a 2-colour LCP-core graph built FROM READS
// (normal vs tumor HiFi reads). READS ONLY; no assembly, no reference. Three resumable,
// checkpointed stages so the expensive LCP pass is paid once:
//   svp graph  <normal.fq.gz> <tumor.fq.gz> <dir>  -> dir/graph.bin    (count + fill: colour + edges + payload)
//   svp refine <dir>                               -> dir/refined.bin  (topology + counts error-correction)
//   svp call   <dir>                               -> dir/calls.*      (somatic variant calls)
//   svp run    <normal.fq.gz> <tumor.fq.gz> <dir>  -> all three, skipping completed stages
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
        "  svp graph  <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5] [--margin 10] [--min-mult 3] [--threads 16]\n"
        "  svp refine <dir>                               (graph.bin -> refined.bin: topology+counts EC)\n"
        "  svp call   <dir>                               (refined.bin -> calls.tsv + calls.altsigs.bin)\n"
        "  svp run    <normal.fq.gz> <tumor.fq.gz> <dir> [all opts]   (graph -> refine -> call; skips completed)\n");
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];
    int level = 5, margin = 10, min_mult = 1, threads = 16;
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
    else if (cmd == "call" && pos.size() == 1) run_call(pos[0]);
    else if (cmd == "run" && pos.size() == 3) {
        const std::string& dir = pos[2];
        if (exists(dir + "/graph.bin")) std::fprintf(stderr, "[skip] graph: %s/graph.bin exists\n", dir.c_str());
        else run_graph_reads(pos[0], pos[1], dir, level, margin, min_mult, threads);
        if (exists(dir + "/refined.bin")) std::fprintf(stderr, "[skip] refine: %s/refined.bin exists\n", dir.c_str());
        else run_refine(dir);
        if (exists(dir + "/calls.tsv")) std::fprintf(stderr, "[skip] call: %s/calls.tsv exists\n", dir.c_str());
        else run_call(dir);
    } else { usage(); return 2; }
    return 0;
}
