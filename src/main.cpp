// svp (lcp-cortex): somatic SV caller over a 2-colour LCP-core graph built FROM READS
// (normal vs tumor HiFi reads). READS ONLY; no assembly in the pipeline. Runs as three
// resumable, checkpointed stages so the expensive LCP pass is paid once:
//   svp lcp     <normal.fq.gz> <tumor.fq.gz> <dir>   -> dir/counts.bin (post-prune sigs +
//                                                       inline bytes), dir/edges.bin
//   svp graph   <dir>                                -> dir/graph.bin   (colour + edges)
//   svp bubbles <dir>                                -> dir/bubbles.tsv
//   svp run     <normal.fq.gz> <tumor.fq.gz> <dir>   -> all, skipping completed stages
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
        "  svp lcp     <normal.fq.gz> <tumor.fq.gz> <dir> [--level 5] [--flank 10] [--margin 10] [--min-mult 3] [--threads 16]\n"
        "  svp graph   <dir>\n"
        "  svp modea   <dir>\n"
        "  svp modeb   <dir>\n"
        "  svp posmap  <ref.fasta> <dir> [--level 5] [--flank 10]\n"
        "  svp anchor  <dir>\n"
        "  svp vcf     <dir> [--ref <ref.fasta>] [--out <calls.vcf>] [--dup-flags <path>] [--inv-flags <path>]\n"
        "  svp bubbles <dir> [--out <dir>/bubbles.tsv]\n"
        "  svp run     <normal.fq.gz> <tumor.fq.gz> <dir> [--ref <ref.fasta>] [all opts]   (skips completed stages)\n");
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];
    int level = 5, flank = 10, margin = 10, min_mult = 3, threads = 16;
    std::string out, ref, dup_flags, inv_flags;
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
        else if (a == "--out") out = v;
        else if (a == "--ref") ref = v;
        else if (a == "--dup-flags") dup_flags = v;
        else if (a == "--inv-flags") inv_flags = v;
        else { std::fprintf(stderr, "[FATAL] arg: %s\n", a.c_str()); std::abort(); }
    }
    LCP_INIT();

    if (cmd == "lcp" && pos.size() == 3) run_lcp(pos[0], pos[1], pos[2], level, flank, margin, min_mult, threads);
    else if (cmd == "graph" && pos.size() == 1) run_graph(pos[0]);
    else if (cmd == "modea" && pos.size() == 1) run_modea(pos[0]);
    else if (cmd == "modeb" && pos.size() == 1) run_modeb(pos[0]);
    else if (cmd == "posmap" && pos.size() == 2) run_posmap(pos[1], pos[0], level, flank);
    else if (cmd == "anchor" && pos.size() == 1) run_anchor(pos[0]);
    else if (cmd == "vcf" && pos.size() == 1) run_vcf(pos[0], ref, out.empty() ? pos[0] + "/calls.vcf" : out, dup_flags, inv_flags);
    else if (cmd == "bubbles" && pos.size() == 1) run_bubbles(pos[0], out.empty() ? pos[0] + "/bubbles.tsv" : out);
    else if (cmd == "run" && pos.size() == 3) {
        const std::string& dir = pos[2];
        if (exists(dir + "/counts.bin") && exists(dir + "/edges.bin")) std::fprintf(stderr, "[skip] lcp: %s/counts.bin exists\n", dir.c_str());
        else run_lcp(pos[0], pos[1], dir, level, flank, margin, min_mult, threads);
        if (exists(dir + "/graph.bin")) std::fprintf(stderr, "[skip] graph: %s/graph.bin exists\n", dir.c_str());
        else run_graph(dir);
        if (exists(dir + "/modea.bin")) std::fprintf(stderr, "[skip] modea: %s/modea.bin exists\n", dir.c_str());
        else run_modea(dir);
        if (exists(dir + "/modeb.bin")) std::fprintf(stderr, "[skip] modeb: %s/modeb.bin exists\n", dir.c_str());
        else run_modeb(dir);
        if (!ref.empty()) {
            if (exists(dir + "/sig_ref_map.bin")) std::fprintf(stderr, "[skip] posmap: %s/sig_ref_map.bin exists\n", dir.c_str());
            else run_posmap(dir, ref, level, flank);
            if (exists(dir + "/bubble_pos.bin")) std::fprintf(stderr, "[skip] anchor: %s/bubble_pos.bin exists\n", dir.c_str());
            else run_anchor(dir);
            if (exists(dir + "/calls.vcf")) std::fprintf(stderr, "[skip] vcf: %s/calls.vcf exists\n", dir.c_str());
            else run_vcf(dir, ref, dir + "/calls.vcf", dup_flags, inv_flags);
        }
        run_bubbles(dir, out.empty() ? dir + "/bubbles.tsv" : out);
    } else { usage(); return 2; }
    return 0;
}
