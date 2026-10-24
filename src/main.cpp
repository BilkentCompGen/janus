#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "cmd.h"
#include "encoding.h"

static const char* USAGE =
    "usage: janus <command> [options]\n"
    "\n"
    "  run       --samples sheet.tsv --ref PREFIX --model FILE --out DIR [--threads N] [--shallow N] [--clonal N]\n"
    "            the whole pipeline: one normal and one or more tumors of a patient (sheet columns: name, role, path);\n"
    "            without --model it stops after `call` (calls.vcf)\n"
    "\n"
    "  pack      --reads FILE --out x.svb [--threads N] [--shuffle]      reads (FASTQ, FASTQ.gz or BAM) packed two bits per base\n"
    "  extract   --reads x.svb [--threads N]                             LCP cores of every read, into x.svb.exc\n"
    "  count     --reads x.svb --out x.cnt [--normal n.cnt] [--support M] [--shallow N] [--clonal N] [--threads N]\n"
    "  graph     --counts n.cnt,t.cnt,... --out graph.gph --support M\n"
    "  prune     --graph graph.gph --reads n.svb,t.svb,... --out PREFIX --support M [--threads N]\n"
    "  select    --graph graph.gph --sites PREFIX.sites --reads n.svb,t.svb,... --out ev [--threads N]\n"
    "  refidx    --ref ref.fa --out ref.cidx [--threads N]\n"
    "  place     --pieces ev --index ref.cidx --ref ref.fa [--threads N]\n"
    "  lift      --paf local.paf --ref ref.fa --out aln.paf\n"
    "  clip      --pieces ev --paf aln.paf --out ev.gaps.fa\n"
    "  patch     --pieces ev --paf ev.gaps.raw.paf --aln aln.paf --out ev.paf\n"
    "  call      --pieces ev --reads n.svb,t.svb,... --ref ref.fa --graph graph.gph --vntr ref.trf.bed --out calls.vcf --support M [--threads N]\n"
    "  classify  --calls calls.vcf --hist t.cnt.hist --model FILE --out janus.vcf\n"
    "  joint     --calls calls.vcf --hist t1.cnt.hist,... --minimum M1,... --model FILE --out janus.vcf\n"
    "            [--raise 0.2] [--floor 0.1] [--mu 0.5] [--no-clade]   the acceptance rules of a joint run (see the paper)\n"
    "\n"
    "--support is the minimum read support (0 = chosen from the count spectrum); --shallow and --clonal are the\n"
    "count-stage depth thresholds (default 20); --threads defaults to 16.\n";

Args::Args(int argc, char** argv) : command(argc > 1 ? argv[1] : ""), reads(nullptr), out(nullptr), normal(nullptr), counts(nullptr), graph(nullptr), sites(nullptr), index(nullptr), ref(nullptr), pieces(nullptr), paf(nullptr), aln(nullptr), vntr(nullptr), model(nullptr), hist(nullptr), minimum(nullptr), calls(nullptr), samples(nullptr), threads(16), support(3), shallow(20), clonal(20), raise_(0.2), floor_(0.1), mu(0.5), shuffle(false), noclade(false) {
    struct { const char* name; const char** value; } strings[] = {{"--reads", &reads}, {"--out", &out}, {"--normal", &normal}, {"--counts", &counts}, {"--graph", &graph}, {"--sites", &sites}, {"--index", &index}, {"--ref", &ref}, {"--pieces", &pieces}, {"--paf", &paf}, {"--aln", &aln}, {"--vntr", &vntr}, {"--model", &model}, {"--hist", &hist}, {"--minimum", &minimum}, {"--calls", &calls}, {"--samples", &samples}};
    struct { const char* name; uint32_t* value; } numbers[] = {{"--support", &support}, {"--shallow", &shallow}, {"--clonal", &clonal}};
    struct { const char* name; double* value; } reals[] = {{"--raise", &raise_}, {"--floor", &floor_}, {"--mu", &mu}};
    for (int i = 2; i < argc; i++) {
        bool found = false;
        if (!std::strcmp(argv[i], "--shuffle")) { shuffle = true; continue; }
        if (!std::strcmp(argv[i], "--no-clade")) { noclade = true; continue; }
        if (i + 1 >= argc) { std::fprintf(stderr, "janus: %s needs a value\n", argv[i]); std::exit(2); }
        if (!std::strcmp(argv[i], "--threads")) { threads = std::atoi(argv[++i]); continue; }
        for (auto& s : strings)
            if (!std::strcmp(argv[i], s.name)) { *s.value = argv[++i]; found = true; break; }
        for (auto& n : numbers)
            if (!found && !std::strcmp(argv[i], n.name)) { *n.value = std::strtoul(argv[++i], nullptr, 10); found = true; break; }
        for (auto& r : reals)
            if (!found && !std::strcmp(argv[i], r.name)) { *r.value = std::strtod(argv[++i], nullptr); found = true; break; }
        if (!found) { std::fprintf(stderr, "janus: unknown option %s\n%s", argv[i], USAGE); std::exit(2); }
    }
}

const char* Args::need(const char* value, const char* option) const {
    if (!value) { std::fprintf(stderr, "janus %s: %s is required\n%s", command, option, USAGE); std::exit(2); }
    return value;
}

int main(int argc, char** argv) {
    Args a(argc, argv);
    const char* c = a.command;
    LCP_INIT();
    if (!std::strcmp(c, "run")) return cmd_run(a);
    if (!std::strcmp(c, "pack")) return cmd_pack(a);
    if (!std::strcmp(c, "extract")) return cmd_extract(a);
    if (!std::strcmp(c, "count")) return cmd_count(a);
    if (!std::strcmp(c, "graph")) return cmd_graph(a);
    if (!std::strcmp(c, "prune")) return cmd_prune(a);
    if (!std::strcmp(c, "select")) return cmd_select(a);
    if (!std::strcmp(c, "refidx")) return cmd_refidx(a);
    if (!std::strcmp(c, "place")) return cmd_place(a);
    if (!std::strcmp(c, "lift")) return cmd_lift(a);
    if (!std::strcmp(c, "clip")) return cmd_clip(a);
    if (!std::strcmp(c, "patch")) return cmd_patch(a);
    if (!std::strcmp(c, "call")) return cmd_call(a);
    if (!std::strcmp(c, "classify")) return cmd_classify(a);
    if (!std::strcmp(c, "joint")) return cmd_joint(a);
    std::fputs(USAGE, stderr);
    return 1;
}
