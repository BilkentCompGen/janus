// svp: 2-colour LCP-core graph built FROM READS (normal vs tumour reads), then the somatic material
// carved out of those reads as junction contigs. READS ONLY -- no assembly, and no reference until
// the contigs are handed to an aligner.
//
//   svp run       <normal> <tumour> <ref.fa> <dir>   THE PIPELINE, end to end
//
//   svp graph     <normal> <tumour> <dir>      reads           -> dir/graph.bin
//   svp bubbles   <graph.bin> <normal> <tumour> graph + reads  -> dir/drop.sigs (sub-SV cores)
//   svp junctions <graph.bin> <reads> <out.fa> graph + reads   -> junction contigs
//   svp call      <ctg.bam> <normal.bam|-> <out.vcf> aligned contigs -> breakpoints
//   svp pick      <ctg.fa> <calls.vcf> <normal> <out.fq> the normal reads worth aligning
//   svp pack      <out.svb>                    FASTQ (stdin)   -> 2-bit .svb
//   svp dump      <in.svb>                     .svb            -> FASTQ (stdout)
//
// Every stage takes .svb or uncompressed FASTQ, so the pipeline needs no shell to convert between
// them. The only outside programs are the aligner and samtools, which svp is not and does not
// reimplement; `run` invokes them and every default it passes is the one measured in
// docs/RESULTS.md.
#include "bubbles.h"
#include "call.h"
#include "graph.h"
#include "junctions.h"
#include "pick.h"
#include "svb.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <vector>

extern "C" {
#include "encoding.h"
}

static bool exists(const std::string& p) { FILE* f = std::fopen(p.c_str(), "rb"); if (f) { std::fclose(f); return true; } return false; }

static void usage() {
    std::fprintf(stderr,
        "usage:\n"
        "  svp run       <normal> <tumour> <ref.fa> <dir> [--threads 16] [--level 5] [--margin 10] [--min-mult 3]\n"
        "  svp graph     <normal> <tumour> <dir> [--level 5] [--margin 10] [--min-mult 3] [--threads 16]\n"
        "  svp junctions <graph.bin> <reads> <out.fa> [threads=16] [pad_cores=50] [per_chunk=5] [min_flank=300] [min_vaf=0]\n"
        "  svp bubbles   <graph.bin> <normal> <tumour> [threads=16] [k=21] [w=22] [maxlen=40] [max_dlen=12] [max_unmatched=20]\n"
        "  svp call      <ctg.bam> <normal.bam|-> <out.vcf> [min_mapq=10] [min_aln=0] [min_sv=50] [cluster=50] [min_supp=4]\n"
        "  svp pick      <ctg.fa> <calls.vcf> <normal reads> <out.fastq> [threads=16] [min_shared=5]\n"
        "  svp pack      <out.svb> < reads.fastq\n"
        "  svp dump      <in.svb> > reads.fastq\n");
}

// ---- `run`: the pipeline ---------------------------------------------------------------------
//
// Each stage is skipped when its output is already there, so a kill and a restart pick up where this
// left off. Two of them are the aligner and samtools; the rest is svp.
//
//   graph.bin        the 2-colour core graph, from both samples' reads
//   drop.sigs        the somatic cores that are sub-SV edits and must never become a contig
//   ctg.fa           one junction contig per somatic chunk per read, 50 cores of anchor each side
//   ctg.bam          those contigs against the reference
//   candidates.vcf   `call` with no normal: the loci the germline question will be asked about
//   normal_at_       the normal reads sharing cores with those candidates' germline anchors --
//     candidates       1.56% of the sample on hcc1954
//   normal.bam       ONLY those, against the reference
//   calls.vcf        `call` again, now able to ask the normal at each locus
//
// The normal is not aligned whole. That decision is worth 1168 false calls on hcc1954 and needs full
// depth to make (3x leaves 290 false calls, 5.5x 172, 11x 136, 44x 116) -- but it is only ever asked
// at the candidate loci, and the graph can say which normal reads cover those. Measured: 177k picked
// reads give 0.887/0.858 against 0.891/0.858 for the whole 8.35M-read normal, in fourteen minutes
// of alignment instead of two hours.

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

// minimap2 cannot read .svb, and it is the only thing in the pipeline that cannot, so this is where
// the conversion happens -- streamed, never written out.
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

static int run_pipeline(const std::string& normal, const std::string& tumour, const std::string& ref,
                        const std::string& dir, int level, int margin, int min_mult, int threads) {
    const std::string self = self_path();
    const std::string graph = dir + "/graph.bin", drop = dir + "/drop.sigs", ctg = dir + "/ctg.fa";
    const std::string cbam = dir + "/ctg.bam", nbam = dir + "/normal.bam", vcf = dir + "/calls.vcf";
    const std::string cand = dir + "/candidates.vcf", nfq = dir + "/normal_at_candidates.fastq";
    const std::string T = std::to_string(threads);

    if (exists(graph)) std::fprintf(stderr, "[run] atlandi graph: %s var\n", graph.c_str());
    else run_graph_reads(normal, tumour, dir, level, margin, min_mult, threads);

    if (exists(drop)) std::fprintf(stderr, "[run] atlandi bubbles: %s var\n", drop.c_str());
    else {
        std::vector<std::string> a{"bubbles", graph, normal, tumour, T};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_bubbles((int)v.size(), v.data())) return rc;
    }

    if (exists(ctg)) std::fprintf(stderr, "[run] atlandi junctions: %s var\n", ctg.c_str());
    else {
        std::vector<std::string> a{"junctions", graph, tumour, ctg, T};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_junctions((int)v.size(), v.data())) return rc;
    }

    if (exists(cbam)) std::fprintf(stderr, "[run] atlandi align: %s var\n", cbam.c_str());
    else {
        if (sh(align(ctg, ref, cbam, dir + "/stc", threads, self))) return 1;
        if (sh("samtools index -@ " + T + " '" + cbam + "'")) return 1;
    }

    // The germline check needs the normal AT the candidate loci and nowhere else, so the candidates
    // are made first, without it, and they are what says which normal reads to align.
    if (exists(cand)) std::fprintf(stderr, "[run] atlandi aday call: %s var\n", cand.c_str());
    else {
        std::vector<std::string> a{"call", cbam, "-", cand};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_call((int)v.size(), v.data())) return rc;
    }

    if (exists(nfq)) std::fprintf(stderr, "[run] atlandi pick: %s var\n", nfq.c_str());
    else {
        std::vector<std::string> a{"pick", ctg, cand, normal, nfq, T};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_pick((int)v.size(), v.data())) return rc;
    }

    if (exists(nbam)) std::fprintf(stderr, "[run] atlandi normal align: %s var\n", nbam.c_str());
    else {
        if (sh(align(nfq, ref, nbam, dir + "/stn", threads, self))) return 1;
        if (sh("samtools index -@ " + T + " '" + nbam + "'")) return 1;
    }

    if (exists(vcf)) std::fprintf(stderr, "[run] atlandi call: %s var\n", vcf.c_str());
    else {
        std::vector<std::string> a{"call", cbam, nbam, vcf};
        std::vector<char*> v; for (auto& s : a) v.push_back(&s[0]);
        if (int rc = cmd_call((int)v.size(), v.data())) return rc;
    }
    std::fprintf(stderr, "[run] bitti -> %s\n", vcf.c_str());
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 2; }
    const std::string cmd = argv[1];

    // These take their arguments positionally, with no --flags to strip, so they are dispatched
    // before the option parser runs. junctions needs the LCP alphabet; pack/dump do not.
    if (cmd == "junctions") { LCP_INIT(); return cmd_junctions(argc - 1, argv + 1); }
    if (cmd == "bubbles")   { LCP_INIT(); return cmd_bubbles(argc - 1, argv + 1); }
    if (cmd == "call")      return cmd_call(argc - 1, argv + 1);
    if (cmd == "pick")      { LCP_INIT(); return cmd_pick(argc - 1, argv + 1); }
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

    if (cmd == "run" && pos.size() == 4)
        return run_pipeline(pos[0], pos[1], pos[2], pos[3], level, margin, min_mult, threads);
    if (cmd == "graph" && pos.size() == 3) run_graph_reads(pos[0], pos[1], pos[2], level, margin, min_mult, threads);
    else { usage(); return 2; }
    return 0;
}
