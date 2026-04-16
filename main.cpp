#include <any>
#include <filesystem>
#include <map>
#include <string>

#include "lib/lcptools/lps.h"

#include "src/graph.h"
#include "src/pipeline.h"
#include "src/ref.h"
#include "src/sv.h"
#include "src/utils.h"
#include "src/vcf.h"

using namespace std;

int main(int argc, char** argv) {
    LCP_INIT();
    ArgParser parser("Build LCP graph from reference, normal, and tumor FASTA files.");
    parser.add("ref", "--ref", "-r", "Reference genome FASTA file", "string", true);
    parser.add("normal", "--normal", "-n", "Normal genome FASTA file", "string", true);
    parser.add("tumor", "--tumor", "-t", "Tumor genome FASTA file", "string", true);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("chr", "--chr", "-c", "Chromosome to process", "int", false);
    parser.add("output", "--output", "-o", "Output prefix (without extension)", "string", true);

    map<string, any> args = parser.parse(argc, argv);
    string ref_file = any_cast<string>(args["ref"]);
    string normal_file = any_cast<string>(args["normal"]);
    string tumor_file = any_cast<string>(args["tumor"]);
    string output_prefix = any_cast<string>(args["output"]);
    int depth = any_cast<int>(args["depth"]);
    int chr = args.count("chr") ? any_cast<int>(args["chr"]) : 19;
    string chr_str = to_string(chr);

    namespace fs = std::filesystem;
    fs::path output_dir = "outputs";
    fs::path gfa_dir = output_dir / "gfa";
    fs::path vcf_dir = output_dir / "vcf";
    fs::create_directories(gfa_dir);
    fs::create_directories(vcf_dir);

    RefIndex ref_index;
    LCPGraph graph;
    uint32_t normal_seq_id = 0, tumor_seq_id = 0;
    string ref_sequence, normal_sequence, tumor_sequence;

    auto ref_hashes = buildRefHashes(ref_file, chr_str, depth, ref_index, ref_sequence);
    auto rcref_hashes = buildRcrefHashes(ref_sequence, depth);

    processGenome(normal_file, 0, normal_seq_id, normal_sequence, graph, chr_str, depth);
    processGenome(tumor_file, 1, tumor_seq_id, tumor_sequence, graph, chr_str, depth);

    auto large_inv_calls = detectLargeINVs(tumor_sequence, depth, ref_hashes, rcref_hashes, ref_index);

    BubbleResult boundaries = graph.detectBubbles();
    auto calls = extractAndCallSVs(graph, ref_index, boundaries, large_inv_calls);
    calls.insert(calls.end(), large_inv_calls.begin(), large_inv_calls.end());

    string vcf_file = (vcf_dir / (output_prefix + ".vcf")).string();
    writeVCF(vcf_file, calls, ref_sequence, chr_str);

    string gfa_file = (gfa_dir / (output_prefix + "_graph.gfa")).string();
    graph.writeGFA(gfa_file);

    return 0;
}
