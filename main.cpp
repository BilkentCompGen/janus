#include <any>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>

#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"

#include "src/graph.h"
#include "src/ref.h"
#include "src/sv.h"
#include "src/utils.h"

using namespace std;

Logger main_logger("main");

void processGenome(const string& file, int color, uint32_t& sid, string& seq, LCPGraph& graph, const string& chr_str, int depth) {
    FastaStreamer streamer(file, false);
    for (const FastaPair& pair : streamer) {
        string chrom = pair.first;
        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
        if (chrom != chr_str) continue;

        sid = graph.seq_names.size();
        graph.seq_names.push_back(pair.first);
        seq = pair.second;

        struct lps lcp_string;
        init_lps(&lcp_string, seq.c_str(), seq.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, color, seq.c_str(), sid);
        free_lps(&lcp_string);
    }
}

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
    string ref_file = any_cast<string>(args["ref"]), normal_file = any_cast<string>(args["normal"]), tumor_file = any_cast<string>(args["tumor"]), output_prefix = any_cast<string>(args["output"]);
    int depth = any_cast<int>(args["depth"]);
    int chr = args.count("chr") ? any_cast<int>(args["chr"]) : 19;
    string chr_str = to_string(chr);

    namespace fs = std::filesystem;
    fs::path output_dir = "outputs";
    fs::path gfa_dir = output_dir / "gfa";
    fs::path vcf_dir = output_dir / "vcf";
    fs::create_directories(output_dir);
    fs::create_directories(gfa_dir);
    fs::create_directories(vcf_dir);

    RefIndex ref_index;
    LCPGraph graph;
    uint32_t normal_seq_id = 0, tumor_seq_id = 0;
    string ref_sequence, normal_sequence, tumor_sequence;

    // Parse reference into RefIndex
    FastaStreamer ref_streamer(ref_file, false);
    for (const FastaPair& pair : ref_streamer) {
        string chrom = pair.first;
        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
        if (chrom != chr_str) continue;

        ref_sequence = pair.second;
        struct lps lcp_string;
        init_lps(&lcp_string, ref_sequence.c_str(), ref_sequence.size());
        lps_deepen(&lcp_string, depth);
        ref_index.addString(&lcp_string, ref_sequence.c_str());
        free_lps(&lcp_string);
    }

    processGenome(normal_file, 0, normal_seq_id, normal_sequence, graph, chr_str, depth);  // Normal = bit 0 (0x1)
    processGenome(tumor_file, 1, tumor_seq_id, tumor_sequence, graph, chr_str, depth);     // Tumor = bit 1 (0x2)

    BubbleResult bubbles = graph.detectBubbles();
    auto extracted_bubbles = graph.extractBubbles(bubbles);

    // Call SVs and write VCF
    string vcf_file = (vcf_dir / (output_prefix + ".vcf")).string();
    ofstream vcf(vcf_file);
    vcf << "##fileformat=VCFv4.2\n"
        << "##source=svp\n"
        << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
        << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of variant\">\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    main_logger.debug("Calling SVs on " + to_string(extracted_bubbles.size()) + " bubbles");
    int sv_count = 0, anchored = 0;
    for (const auto& b : extracted_bubbles) {
        int64_t ref_offset = ref_index.anchorBubble(b, graph);
        if (ref_offset < 0) continue;
        anchored++;

        for (const auto& c : callSV(b, ref_offset)) {
            sv_count++;
            string ref_allele = c.normal_seq.empty() ? "." : c.normal_seq;
            string alt_allele = c.tumor_seq.empty() ? "." : c.tumor_seq;
            vcf << "chr" << chr_str << "\t" << c.ref_pos << "\t.\t" << ref_allele << "\t" << alt_allele << "\t.\tPASS\t" << "SVTYPE=" << c.svtype << ";SVLEN=" << c.length << "\n";
        }
    }
    vcf.close();
    main_logger.debug("Anchored " + to_string(anchored) + "/" + to_string(extracted_bubbles.size()) + " bubbles, called " + to_string(sv_count) + " SVs → " + vcf_file);

    string gfa_file = (gfa_dir / (output_prefix + "_graph.gfa")).string();
    main_logger.debug("Writing GFA to " + gfa_file);
    graph.writeGFA(gfa_file);

    main_logger.debug("Done!");
    return 0;
}
