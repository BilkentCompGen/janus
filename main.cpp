#include <any>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <map>
#include <string>

#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"

#include "src/align.h"
#include "src/graph.h"
#include "src/ref.h"
#include "src/utils.h"

using namespace std;

Logger main_logger("main");

int main(int argc, char** argv) {
    LCP_INIT();

    ArgParser parser("Build LCP graph from reference, normal, and tumor FASTA files.");

    parser.add("ref", "--ref", "-r", "Reference genome FASTA file", "string", true);
    parser.add("normal", "--normal", "-n", "Normal genome FASTA file", "string", true);
    parser.add("tumor", "--tumor", "-t", "Tumor genome FASTA file", "string", true);
    parser.add("vcf", "--vcf", "-v", "VCF file with ground truth deletions", "string", false);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("output", "--output", "-o", "Output prefix (without extension)", "string", true);

    map<string, any> args = parser.parse(argc, argv);
    string ref_file = any_cast<string>(args["ref"]);
    string normal_file = any_cast<string>(args["normal"]);
    string tumor_file = any_cast<string>(args["tumor"]);
    int depth = any_cast<int>(args["depth"]);
    string output_prefix = any_cast<string>(args["output"]);
    
    string vcf_file = "";
    if (args.count("vcf")) {
        vcf_file = any_cast<string>(args["vcf"]);
    }

    main_logger.debug("Reference file: " + ref_file);
    main_logger.debug("Normal file: " + normal_file);
    main_logger.debug("Tumor file: " + tumor_file);
    if (!vcf_file.empty()) main_logger.debug("VCF file: " + vcf_file);
    main_logger.debug("Depth: " + to_string(depth));
    main_logger.debug("Output prefix: " + output_prefix);

    namespace fs = std::filesystem;
    fs::path output_dir = "outputs";
    fs::create_directories(output_dir);
    main_logger.debug("Output directory: " + output_dir.string());

    LCPGraph graph;
    uint32_t ref_seq_id = 0, normal_seq_id = 0, tumor_seq_id = 0;
    string ref_sequence, normal_sequence, tumor_sequence;

    // Process sequences
    auto process = [&](const string& file, int color, uint32_t& sid, string& seq) {
        FastaStreamer streamer(file, false);
        for (const FastaPair& pair : streamer) {
            string chrom = pair.first;
            size_t hap_pos = chrom.find("_hap");
            if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
            if (chrom != "19") continue;
            sid = graph.seq_names.size();
            graph.seq_names.push_back(pair.first);
            seq = pair.second;
            struct lps lcp_string;
            init_lps(&lcp_string, seq.c_str(), seq.size());
            lps_deepen(&lcp_string, depth);
            graph.addString(&lcp_string, color, seq.c_str(), sid);
            free_lps(&lcp_string);
        }
    };

    process(ref_file, 0, ref_seq_id, ref_sequence);
    process(normal_file, 1, normal_seq_id, normal_sequence);
    process(tumor_file, 2, tumor_seq_id, tumor_sequence);

    main_logger.debug("Graph built with " + to_string(graph.nodeCount()) + " nodes");

    // Detect and mark
    BubbleResult bubbles = graph.detectBubbles();
    graph.markBubbles(bubbles);
    if (!vcf_file.empty()) {
        auto deletions = parseVCFDeletions(vcf_file);
        graph.markDeletions(findDeletionNodeHashes(deletions, graph, ref_seq_id));
    }

    // Alignment loop
    auto extracted_bubbles = graph.extractBubbles(bubbles, ref_seq_id, normal_seq_id, tumor_seq_id);
    main_logger.debug("Extracted " + to_string(extracted_bubbles.size()) + " bubbles");
    
    vector<VariantCall> all_variants;
    for (const auto& b : extracted_bubbles) {
        if (b.ref_s == UINT64_MAX || b.norm_s == UINT64_MAX || b.tum_s == UINT64_MAX) continue;
        if (b.norm_s >= b.norm_e || b.norm_e > normal_sequence.size()) continue;
        if (b.tum_s >= b.tum_e || b.tum_e > tumor_sequence.size()) continue;

        string n_seq = normal_sequence.substr(b.norm_s, b.norm_e - b.norm_s);
        string t_seq = tumor_sequence.substr(b.tum_s, b.tum_e - b.tum_s);
        if (n_seq.size() > 500000) continue; 

        auto vars = alignAndCall(n_seq, t_seq, b.ref_s + 1, "19");
        for (auto& v : vars) all_variants.push_back(std::move(v));
    }
    
    exportAlignedVCF((output_dir / (output_prefix + "_variants.vcf")).string(), all_variants);

    map<uint32_t, uint32_t> color_map = {{0x01, 0x808080}, {0x02, 0x0000FF}, {0x04, 0xFF0000}, {0x08, 0x00FF00}, {0x10, 0xFF8C00}, {0x20, 0xFFFF00}};
    graph.exportGFA(output_dir.string(), output_prefix + "_raw", &color_map);
    graph.compressGraph();
    graph.exportGFA(output_dir.string(), output_prefix + "_compressed", &color_map);

    main_logger.debug("Done!");
    return 0;
}
