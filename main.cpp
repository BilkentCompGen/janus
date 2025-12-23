#include <iostream>
#include <string>
#include <map>
#include <any>
#include <filesystem>

#include "src/utils.h"
#include "src/graph.h"
#include "src/vcf.h"
#include "lib/lcptools/lps.h"
#include "lib/lcptools/encoding.h"

using namespace std;

Logger main_logger("main");

int main(int argc, char** argv) {
    LCP_INIT();

    ArgParser parser("Build LCP graph from reference (normal) and tumor FASTA files.");

    parser.add("normal", "--normal", "-n", "Reference genome FASTA file (normal)", "string", true);
    parser.add("tumor", "--tumor", "-t", "Tumor genome FASTA file", "string", true);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("output", "--output", "-o", "Output prefix (without extension)", "string", true);

    map<string, any> args = parser.parse(argc, argv);
    string normal_file = any_cast<string>(args["normal"]);
    string tumor_file = any_cast<string>(args["tumor"]);
    int depth = any_cast<int>(args["depth"]);
    string output_prefix = any_cast<string>(args["output"]);

    main_logger.debug("Normal (reference) file: " + normal_file);
    main_logger.debug("Tumor file: " + tumor_file);
    main_logger.debug("Depth: " + to_string(depth));
    main_logger.debug("Output prefix: " + output_prefix);

    namespace fs = std::filesystem;
    fs::path output_dir = "outputs";
    fs::create_directories(output_dir);
    main_logger.debug("Output directory: " + output_dir.string());

    LCPGraph graph;

    main_logger.debug("Processing reference genome (normal, color 0) - only chromosome 19...");
    FastaStreamer ref_streamer(normal_file, false);
    int processed_count = 0;
    int skipped_count = 0;
    for (const FastaPair& pair : ref_streamer) {
        string chrom = pair.first;
        if (chrom != "19") continue;

        processed_count++;
        main_logger.debug("Processing chromosome: " + chrom);

        graph.setSequence(0, chrom, pair.second);

        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, 0, pair.second.c_str(), chrom);
        free_lps(&lcp_string);
    }
    main_logger.debug("Processed " + to_string(processed_count) + " reference sequences (skipped " + to_string(skipped_count) + " non-chr19 sequences)");

    main_logger.debug("Processing tumor genome (color 1) - only chromosome 19...");
    FastaStreamer tumor_streamer(tumor_file, false);    
    for (const FastaPair& pair : tumor_streamer) {
        string chrom = pair.first;

        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
        if (chrom != "19") continue;

        main_logger.debug("Processing tumor chromosome: " + chrom);

        graph.setSequence(1, chrom, pair.second);

        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, 1, pair.second.c_str(), chrom);
        free_lps(&lcp_string);
    }
    main_logger.debug("Processed " + to_string(processed_count) + " tumor sequences (skipped " + to_string(skipped_count) + " non-chr19 sequences)");

    main_logger.debug("Graph built with " + to_string(graph.nodeCount()) + " nodes");

    main_logger.debug("Parsing VCF for deletion regions...");
    std::string vcf_path = "simulation/sim.vcf";
    auto deletions = parseVCFDeletions(vcf_path);
    main_logger.debug("Found " + to_string(deletions.size()) + " deletion regions in VCF");
    
    std::unordered_set<node_hash> deletion_node_hashes = findDeletionNodeHashes(deletions, graph);
    main_logger.debug("Found " + to_string(deletion_node_hashes.size()) + " node hashes overlapping with deletions");
    
    main_logger.debug("Compressing graph...");
    graph.compressGraph(&deletion_node_hashes);
    main_logger.debug("Graph compressed. Deletion nodes after compression: " + to_string(deletion_node_hashes.size()));

    main_logger.debug("Exporting graph to GFA and CSV...");
    graph.exportGFA(output_dir.string(), output_prefix);
    main_logger.debug("Graph exported");

    main_logger.debug("Done!");

    return 0;
}
