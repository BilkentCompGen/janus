#include <any>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <map>
#include <string>

#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"

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
    uint32_t ref_seq_id = 0;
    uint32_t normal_seq_id = 0;
    uint32_t tumor_seq_id = 0;
    
    // Store sequences for later alignment
    string ref_sequence;
    string normal_sequence;
    string tumor_sequence;

    // Process reference genome (color 0)
    main_logger.debug("Processing reference genome (color 0) - only chromosome 19...");
    FastaStreamer ref_streamer(ref_file, false);
    for (const FastaPair& pair : ref_streamer) {
        string chrom = pair.first;
        if (chrom != "19") continue;

        main_logger.debug("Processing reference chromosome: " + chrom);

        ref_seq_id = graph.seq_names.size();
        graph.seq_names.push_back(pair.first);
        ref_sequence = pair.second;  // Store for alignment

        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, 0, pair.second.c_str(), ref_seq_id);
        free_lps(&lcp_string);
    }
    main_logger.debug("Reference genome processed. Length: " + to_string(ref_sequence.size()));

    // Process normal genome (color 1)
    main_logger.debug("Processing normal genome (color 1) - only chromosome 19...");
    FastaStreamer normal_streamer(normal_file, false);
    for (const FastaPair& pair : normal_streamer) {
        string chrom = pair.first;
        if (chrom != "19") continue;

        main_logger.debug("Processing normal chromosome: " + chrom);

        normal_seq_id = graph.seq_names.size();
        graph.seq_names.push_back(pair.first);
        normal_sequence = pair.second;  // Store for alignment

        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, 1, pair.second.c_str(), normal_seq_id);
        free_lps(&lcp_string);
    }
    main_logger.debug("Normal genome processed. Length: " + to_string(normal_sequence.size()));

    // Process tumor genome (color 2)
    main_logger.debug("Processing tumor genome (color 2) - only chromosome 19...");
    FastaStreamer tumor_streamer(tumor_file, false);
    for (const FastaPair& pair : tumor_streamer) {
        string chrom = pair.first;

        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
        if (chrom != "19") continue;

        main_logger.debug("Processing tumor chromosome: " + chrom);

        tumor_seq_id = graph.seq_names.size();
        graph.seq_names.push_back(pair.first);
        tumor_sequence = pair.second;  // Store for alignment

        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, 2, pair.second.c_str(), tumor_seq_id);
        free_lps(&lcp_string);
    }
    main_logger.debug("Tumor genome processed. Length: " + to_string(tumor_sequence.size()));

    main_logger.debug("Graph built with " + to_string(graph.nodeCount()) + " nodes");

    // Detect bubbles
    main_logger.debug("Detecting bubbles...");
    BubbleResult bubbles = graph.detectBubbles();
    graph.markBubbles(bubbles);

    // Mark deletions from VCF if provided
    if (!vcf_file.empty()) {
        main_logger.debug("Parsing VCF for deletion regions...");
        auto deletions = parseVCFDeletions(vcf_file);
        main_logger.debug("Found " + to_string(deletions.size()) + " deletion regions in VCF");

        std::unordered_set<node_hash> deletion_node_hashes = findDeletionNodeHashes(deletions, graph, ref_seq_id);
        main_logger.debug("Found " + to_string(deletion_node_hashes.size()) + " nodes overlapping with deletions");
        
        graph.markDeletions(deletion_node_hashes);
    }

    // Extended color map for visualization
    // Bits: 0=Ref, 1=Normal, 2=Tumor, 3=BubbleStart, 4=BubbleEnd, 5=Deletion
    map<uint32_t, uint32_t> color_map = {
        // Basic colors (no special markers)
        {0x01, 0x808080},  // 000001: Ref only -> Gray
        {0x02, 0x0000FF},  // 000010: Normal only -> Blue
        {0x04, 0xFF0000},  // 000100: Tumor only -> Red
        {0x03, 0x00CED1},  // 000011: Ref + Normal -> Dark Cyan
        {0x05, 0xFF00FF},  // 000101: Ref + Tumor -> Magenta
        {0x06, 0x800080},  // 000110: Normal + Tumor -> Purple
        {0x07, 0x000000},  // 000111: Ref + Normal + Tumor -> Black
        
        // Bubble Start (bit 3 set) - Green variants
        {0x08, 0x00FF00},  // BubbleStart only -> Bright Green
        {0x09, 0x228B22},  // BubbleStart + Ref -> Forest Green
        {0x0A, 0x32CD32},  // BubbleStart + Normal -> Lime Green
        {0x0B, 0x006400},  // BubbleStart + Ref + Normal -> Dark Green
        {0x0C, 0x7CFC00},  // BubbleStart + Tumor -> Lawn Green
        {0x0D, 0x90EE90},  // BubbleStart + Ref + Tumor -> Light Green
        {0x0E, 0x98FB98},  // BubbleStart + Normal + Tumor -> Pale Green
        {0x0F, 0x00FF7F},  // BubbleStart + All -> Spring Green
        
        // Bubble End (bit 4 set) - Orange variants
        {0x10, 0xFF8C00},  // BubbleEnd only -> Dark Orange
        {0x11, 0xFFA500},  // BubbleEnd + Ref -> Orange
        {0x12, 0xFF7F50},  // BubbleEnd + Normal -> Coral
        {0x13, 0xE9967A},  // BubbleEnd + Ref + Normal -> Dark Salmon
        {0x14, 0xFF6347},  // BubbleEnd + Tumor -> Tomato
        {0x15, 0xFA8072},  // BubbleEnd + Ref + Tumor -> Salmon
        {0x16, 0xF08080},  // BubbleEnd + Normal + Tumor -> Light Coral
        {0x17, 0xCD853F},  // BubbleEnd + All -> Peru
        
        // Deletion (bit 5 set) - Yellow variants
        {0x20, 0xFFFF00},  // Deletion only -> Yellow
        {0x21, 0xFFD700},  // Deletion + Ref -> Gold
        {0x22, 0xDAA520},  // Deletion + Normal -> Goldenrod
        {0x23, 0xF0E68C},  // Deletion + Ref + Normal -> Khaki
        {0x24, 0xBDB76B},  // Deletion + Tumor -> Dark Khaki
        {0x25, 0xEEE8AA},  // Deletion + Ref + Tumor -> Pale Goldenrod
        {0x26, 0xFFE4B5},  // Deletion + Normal + Tumor -> Moccasin
        {0x27, 0xFFDAB9},  // Deletion + All -> Peach Puff
    };

    // Extract bubble structures with coordinates
    main_logger.debug("Extracting bubble structures...");
    auto extracted_bubbles = graph.extractBubbles(bubbles, ref_seq_id, normal_seq_id, tumor_seq_id);
    main_logger.debug("Extracted " + to_string(extracted_bubbles.size()) + " bubbles");
    
    // Collect deletion regions from bubbles (simplified - no alignment)
    main_logger.debug("Collecting deletion regions from bubbles...");
    vector<VariantRegion> bubble_variants;
    
    for (const auto& bubble : extracted_bubbles) {
        uint64_t normal_start = bubble.ref_start;
        uint64_t normal_end = bubble.ref_end;
        
        // Validate coordinates
        if (normal_start >= normal_end) continue;
        if (normal_end > normal_sequence.size()) continue;
        
        // Calculate deletion size based on coordinate difference
        uint64_t normal_len = normal_end - normal_start;
        uint64_t tumor_len = 0;
        if (bubble.tumor_start < bubble.tumor_end) {
            tumor_len = bubble.tumor_end - bubble.tumor_start;
        }
        
        // If normal region is larger than tumor region, it's a deletion
        if (normal_len > tumor_len) {
            bubble_variants.push_back({"19", normal_start + 1, normal_end, "DEL"});
        }
    }
    
    main_logger.debug("Found " + to_string(bubble_variants.size()) + " deletions from bubbles");
    
    // Export VCF
    string vcf_output = (output_dir / (output_prefix + "_variants.vcf")).string();
    main_logger.debug("Exporting variants to: " + vcf_output);
    exportVCF(vcf_output, bubble_variants);
    main_logger.debug("VCF exported");

    // Export raw (uncompressed) graph first
    main_logger.debug("Exporting raw graph to GFA and CSV...");
    graph.exportGFA(output_dir.string(), output_prefix + "_raw", &color_map);
    main_logger.debug("Raw graph exported (" + to_string(graph.nodeCount()) + " nodes)");

    // Compress graph for visualization (merge linear chains with same color)
    main_logger.debug("Compressing graph for visualization...");
    graph.compressGraph();
    main_logger.debug("Graph compressed to " + to_string(graph.nodeCount()) + " nodes");

    // Export compressed graph
    main_logger.debug("Exporting compressed graph to GFA and CSV...");
    graph.exportGFA(output_dir.string(), output_prefix + "_compressed", &color_map);
    main_logger.debug("Compressed graph exported");

    main_logger.debug("Done!");

    return 0;
}
