#include <iostream>
#include <string>
#include <map>
#include <any>
#include <filesystem>

#include "src/utils.h"
#include "src/graph.h"
#include "src/sv.h"
#include "src/analysis.h"
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

    // Graph and reference
    LCPGraph graph;
    LinearReference reference;
    CollisionTracker collision_tracker;

    main_logger.debug("Processing reference genome (normal, color 0) - only chromosome 19...");
    FastaStreamer ref_streamer(normal_file, false);
    int processed_count = 0;
    int skipped_count = 0;
    for (const FastaPair& pair : ref_streamer) {
        string chrom = pair.first;
        if (!chrom.empty() && chrom[0] == '>') {
            chrom = chrom.substr(1);
        }
        // Remove anything after first space
        size_t space_pos = chrom.find(' ');
        if (space_pos != string::npos) {
            chrom = chrom.substr(0, space_pos);
        }
        
        // Only process chromosome 19 (exact match, skip alt/random contigs)
        if (chrom != "19") {
            skipped_count++;
            continue;
        }

        processed_count++;
        main_logger.debug("Processing chromosome: " + chrom);

        // Store global sequence for this chromosome
        collision_tracker.setSequence(chrom, pair.second);

        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, 0, pair.second.c_str());

        for (int j = 0; j < lcp_string.size; ++j) {
            node_hash_t hash = hash_core(&lcp_string.cores[j], pair.second.c_str());
            uint64_t start = lcp_string.cores[j].start;
            uint64_t end = lcp_string.cores[j].end;
            reference.addCore(hash, chrom, start, end);
            collision_tracker.addHash(hash, chrom, start, end);
        }
        free_lps(&lcp_string);
    }
    main_logger.debug("Processed " + to_string(processed_count) + " reference sequences (skipped " + to_string(skipped_count) + " non-chr19 sequences)");

    reference.buildIndex();
    main_logger.debug("Linear reference built with " + to_string(reference.size()) + " positions");

    main_logger.debug("Analyzing hash collisions...");
    collision_tracker.generateTSVReport(output_dir, output_prefix);

    main_logger.debug("Processing tumor genome (color 1) - only chromosome 19...");
    FastaStreamer tumor_streamer(tumor_file, false);    
    for (const FastaPair& pair : tumor_streamer) {
        // Extract chromosome from FASTA header
        string chrom = pair.first;
        if (!chrom.empty() && chrom[0] == '>') {
            chrom = chrom.substr(1);
        }
        // Remove anything after first space
        size_t space_pos = chrom.find(' ');
        if (space_pos != string::npos) {
            chrom = chrom.substr(0, space_pos);
        }
        // Remove _hapA or _hapB suffix if present
        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) {
            chrom = chrom.substr(0, hap_pos);
        }
        
        // Only process chromosome 19 (exact match)
        if (chrom != "19") {
            continue;
        }

        main_logger.debug("Processing tumor chromosome: " + chrom);

        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, 1, pair.second.c_str());
        free_lps(&lcp_string);
    }
    main_logger.debug("Processed " + to_string(processed_count) + " tumor sequences (skipped " + to_string(skipped_count) + " non-chr19 sequences)");

    int num_tumor_samples = 1;
    main_logger.debug("Graph built with " + to_string(graph.nodeCount()) + " nodes");

    main_logger.debug("Compressing graph...");
    graph.compressGraph();
    main_logger.debug("Graph compressed");

    main_logger.debug("Exporting graph to GFA and CSV...");
    graph.exportGFA(output_dir.string(), output_prefix);
    main_logger.debug("Graph exported");

    main_logger.debug("Detecting bubbles...");
    std::vector<BubblePathReport> tumor_paths_report;
    std::vector<BubblePathReport> normal_paths_report;
    auto bubbles = detectBubbles(graph, reference, &tumor_paths_report, &normal_paths_report);
    main_logger.debug("Found " + to_string(bubbles.size()) + " bubbles");

    main_logger.debug("Exporting bubbles...");
    exportBubblesTSV(bubbles, output_dir.string(), output_prefix);
    exportBubblesDetailedTSV(bubbles, graph, output_dir.string(), output_prefix);

    // Export all tumor-only and normal-only paths for debugging
    std::vector<BubblePathReport> all_paths;
    all_paths.reserve(tumor_paths_report.size() + normal_paths_report.size());
    all_paths.insert(all_paths.end(), tumor_paths_report.begin(), tumor_paths_report.end());
    all_paths.insert(all_paths.end(), normal_paths_report.begin(), normal_paths_report.end());
    exportBubblePathsTSV(all_paths, output_dir.string(), output_prefix);

    main_logger.debug("Bubbles exported");

    main_logger.debug("Detecting structural variants...");
    auto sv_calls = detectDeletions(graph, reference, num_tumor_samples);

    int deletion_count = 0;
    for (const auto& sv : sv_calls) {
        if (sv.type == SVType::DELETION) deletion_count++;
    }

    main_logger.debug("Found " + to_string(sv_calls.size()) + " structural variants:");
    main_logger.debug("  Deletions: " + to_string(deletion_count));

    fs::path tsv_dir = output_dir / "tsv";
    fs::create_directories(tsv_dir);
    fs::path sv_filename = tsv_dir / (output_prefix + ".tsv");
    exportTSV(sv_calls, sv_filename.string());
    main_logger.debug("SVs saved to " + sv_filename.string());

    main_logger.debug("Done!");

    return 0;
}
