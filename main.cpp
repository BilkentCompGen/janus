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

        struct lps lcp_string;
        init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, 0);

        for (int j = 0; j < lcp_string.size; ++j) {
            ublock hash = hash_core(&lcp_string.cores[j]);
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
        graph.addString(&lcp_string, 1);
        free_lps(&lcp_string);
    }
    main_logger.debug("Processed " + to_string(processed_count) + " tumor sequences (skipped " + to_string(skipped_count) + " non-chr19 sequences)");

    int num_tumor_samples = 1;
    main_logger.debug("Graph built with " + to_string(graph.nodeCount()) + " nodes");

    main_logger.debug("Exporting graph to GFA and CSV...");
    graph.exportGFA(output_dir.string(), output_prefix);
    main_logger.debug("Graph exported to " + (output_dir / "gfa" / (output_prefix + ".gfa")).string() + " and " + (output_dir / "csv" / (output_prefix + ".csv")).string());

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
