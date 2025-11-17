#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <any>
#include <unordered_set>
#include <filesystem>

#include "src/utils.h"
#include "src/graph.h"
#include "lps.h"

using namespace std;

Logger main_logger("main_logger");
Logger parser_logger("parser_logger");
Logger lcp_logger("lcp_logger");
Logger graph_logger("graph_logger");


int main(int argc, char** argv) {
    LCP_INIT();

    ArgParser parser("Build LCP graph from normal and tumor FASTA directories.");

    parser.add("normal_dir", "--normal", "-n", "Directory containing normal FASTA files", "string", true);
    parser.add("tumor_dir", "--tumor", "-t", "Directory containing tumor FASTA files", "string", true);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("k", "--k", "-k", "K-mer size", "int", true);
    parser.add("output", "--output", "-o", "Output graph name (without extension)", "string", true);

    map<string, any> args = parser.parse(argc, argv);
    string normal_dir = any_cast<string>(args["normal_dir"]);
    string tumor_dir = any_cast<string>(args["tumor_dir"]);
    int depth = any_cast<int>(args["depth"]);
    int k = any_cast<int>(args["k"]);
    string output_name = any_cast<string>(args["output"]);

    main_logger.debug("Normal directory: " + normal_dir);
    main_logger.debug("Tumor directory: " + tumor_dir);
    main_logger.debug("Depth: " + to_string(depth));
    main_logger.debug("K: " + to_string(k));
    main_logger.debug("Output: " + output_name);
    
    main_logger.debug("Scanning directories for FASTA files...");
    vector<string> normal_files, tumor_files;
    for (const auto& entry : std::filesystem::directory_iterator(normal_dir)) {
        if (entry.is_regular_file()) {
            string ext = entry.path().extension().string();
            if (ext == ".fasta" || ext == ".fa" || ext == ".fna") {
                normal_files.push_back(entry.path().string());
            }
        }
    }
    for (const auto& entry : std::filesystem::directory_iterator(tumor_dir)) {
        if (entry.is_regular_file()) {
            string ext = entry.path().extension().string();
            if (ext == ".fasta" || ext == ".fa" || ext == ".fna") {
                tumor_files.push_back(entry.path().string());
            }
        }
    }
    
    main_logger.debug("Found " + to_string(normal_files.size()) + " normal files");
    main_logger.debug("Found " + to_string(tumor_files.size()) + " tumor files");
    
    graph_logger.debug("Initializing graph with k=" + to_string(k));
    LCPGraph graph(k);
    
    main_logger.debug("Processing normal samples (color 0)...");
    int normal_seq_count = 0;
    for (size_t i = 0; i < normal_files.size(); ++i) {
        const auto& file = normal_files[i];
        main_logger.debug("[" + to_string(i+1) + "/" + to_string(normal_files.size()) + "] Processing: " + file);
        FastaStreamer streamer(file, false);
        for (const FastaPair& pair : streamer) {
            struct lps lcp_string;
            init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
            lps_deepen(&lcp_string, depth);
            graph.addString(&lcp_string, 0);
            free_lps(&lcp_string);
            normal_seq_count++;
            break;
        }
    }
    main_logger.debug("Processed " + to_string(normal_seq_count) + " normal sequences");
    
    main_logger.debug("Processing tumor samples (colors 1-" + to_string(tumor_files.size()) + ")...");
    int tumor_color = 1;
    int tumor_seq_count = 0;
    for (size_t i = 0; i < tumor_files.size(); ++i) {
        const auto& file = tumor_files[i];
        main_logger.debug("[" + to_string(i+1) + "/" + to_string(tumor_files.size()) + "] Processing: " + file + " (color " + to_string(tumor_color) + ")");
        FastaStreamer streamer(file, false);
        for (const FastaPair& pair : streamer) {
            struct lps lcp_string;
            init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
            lps_deepen(&lcp_string, depth);
            graph.addString(&lcp_string, tumor_color);
            free_lps(&lcp_string);
            tumor_seq_count++;
            break;
        }
        tumor_color++;
    }
    main_logger.debug("Processed " + to_string(tumor_seq_count) + " tumor sequences");
    
    graph_logger.debug("Exporting graph to " + output_name + "...");
    graph.exportGFA(output_name);
    graph_logger.debug("Graph exported to " + output_name + ".gfa and " + output_name + ".csv");

    main_logger.debug("Done!");

    return 0;
}
