#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <any>
#include <unordered_set>
#include <filesystem>
#include <regex>
#include <algorithm>

#include "src/utils.h"
#include "src/graph.h"
#include "src/sv.h"
#include "lps.h"

using namespace std;

Logger main_logger("main_logger");
Logger parser_logger("parser_logger");
Logger lcp_logger("lcp_logger");
Logger graph_logger("graph_logger");

// Read tumor.tipnode.map and return map: nodeX -> carrier_id
map<string, int> readTipNodeMap(const string& map_file) {
    map<string, int> node_to_carrier;
    ifstream fin(map_file);
    if (!fin) {
        main_logger.debug("Warning: Could not open map file: " + map_file);
        return node_to_carrier;
    }
    string line;
    while (getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        // Format: tip_node<TAB>cell_count<TAB>cells
        // e.g., node11	1	3
        size_t tab1 = line.find('\t');
        size_t tab2 = line.find('\t', tab1 + 1);
        if (tab1 != string::npos && tab2 != string::npos) {
            string node_name = line.substr(0, tab1);
            string carrier_str = line.substr(tab2 + 1);
            int carrier_id = stoi(carrier_str);
            node_to_carrier[node_name] = carrier_id;
        }
    }
    return node_to_carrier;
}

// Extract nodeX from filename like "node11.parental_0.fa" -> "node11"
string extractNodeName(const string& filepath) {
    string filename = std::filesystem::path(filepath).filename().string();
    // Match nodeX pattern
    regex node_regex("(node\\d+)");
    smatch match;
    if (regex_search(filename, match, node_regex)) {
        return match[1].str();
    }
    return "";
}

int main(int argc, char** argv) {
    LCP_INIT();

    ArgParser parser("Build LCP graph from normal and tumor FASTA directories.");

    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("k", "--k", "-k", "K-mer size", "int", true);
    parser.add("output", "--output", "-o", "Output graph name (without extension)", "string", true);

    map<string, any> args = parser.parse(argc, argv);
    int depth = any_cast<int>(args["depth"]);
    int k = any_cast<int>(args["k"]);
    string output_name = any_cast<string>(args["output"]);
    
    // Fixed paths
    string normal_dir = "simulation/normal";
    string tumor_dir = "simulation/tumors";
    string map_file = "simulation/tumor.tipnode.map";
    
    std::filesystem::path base_output_dir = "output";
    std::filesystem::path tsv_dir = base_output_dir / "tsv"; 

    std::filesystem::create_directories(tsv_dir);
    main_logger.debug("Output dir: " + base_output_dir.string());
    main_logger.debug("Output: " + output_name);

    main_logger.debug("Normal directory: " + normal_dir);
    main_logger.debug("Tumor directory: " + tumor_dir);
    main_logger.debug("Depth: " + to_string(depth));
    main_logger.debug("K: " + to_string(k));
    
    // Read node -> carrier mapping
    main_logger.debug("Reading tip node map from: " + map_file);
    map<string, int> node_to_carrier = readTipNodeMap(map_file);
    main_logger.debug("Loaded " + to_string(node_to_carrier.size()) + " node->carrier mappings");
    
    main_logger.debug("Scanning directories for FASTA files (parental_0 only)...");
    vector<string> normal_files, tumor_files;
    
    // Only pick parental_0 files for normal
    for (const auto& entry : std::filesystem::directory_iterator(normal_dir)) {
        if (entry.is_regular_file()) {
            string filename = entry.path().filename().string();
            string ext = entry.path().extension().string();
            if ((ext == ".fasta" || ext == ".fa" || ext == ".fna") && 
                filename.find("parental_0") != string::npos) {
                normal_files.push_back(entry.path().string());
            }
        }
    }
    
    // Only pick parental_0 files for tumor
    for (const auto& entry : std::filesystem::directory_iterator(tumor_dir)) {
        if (entry.is_regular_file()) {
            string filename = entry.path().filename().string();
            string ext = entry.path().extension().string();
            if ((ext == ".fasta" || ext == ".fa" || ext == ".fna") && 
                filename.find("parental_0") != string::npos) {
                tumor_files.push_back(entry.path().string());
            }
        }
    }
    
    // Sort for deterministic order
    sort(normal_files.begin(), normal_files.end());
    sort(tumor_files.begin(), tumor_files.end());
    
    main_logger.debug("Found " + to_string(normal_files.size()) + " normal parental_0 files");
    main_logger.debug("Found " + to_string(tumor_files.size()) + " tumor parental_0 files");
    
    graph_logger.debug("Initializing graph with k=" + to_string(k));
    LCPGraph graph(k);
    LinearReference reference;
    
    main_logger.debug("Processing normal samples (color 0)...");
    int normal_seq_count = 0;
    for (size_t i = 0; i < normal_files.size(); ++i) {
        const auto& file = normal_files[i];
        main_logger.debug("[" + to_string(i+1) + "/" + to_string(normal_files.size()) + "] Processing: " + file);
        FastaStreamer streamer(file, false);
        for (const FastaPair& pair : streamer) {
            // Extract chromosome from FASTA header (e.g., ">1" -> "1", ">chr1" -> "chr1")
            string chrom = pair.first;
            if (!chrom.empty() && chrom[0] == '>') {
                chrom = chrom.substr(1);
            }
            // Remove anything after first space
            size_t space_pos = chrom.find(' ');
            if (space_pos != string::npos) {
                chrom = chrom.substr(0, space_pos);
            }
            
            struct lps lcp_string;
            init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
            lps_deepen(&lcp_string, depth);
            graph.addString(&lcp_string, 0);
            
            // Build linear reference with real genomic positions
            for (int j = 0; j <= lcp_string.size - k; ++j) {
                ublock hash = hash_k_cores(&lcp_string.cores[j], k);
                uint64_t start = lcp_string.cores[j].start;
                uint64_t end = lcp_string.cores[j + k - 1].end;
                reference.addCore(hash, chrom, start, end);
            }
            
            free_lps(&lcp_string);
            normal_seq_count++;
            break;  // Only first sequence per file
        }
    }
    main_logger.debug("Processed " + to_string(normal_seq_count) + " normal sequences");
    
    // Build index for fast hash lookup
    reference.buildIndex();
    main_logger.debug("Linear reference built with " + to_string(reference.size()) + " positions");
    
    main_logger.debug("Processing tumor samples (using carrier IDs from map)...");
    int tumor_seq_count = 0;
    int max_carrier_id = 0;
    for (size_t i = 0; i < tumor_files.size(); ++i) {
        const auto& file = tumor_files[i];
        string node_name = extractNodeName(file);
        
        int carrier_id = 0;
        if (node_to_carrier.count(node_name)) {
            carrier_id = node_to_carrier[node_name];
            if (carrier_id > max_carrier_id) max_carrier_id = carrier_id;
        } else {
            main_logger.debug("Warning: No carrier mapping for " + node_name + ", skipping");
            continue;
        }
        
        main_logger.debug("[" + to_string(i+1) + "/" + to_string(tumor_files.size()) + "] Processing: " + file + " (" + node_name + " -> carrier " + to_string(carrier_id) + ")");
        FastaStreamer streamer(file, false);
        for (const FastaPair& pair : streamer) {
            struct lps lcp_string;
            init_lps(&lcp_string, pair.second.c_str(), pair.second.size());
            lps_deepen(&lcp_string, depth);
            graph.addString(&lcp_string, carrier_id);
            free_lps(&lcp_string);
            tumor_seq_count++;
            break;
        }
    }
    main_logger.debug("Processed " + to_string(tumor_seq_count) + " tumor sequences (max carrier ID: " + to_string(max_carrier_id) + ")");
    
    main_logger.debug("Compressing graph...");
    graph.compress();
    main_logger.debug("Graph compressed: " + to_string(graph.originalNodeCount()) + " -> " + to_string(graph.nodeCount()) + " nodes");
    
    graph_logger.debug("Exporting graph to " + base_output_dir.string() + "...");
    graph.exportGFA(base_output_dir.string(), output_name);
    graph_logger.debug("Graph exported.");


    main_logger.debug("Detecting structural variants...");
    auto sv_calls = detectDeletions(graph, reference, max_carrier_id);
    
    int deletion_count = 0;
    for (const auto& sv : sv_calls) {
        if (sv.type == SVType::DELETION) deletion_count++;
    }
    
    main_logger.debug("Found " + to_string(sv_calls.size()) + " structural variants:");
    main_logger.debug("  Deletions: " + to_string(deletion_count));

    // Save SVs to TSV file
    std::filesystem::path sv_filename = tsv_dir / (output_name + ".tsv");
    exportTSV(sv_calls, sv_filename.string());
    main_logger.debug("SVs saved to " + sv_filename.string());


    main_logger.debug("Done!");

    return 0;
}
