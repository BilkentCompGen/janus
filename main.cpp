#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <any>
#include <unordered_set>
#include <filesystem>

#include "src/utils.h"
#include "src/graph.h"
#include "src/sv.h"
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
    
    main_logger.debug("Compressing graph...");
    graph.compress();
    main_logger.debug("Graph compressed: " + to_string(graph.originalNodeCount()) + " -> " + to_string(graph.nodeCount()) + " nodes");
    
    graph_logger.debug("Exporting graph to " + output_name + "...");
    graph.exportGFA(output_name);
    graph_logger.debug("Graph exported to " + output_name + ".gfa and " + output_name + ".csv");


    main_logger.debug("Detecting structural variants...");
    auto sv_calls = detectPresenceAbsenceSVs(graph, reference, static_cast<int>(tumor_files.size()));
    
    int deletion_count = 0;
    int insertion_count = 0;
    int duplication_count = 0;
    for (const auto& sv : sv_calls) {
        if (sv.type == SVType::DELETION) {
            deletion_count++;
        } else if (sv.type == SVType::INSERTION) {
            insertion_count++;
        } else if (sv.type == SVType::DUPLICATION) {
            duplication_count++;
        }
    }
    
    main_logger.debug("Found " + to_string(sv_calls.size()) + " structural variants:");
    main_logger.debug("  Deletions: " + to_string(deletion_count));
    main_logger.debug("  Insertions: " + to_string(insertion_count));
    main_logger.debug("  Duplications: " + to_string(duplication_count));
    
    int print_count = min(10, static_cast<int>(sv_calls.size()));
    if (print_count > 0) {
        main_logger.debug("First " + to_string(print_count) + " SVs:");
        
        for (int i = 0; i < print_count; ++i) {
            const auto& sv = sv_calls[i];
            string type_str;
            switch (sv.type) {
                case SVType::DELETION: type_str = "DELETION"; break;
                case SVType::INSERTION: type_str = "INSERTION"; break;
                case SVType::DUPLICATION: type_str = "DUPLICATION"; break;
            }
            main_logger.debug("  " + type_str + " at node_hash: " + to_string(sv.node_hash) + 
                             " (normal:" + to_string(sv.normal_count) + " tumor:" + to_string(sv.tumor_count) + ")");
            
            // Print path info
            string path_string = "Path (" + to_string(sv.path.size()) + " hashes): ";
            size_t max_hashes = min(size_t(5), sv.path.size());
            for (size_t j = 0; j < max_hashes; ++j) {
                path_string += to_string(sv.path[j]) + " ";
            }
            if (sv.path.size() > max_hashes) {
                path_string += "...";
            }
            main_logger.debug("    " + path_string);
            
            // Print flanking info for insertions
            if (sv.type == SVType::INSERTION) {
                if (!sv.flanking_before.empty()) {
                    string flank_str = "Flanking before (" + to_string(sv.flanking_before.size()) + " hashes): ";
                    size_t max_f = min(size_t(3), sv.flanking_before.size());
                    for (size_t j = 0; j < max_f; ++j) {
                        flank_str += to_string(sv.flanking_before[j]) + " ";
                    }
                    if (sv.flanking_before.size() > max_f) flank_str += "...";
                    main_logger.debug("    " + flank_str);
                }
                if (!sv.flanking_after.empty()) {
                    string flank_str = "Flanking after (" + to_string(sv.flanking_after.size()) + " hashes): ";
                    size_t max_f = min(size_t(3), sv.flanking_after.size());
                    for (size_t j = 0; j < max_f; ++j) {
                        flank_str += to_string(sv.flanking_after[j]) + " ";
                    }
                    if (sv.flanking_after.size() > max_f) flank_str += "...";
                    main_logger.debug("    " + flank_str);
                }
            }
            
            // Print genomic positions
            if (!sv.genomic_starts.empty()) {
                string pos_string = "Genomic positions: ";
                size_t max_pos = min(size_t(3), sv.genomic_starts.size());
                for (size_t j = 0; j < max_pos; ++j) {
                    pos_string += to_string(sv.genomic_starts[j]) + "-" + to_string(sv.genomic_ends[j]) + " ";
                }
                if (sv.genomic_starts.size() > max_pos) {
                    pos_string += "... (" + to_string(sv.genomic_starts.size()) + " total)";
                }
                main_logger.debug("    " + pos_string);
            } else {
                main_logger.debug("    No genomic positions found");
            }
        }
    }

    // Save SVs to TSV file
    string sv_filename = output_name + "_svs.tsv";
    ofstream sv_file(sv_filename);
    if (sv_file) {
        sv_file << "type\tchrom\tstart\tend\ttumor_samples\tnormal_count\ttumor_count\n";
        for (const auto& sv : sv_calls) {
            string type_str;
            switch (sv.type) {
                case SVType::DELETION: type_str = "DEL"; break;
                case SVType::INSERTION: type_str = "INS"; break;
                case SVType::DUPLICATION: type_str = "DUP"; break;
            }
            
            // Build tumor samples string
            string samples_str;
            for (size_t i = 0; i < sv.tumor_samples.size(); ++i) {
                if (i > 0) samples_str += ",";
                samples_str += to_string(sv.tumor_samples[i]);
            }
            if (samples_str.empty()) samples_str = ".";
            
            for (size_t i = 0; i < sv.genomic_starts.size(); ++i) {
                sv_file << type_str << "\t" 
                        << (sv.chrom.empty() ? "." : sv.chrom) << "\t"
                        << sv.genomic_starts[i] << "\t" 
                        << sv.genomic_ends[i] << "\t"
                        << samples_str << "\t"
                        << sv.normal_count << "\t"
                        << sv.tumor_count << "\n";
            }
        }
        sv_file.close();
        main_logger.debug("SVs saved to " + sv_filename);
    }

    main_logger.debug("Done!");

    return 0;
}
