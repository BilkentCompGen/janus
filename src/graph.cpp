#include "graph.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <filesystem>

static Logger graph_logger("graph");

void LCPGraph::addString(struct lps* lcp_string, int color) {
    if (lcp_string->size == 0) return;
    
    // For k=1: each core becomes a single node
    for (int i = 0; i < lcp_string->size; ++i) {
        ublock hash = hash_core(&lcp_string->cores[i]);
        
        // Get or create node
        Node& node = nodes[hash];
        if (node.hash == 0) node.hash = hash;
        node.count++;
        node.setColor(color);
        
        // Add edge to next core (if exists)
        if (i + 1 < lcp_string->size) {
            ublock next_hash = hash_core(&lcp_string->cores[i + 1]);
            next_edges[hash].insert(next_hash);
            prev_edges[next_hash].insert(hash);
        }
    }
}

void LCPGraph::exportGFA(const std::string& output_dir, const std::string& experiment_name) const {
    namespace fs = std::filesystem;
    fs::path base_dir = output_dir;
    fs::path gfa_path = base_dir / "gfa" / (experiment_name + ".gfa");
    fs::path csv_path = base_dir / "csv" / (experiment_name + ".csv");
    
    fs::create_directories(gfa_path.parent_path());
    fs::create_directories(csv_path.parent_path());
    
    std::ofstream gfa_out(gfa_path);
    if (!gfa_out) return;
    
    std::ofstream csv_out(csv_path);
    if (!csv_out) return;
    
    gfa_out << "H\tVN:Z:1.0\n";
    csv_out << "Name,Colour\n";
    
    // Map hash to sequential ID
    std::unordered_map<ublock, size_t> hash_to_id;
    size_t node_id = 1;
    for (const auto& [hash, _] : nodes) {
        hash_to_id[hash] = node_id++;
    }
    
    // Write nodes
    for (const auto& [hash, node] : nodes) {
        size_t id = hash_to_id[hash];
        
        // GFA: S line
        gfa_out << "S\t" << id << "\t*\tLN:i:1\tCN:i:" << node.count << "\tHASH:i:" << hash;
        // Add color info
        for (int i = 0; i < 32; ++i) {
            if (node.hasColor(i)) {
                gfa_out << "\tC" << i << ":i:1";
            }
        }
        gfa_out << "\n";
        
        // CSV: color for visualization
        uint32_t c = 0;
        if (node.hasColor(0)) c = 0x000000;  // Black for normal
        else {
            // Find first tumor color
            for (int i = 1; i < 32; ++i) {
                if (node.hasColor(i)) {
                    c = 0xFF0000 + (i * 0x010000);  // Simple color scheme
                    break;
                }
            }
        }
        char color_buf[9];  // "#" + 6 hex digits + null terminator
        snprintf(color_buf, sizeof(color_buf), "#%06X", c);
        csv_out << id << "," << color_buf << "\n";
    }
    
    // Write edges
    for (const auto& [hash, next_set] : next_edges) {
        size_t from_id = hash_to_id[hash];
        for (ublock next_hash : next_set) {
            if (hash_to_id.count(next_hash) && next_hash != hash) {
                size_t to_id = hash_to_id[next_hash];
                gfa_out << "L\t" << from_id << "\t+\t" << to_id << "\t+\t0M\n";
            }
        }
    }
}
