#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <filesystem>

#include "utils.h"
#include "graph.h"
#include "lps.h"

static Logger compression_logger("compression");
static Logger statistics_logger("statistics");


void LCPGraph::addString(struct lps* lcp_string, int color) {
    if (is_compressed) {
        throw std::runtime_error("Cannot add strings to a compressed graph");
    }
    
    if (lcp_string->size == 0) return;
    
    ublock curr = hash_k_cores(&lcp_string->cores[0], k);
    
    for (int i = 0; i <= static_cast<int>(lcp_string->size) - k; ++i) {
        if (i > 0) {
            curr = update_k_cores_hash(curr, &lcp_string->cores[i - 1], &lcp_string->cores[i + k - 1], k);
        }
        
        // Update node data
        NodeData& data = node_dict.try_emplace(curr, curr).first->second;
        data.count++;
        data.colors[color]++;
        
        // Update node edges
        NodeEdges& edges = nodes[curr];
        
        if (i + 1 <= static_cast<int>(lcp_string->size) - k) {
            ublock next_hash = update_k_cores_hash(curr, &lcp_string->cores[i], &lcp_string->cores[i + k], k);
            edges.next.emplace(next_hash);
            
            // Update next node
            NodeData& next_data = node_dict.try_emplace(next_hash, next_hash).first->second;
            next_data.count++;
            next_data.colors[color]++;
            
            NodeEdges& next_edges = nodes[next_hash];
            next_edges.prev.emplace(curr);
        }
    }
}

CompressedNode LCPGraph::createCompressedNode(const std::vector<ublock>& path) const {
    CompressedNode comp_node;
    comp_node.path = path;
    
    if (!path.empty()) {
        // Get outgoing edges from the last node in path
        auto edges_it = nodes.find(path.back());
        if (edges_it != nodes.end()) {
            comp_node.next = edges_it->second.next;
        }
        
        // Aggregate counts and colors from all nodes in path
        for (ublock h : path) {
            auto data_it = node_dict.find(h);
            if (data_it != node_dict.end()) {
                const NodeData& data = data_it->second;
                comp_node.count += data.count;
                for (size_t i = 0; i < MAX_COLORS; ++i) {
                    comp_node.colors[i] += data.colors[i];
                }
            }
        }
    }
    return comp_node;
}

void LCPGraph::compress() {
    if (is_compressed) {
        return; // Already compressed
    }
    
    compression_logger.debug("Starting graph compression...");
    
    std::unordered_set<ublock> processed;
    
    for (const auto& [hash, edges] : nodes) {
        if (processed.count(hash)) continue;
        
        if (!edges.isCompressible()) {
            compressed_nodes[hash] = createCompressedNode({hash});
            processed.insert(hash);
            continue;
        }
        
        // Build chain backward
        ublock chain_start = hash;
        ublock current = hash;
        ublock prev_node = 0;
        
        while (true) {
            auto curr_edges_it = nodes.find(current);
            if (curr_edges_it == nodes.end() || curr_edges_it->second.prev.size() != 1) break;
            
            ublock pred = *curr_edges_it->second.prev.begin();
            auto pred_edges_it = nodes.find(pred);
            if (pred == prev_node || processed.count(pred) || 
                pred_edges_it == nodes.end() || !pred_edges_it->second.isCompressible()) break;
            
            prev_node = current;
            current = pred;
            chain_start = current;
        }
        
        // Build chain forward
        std::vector<ublock> chain;
        current = chain_start;
        
        while (true) {
            if (processed.count(current)) break;
            
            chain.push_back(current);
            processed.insert(current);
            
            auto curr_edges_it = nodes.find(current);
            if (curr_edges_it == nodes.end() || !curr_edges_it->second.isCompressible() || 
                curr_edges_it->second.next.size() != 1) break;
            
            ublock next = *curr_edges_it->second.next.begin();
            auto next_edges_it = nodes.find(next);
            if (next_edges_it == nodes.end() || processed.count(next) || 
                !next_edges_it->second.isCompressible()) break;
            
            current = next;
        }
        
        compressed_nodes[chain_start] = createCompressedNode(chain);
    }
    
    // Update edges to point to compressed nodes
    std::unordered_map<ublock, ublock> node_to_compressed;
    for (const auto& [comp_hash, comp_node] : compressed_nodes) {
        for (ublock h : comp_node.path) {
            node_to_compressed[h] = comp_hash;
        }
    }
    
    for (auto& [comp_hash, comp_node] : compressed_nodes) {
        std::unordered_set<ublock> updated_next;
        for (ublock next_hash : comp_node.next) {
            if (node_to_compressed.count(next_hash)) {
                updated_next.insert(node_to_compressed[next_hash]);
            }
        }
        comp_node.next = updated_next;
    }
    
    // Clear the node edges to free memory
    nodes.clear();
    is_compressed = true;
    
    compression_logger.debug("Compression complete: " + std::to_string(node_dict.size()) + 
                            " -> " + std::to_string(compressed_nodes.size()) + " nodes");
}

std::vector<ublock> LCPGraph::getPath(ublock compressed_node_id) const {
    auto it = compressed_nodes.find(compressed_node_id);
    if (it != compressed_nodes.end()) {
        return it->second.path;
    }
    return {};
}

static const uint32_t COLOR_PALETTE[] = {
    0x000000, 0xFF0000, 0x00FF00, 0x0000FF, 0xFFFF00, 0xFF00FF, 0x00FFFF, 0xFF8000,
    0x8000FF, 0xFF0080, 0x80FF00, 0x00FF80, 0x0080FF, 0xFF4040, 0x40FF40, 0x4040FF,
    0xFFFF80, 0xFF80FF, 0x80FFFF, 0xFF8040, 0x8040FF, 0xFF4080, 0x80FF40, 0x40FF80,
    0x4080FF, 0xC00000, 0x00C000, 0x0000C0, 0xC0C000, 0xC000C0, 0x00C0C0, 0xC06000
};

static std::string getNodeColor(const CompressedNode& node) {
    if (node.colors[0] > 0) return "#000000";
    
    for (size_t i = 1; i < MAX_COLORS && i < 32; ++i) {
        if (node.colors[i] > 0) {
            uint32_t c = COLOR_PALETTE[i];
            char buf[8];
            snprintf(buf, sizeof(buf), "#%02X%02X%02X", 
                (c >> 16) & 0xFF, (c >> 8) & 0xFF, c & 0xFF);
            return buf;
        }
    }
    
    return "#969696";
}

void LCPGraph::exportGFA(const std::string& output_dir, const std::string& experiment_name) const {
    if (!is_compressed) {
        throw std::runtime_error("Graph must be compressed before exporting");
    }
    
    namespace fs = std::filesystem;
    fs::path base_dir = output_dir;
    fs::path gfa_path = base_dir / "gfa" / (experiment_name + ".gfa");
    fs::path csv_path = base_dir / "csv" / (experiment_name + ".csv");
    
    fs::create_directories(gfa_path.parent_path());
    fs::create_directories(csv_path.parent_path());
    
    std::ofstream out(gfa_path);
    if (!out) return;
    out << "H\tVN:Z:1.0\n";
    
    std::unordered_map<uint64_t, size_t> hash_to_id;
    size_t node_id = 1;
    for (const auto& [hash, _] : compressed_nodes) hash_to_id[hash] = node_id++;
    
    for (const auto& [hash, node] : compressed_nodes) {
        out << "S\t" << hash_to_id[hash] << "\t*\tLN:i:" << node.length() << "\tCN:i:" << node.count << "\tHASH:i:" << hash;
        for (size_t i = 0; i < MAX_COLORS; ++i) {
            if (node.colors[i] > 0) out << "\tC" << i << ":i:" << node.colors[i];
        }
        out << "\n";
    }
    
    for (const auto& [hash, node] : compressed_nodes) {
        size_t from_id = hash_to_id[hash];
        for (uint64_t next_hash : node.next) {
            if (hash_to_id.count(next_hash) && next_hash != hash) {
                out << "L\t" << from_id << "\t+\t" << hash_to_id[next_hash] << "\t+\t0M\n";
            }
        }
    }
    
    std::ofstream col(csv_path);
    if (!col) return;
    col << "Name,Colour\n";
    for (const auto& [hash, node] : compressed_nodes) {
        col << hash_to_id.at(hash) << "," << getNodeColor(node) << "\n";
    }
}
