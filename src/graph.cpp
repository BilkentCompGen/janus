#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "utils.h"
#include "graph.h"
#include "lps.h"

static Logger compression_logger("compression");
static Logger statistics_logger("statistics");


void LCPGraph::addString(struct lps* lcp_string, int color) {
    if (lcp_string->size == 0) return;
    
    ublock curr = hash_k_cores(&lcp_string->cores[0], k);
    
    for (int i = 0; i <= static_cast<int>(lcp_string->size) - k; ++i) {
        if (i > 0) {
            curr = update_k_cores_hash(curr, &lcp_string->cores[i - 1], &lcp_string->cores[i + k - 1], k);
        }
        
        Node& node = nodes.try_emplace(curr, curr).first->second;
        node.count++;
        node.colors[color]++;
        
        if (i + 1 <= static_cast<int>(lcp_string->size) - k) {
            ublock next_hash = update_k_cores_hash(curr, &lcp_string->cores[i], &lcp_string->cores[i + k], k);
            node.next.emplace(next_hash);
            
            Node& next_node = nodes.try_emplace(next_hash, next_hash).first->second;
            next_node.prev.emplace(curr);
            next_node.colors[color]++;
            next_node.count++;
        }
    }
}

static CompressedNode createCompressedNode(const std::vector<ublock>& path, const std::unordered_map<ublock, Node>& nodes) {
    CompressedNode comp_node;
    comp_node.path = path;
    
    if (!path.empty()) {
        const Node& last_node = nodes.at(path.back());
        comp_node.next = last_node.next;
        
        for (ublock h : path) {
            const Node& node = nodes.at(h);
            comp_node.count += node.count;
            for (size_t i = 0; i < MAX_COLORS; ++i) {
                comp_node.colors[i] += node.colors[i];
            }
        }
    }
    return comp_node;
}

CompressedGraph LCPGraph::compress() const {
    compression_logger.debug("Starting graph compression...");
    
    std::unordered_map<ublock, CompressedNode> compressed;
    std::unordered_set<ublock> processed;
    
    for (const auto& [hash, node] : nodes) {
        if (processed.count(hash)) continue;
        
        if (!node.isCompressible()) {
            compressed[hash] = createCompressedNode({hash}, nodes);
            processed.insert(hash);
            continue;
        }
        
        // Build chain backward
        ublock chain_start = hash;
        ublock current = hash;
        ublock prev_node = 0;
        
        while (true) {
            const Node& curr_node = nodes.at(current);
            if (curr_node.prev.size() != 1) break;
            
            ublock pred = *curr_node.prev.begin();
            if (pred == prev_node || processed.count(pred) || !nodes.at(pred).isCompressible()) break;
            
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
            
            const Node& curr_node = nodes.at(current);
            if (!curr_node.isCompressible() || curr_node.next.size() != 1) break;
            
            ublock next = *curr_node.next.begin();
            if (!nodes.count(next) || processed.count(next) || !nodes.at(next).isCompressible()) break;
            
            current = next;
        }
        
        compressed[chain_start] = createCompressedNode(chain, nodes);
    }
    
    // Update edges to point to compressed nodes
    std::unordered_map<ublock, ublock> node_to_compressed;
    for (const auto& [comp_hash, comp_node] : compressed) {
        for (ublock h : comp_node.path) {
            node_to_compressed[h] = comp_hash;
        }
    }
    
    for (auto& [comp_hash, comp_node] : compressed) {
        std::unordered_set<ublock> updated_next;
        for (ublock next_hash : comp_node.next) {
            if (node_to_compressed.count(next_hash)) {
                updated_next.insert(node_to_compressed[next_hash]);
            }
        }
        comp_node.next = updated_next;
    }
    
    compression_logger.debug("Compression complete: " + std::to_string(nodes.size()) + " -> " + std::to_string(compressed.size()) + " nodes");
    return CompressedGraph(std::move(compressed));
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

void CompressedGraph::exportGFA(const std::string& graph_name) const {
    std::ofstream out(graph_name + ".gfa");
    if (!out) return;
    out << "H\tVN:Z:1.0\n";
    
    std::unordered_map<ublock, size_t> hash_to_id;
    size_t node_id = 1;
    for (const auto& [hash, _] : nodes) hash_to_id[hash] = node_id++;
    
    for (const auto& [hash, node] : nodes) {
        out << "S\t" << hash_to_id[hash] << "\t*\tLN:i:" << node.length() << "\tCN:i:" << node.count << "\tHASH:i:" << hash;
        for (size_t i = 0; i < MAX_COLORS; ++i) {
            if (node.colors[i] > 0) out << "\tC" << i << ":i:" << node.colors[i];
        }
        out << "\n";
    }
    
    for (const auto& [hash, node] : nodes) {
        size_t from_id = hash_to_id[hash];
        for (ublock next_hash : node.next) {
            if (hash_to_id.count(next_hash) && next_hash != hash) {
                out << "L\t" << from_id << "\t+\t" << hash_to_id[next_hash] << "\t+\t0M\n";
            }
        }
    }
    
    std::ofstream col(graph_name + ".csv");
    if (!col) return;
    col << "Name,Colour\n";
    for (const auto& [hash, node] : nodes) {
        col << hash_to_id.at(hash) << "," << getNodeColor(node) << "\n";
    }
}