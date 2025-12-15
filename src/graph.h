#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "lps.h"
#include "core.h"

typedef unsigned int ublock;

// Simplified node: just hash, color bitmask, and global count
struct Node {
    ublock hash;
    uint32_t colors;  // Bitmask: bit i set = color i present
    unsigned int count;  // Global count across all colors
    
    Node() : hash(0), colors(0), count(0) {}
    Node(ublock h) : hash(h), colors(0), count(0) {}
    
    // Check if color is present
    bool hasColor(int color) const {
        return (colors >> color) & 1;
    }
    
    // Set color bit
    void setColor(int color) {
        colors |= (1U << color);
    }
};

// Simple graph - no compression
class LCPGraph {
private:
    // Node data: hash -> Node
    std::unordered_map<ublock, Node> nodes;
    
    // Edges: hash -> {next hashes}
    std::unordered_map<ublock, std::unordered_set<ublock>> next_edges;
    std::unordered_map<ublock, std::unordered_set<ublock>> prev_edges;

public:
    LCPGraph() = default;
    ~LCPGraph() = default;
    
    // Add string to graph (k=1: each core is a node)
    void addString(struct lps* lcp_string, int color);
    
    // Export to GFA format
    void exportGFA(const std::string& output_dir, const std::string& experiment_name) const;
    
    // Accessors
    const std::unordered_map<ublock, Node>& getNodes() const { return nodes; }
    size_t nodeCount() const { return nodes.size(); }
    
    // Get next/prev edges for a node
    const std::unordered_set<ublock>& getNext(ublock hash) const {
        static const std::unordered_set<ublock> empty;
        auto it = next_edges.find(hash);
        return (it != next_edges.end()) ? it->second : empty;
    }
};

#endif // GRAPH_H
