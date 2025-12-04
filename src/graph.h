#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

#include "lps.h"
#include "core.h"

constexpr std::size_t MAX_COLORS = 32;
typedef unsigned int ublock;

// Node data: core hash, colors, and count
struct NodeData {
    ublock core_hash;
    unsigned int colors[MAX_COLORS] = {0};
    unsigned int count = 0;
    
    NodeData() : core_hash(0) {}
    NodeData(ublock hash) : core_hash(hash) {}
};

// Node edges: next and prev connections (cleared after compression)
struct NodeEdges {
    std::unordered_set<ublock> next;
    std::unordered_set<ublock> prev;
    
    bool isCompressible() const {
        return prev.size() == 1 && next.size() == 1;
    }
};

// Compressed node: path of core hashes and outgoing edges
struct CompressedNode {
    std::vector<ublock> path;
    std::unordered_set<ublock> next;
    unsigned int count = 0;
    unsigned int colors[MAX_COLORS] = {0};
    
    size_t length() const { return path.size(); }
};

// Unified LCP Graph - handles both uncompressed and compressed states
class LCPGraph {
private:
    bool is_compressed = false;
    int k;
    
    // Node data (kept after compression for color/count info)
    std::unordered_map<ublock, NodeData> node_dict;
    
    // Node edges (cleared after compression)
    std::unordered_map<ublock, NodeEdges> nodes;
    
    // Compressed nodes (populated after compression)
    std::unordered_map<ublock, CompressedNode> compressed_nodes;
    
    // Helper for compression
    CompressedNode createCompressedNode(const std::vector<ublock>& path) const;

public:
    LCPGraph(int k = 1) : k(k) {}
    ~LCPGraph() = default;
    
    // Add string to graph (throws if already compressed)
    void addString(struct lps* lcp_string, int color);
    
    // Compress the graph (frees node edges, sets is_compressed)
    void compress();
    
    // Export to GFA format
    void exportGFA(const std::string& graph_name) const;
    
    // Get path for a compressed node
    std::vector<ublock> getPath(ublock compressed_node_id) const;
    
    // Accessors
    bool isCompressed() const { return is_compressed; }
    int getK() const { return k; }
    
    const std::unordered_map<ublock, NodeData>& getNodeDict() const { return node_dict; }
    const std::unordered_map<ublock, CompressedNode>& getCompressedNodes() const { return compressed_nodes; }
    
    // Get node count
    size_t nodeCount() const { 
        return is_compressed ? compressed_nodes.size() : node_dict.size(); 
    }
    size_t originalNodeCount() const { return node_dict.size(); }
};

#endif // GRAPH_H
