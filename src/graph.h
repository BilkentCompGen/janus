#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "lps.h"
#include "core.h"

// Use uint64_t for hash+label combination (avoid conflict with lcptools' ublock)
typedef uint64_t node_hash_t;

struct Node {
    node_hash_t hash;
    uint32_t colors;
    unsigned int count;
    
    Node() : hash(0), colors(0), count(0) {}
    Node(node_hash_t h) : hash(h), colors(0), count(0) {}
    
    bool hasColor(int color) const {
        return (colors >> color) & 1;
    }
    
    void setColor(int color) {
        colors |= (1U << color);
    }
};

class LCPGraph {
private:
    std::unordered_map<node_hash_t, Node> nodes;
    
    std::unordered_map<node_hash_t, std::unordered_set<node_hash_t>> next_edges;
    std::unordered_map<node_hash_t, std::unordered_set<node_hash_t>> prev_edges;
    
    // Store positions for each node: hash -> vector of (start, end) pairs
    std::unordered_map<node_hash_t, std::vector<std::pair<uint64_t, uint64_t>>> node_positions;

public:
    LCPGraph() = default;
    ~LCPGraph() = default;
    
    void addString(struct lps* lcp_string, int color, const char* sequence = nullptr);
    
    void exportGFA(const std::string& output_dir, const std::string& experiment_name) const;
    
    const std::unordered_map<node_hash_t, Node>& getNodes() const { return nodes; }
    size_t nodeCount() const { return nodes.size(); }
    
    const std::unordered_set<node_hash_t>& getNext(node_hash_t hash) const {
        static const std::unordered_set<node_hash_t> empty;
        auto it = next_edges.find(hash);
        return (it != next_edges.end()) ? it->second : empty;
    }
    
    const std::unordered_set<node_hash_t>& getPrev(node_hash_t hash) const {
        static const std::unordered_set<node_hash_t> empty;
        auto it = prev_edges.find(hash);
        return (it != prev_edges.end()) ? it->second : empty;
    }
    
    const std::vector<std::pair<uint64_t, uint64_t>>& getPositions(node_hash_t hash) const {
        static const std::vector<std::pair<uint64_t, uint64_t>> empty;
        auto it = node_positions.find(hash);
        return (it != node_positions.end()) ? it->second : empty;
    }
    
    void compressGraph();
};

#endif // GRAPH_H
