#ifndef GRAPH_H
#define GRAPH_H

#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core.h"
#include "lps.h"

typedef uint64_t node_hash;

struct Node {
    node_hash hash;
    uint32_t colors;
    unsigned int count;
    
    Node() : hash(0), colors(0), count(0) {}
    Node(node_hash h) : hash(h), colors(0), count(0) {}
    
    bool hasColor(int color) const {
        return (colors >> color) & 1;
    }
    
    void setColor(int color) {
        colors |= (1U << color);
    }
};

class LCPGraph {
public:
    std::unordered_map<node_hash, Node> nodes;
    
    std::unordered_map<node_hash, std::unordered_set<node_hash>> next_edges;
    std::unordered_map<node_hash, std::unordered_set<node_hash>> prev_edges;

    LCPGraph() = default;
    ~LCPGraph() = default;
    
    void addString(struct lps* lcp_string, int color, const char* sequence);
    void exportGFA(const std::string& output_dir, const std::string& experiment_name, const std::map<uint32_t, uint32_t>* color_map = nullptr) const;
    size_t nodeCount() const { return nodes.size(); }
    void compressGraph(std::unordered_set<node_hash>* deletion_nodes = nullptr);

private:
    bool nodeExists(node_hash h) const;
    uint32_t getNodeColor(node_hash h) const;
    const std::unordered_set<node_hash>& getNextNodes(node_hash h) const;
    const std::unordered_set<node_hash>& getPrevNodes(node_hash h) const;
    bool isChainStart(node_hash h, const std::unordered_set<node_hash>& claimed) const;
};

#endif // GRAPH_H
