#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "lps.h"
#include "core.h"

constexpr std::size_t MAX_COLORS = 32;
typedef unsigned int ublock;

struct Node {
    ublock core_hash;
    std::unordered_set<ublock> next;
    std::unordered_set<ublock> prev;
    unsigned int colors[MAX_COLORS] = {0};
    unsigned int count = 0;
    
    Node(ublock core_hash) : core_hash(core_hash) {}
    
    bool isCompressible() const {
        return prev.size() == 1 && next.size() == 1;
    }
};

struct CompressedNode {
    std::vector<ublock> path;
    std::unordered_set<ublock> next;
    unsigned int count = 0;
    unsigned int colors[MAX_COLORS] = {0};
    
    size_t length() const { return path.size(); }
};

class LCPGraph {
private:
    std::unordered_map<ublock, Node> nodes;
    int k;
    
    std::unordered_map<ublock, CompressedNode> compressGraph() const;
    void printStatistics(const std::string& label) const;
    void printCompressedStatistics(const std::unordered_map<ublock, CompressedNode>& compressed, const std::string& label) const;
    
public:
    LCPGraph(int k = 1) : k(k) {}
    ~LCPGraph() = default;

    void addString(struct lps* lcp_string, int color);
    void exportGFA(const std::string& graph_name) const;
};

#endif // GRAPH_H
