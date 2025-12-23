#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

#include "lps.h"
#include "core.h"

typedef uint64_t node_hash;

struct PositionInfo {
    int cell_id;
    std::string chrom;
    uint64_t start;
    uint64_t end;
    
    PositionInfo() : cell_id(0), chrom(""), start(0), end(0) {}
    PositionInfo(int cid, const std::string& ch, uint64_t s, uint64_t e) : cell_id(cid), chrom(ch), start(s), end(e) {}
};

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

    std::unordered_map<node_hash, std::vector<PositionInfo>> node_positions;
    std::unordered_map<int, std::unordered_map<std::string, std::string>> sequences;

    LCPGraph() = default;
    ~LCPGraph() = default;
    
    void addString(struct lps* lcp_string, int color, const char* sequence, const std::string& chrom);
    void exportGFA(const std::string& output_dir, const std::string& experiment_name) const;
    size_t nodeCount() const { return nodes.size(); }
    std::vector<PositionInfo> getPositions(node_hash hash, int cell_id = -1) const;
    std::string getSequence(int cell_id, const std::string& chrom) const;
    void setSequence(int cell_id, const std::string& chrom, const std::string& sequence);
    void compressGraph(std::unordered_set<node_hash>* deletion_nodes = nullptr);

private:
    bool nodeExists(node_hash h) const;
    uint32_t getNodeColor(node_hash h) const;
    const std::unordered_set<node_hash>& getNextNodes(node_hash h) const;
    const std::unordered_set<node_hash>& getPrevNodes(node_hash h) const;
    bool isChainStart(node_hash h, const std::unordered_set<node_hash>& claimed) const;
};

#endif // GRAPH_H
