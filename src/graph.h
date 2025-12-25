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

struct Occurrence {
    uint32_t seq_id;
    uint64_t start;
    uint64_t end;
};

struct BubbleResult {
    std::unordered_set<node_hash> bubble_starts;
    std::unordered_set<node_hash> bubble_ends;
};

struct BubblePath {
    std::vector<node_hash> nodes;  // Path'teki node hash'leri (start ve end dahil degil)
    uint32_t seq_id = 0;           // Bu path hangi sequence'ta
    uint64_t start_coord = 0;      // Path baslangic koordinati
    uint64_t end_coord = 0;        // Path bitis koordinati
};

struct Bubble {
    node_hash start_node;
    node_hash end_node;
    
    // Reference koordinatlari (start/end node'larindan)
    uint64_t ref_start;
    uint64_t ref_end;
    
    // Tumor koordinatlari
    uint64_t tumor_start;
    uint64_t tumor_end;
    
    // Path'ler
    BubblePath normal_path;  // Normal-only node'lardan olusan path
    BubblePath tumor_path;   // Tumor-only node'lardan olusan path
};

struct Node {
    node_hash hash;
    uint32_t colors;
    unsigned int count;
    std::vector<Occurrence> occurrences;
    
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
    std::vector<std::string> seq_names;
    
    std::unordered_map<node_hash, std::unordered_set<node_hash>> next_edges;
    std::unordered_map<node_hash, std::unordered_set<node_hash>> prev_edges;

    LCPGraph() = default;
    ~LCPGraph() = default;
    
    void addString(struct lps* lcp_string, int color, const char* sequence, uint32_t seq_id);
    void exportGFA(const std::string& output_dir, const std::string& experiment_name, const std::map<uint32_t, uint32_t>* color_map = nullptr) const;
    size_t nodeCount() const { return nodes.size(); }
    
    BubbleResult detectBubbles() const;
    void markBubbles(const BubbleResult& bubbles);
    void markDeletions(const std::unordered_set<node_hash>& deletion_nodes);
    void compressGraph();
    
    std::vector<Bubble> extractBubbles(
        const BubbleResult& bubble_markers,
        uint32_t ref_seq_id,
        uint32_t normal_seq_id,
        uint32_t tumor_seq_id
    ) const;

private:
    bool nodeExists(node_hash h) const;
    uint32_t getNodeColor(node_hash h) const;
    const std::unordered_set<node_hash>& getNextNodes(node_hash h) const;
    const std::unordered_set<node_hash>& getPrevNodes(node_hash h) const;
    bool isChainStart(node_hash h, const std::unordered_set<node_hash>& claimed) const;
};

#endif // GRAPH_H
