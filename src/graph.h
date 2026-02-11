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

struct BubbleResult {
    std::unordered_set<node_hash> bubble_starts;
    std::unordered_set<node_hash> bubble_ends;
};

struct Bubble {
    node_hash start_node, end_node;
    std::string normal_seq;
    std::string tumor_seq;
};

struct Node {
    node_hash hash;
    uint32_t colors;
    std::string sequence;
    
    Node() : hash(0), colors(0) {}
    Node(node_hash h) : hash(h), colors(0) {}
    
    void merge(const Node& next) { // TODO: I am gonna check this implementation for efficieny
        const std::string& pat = next.sequence;
        size_t n = pat.size();
        if (n == 0) return;
        
        std::vector<int> pi(n);
        for (size_t i = 1, j = 0; i < n; i++) {
            while (j > 0 && pat[i] != pat[j]) j = pi[j - 1];
            if (pat[i] == pat[j]) j++;
            pi[i] = j;
        }
        
        size_t m = sequence.size();
        size_t check_len = std::min(m, n);
        
        size_t j = 0;
        for (size_t i = m - check_len; i < m; i++) {
            while (j > 0 && sequence[i] != pat[j]) j = pi[j - 1];
            if (sequence[i] == pat[j]) j++;
        }
        
        sequence += pat.substr(j);
    }
    
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
    node_hash findColoredPath(node_hash start, const BubbleResult& markers, int color_bit, Bubble& bubble) const;
    size_t nodeCount() const { return nodes.size(); }
    
    BubbleResult detectBubbles() const;
    std::vector<Bubble> extractBubbles(const BubbleResult& bubble_markers) const;
    
    void writeGFA(const std::string& filename) const;
    void writeGFA(std::ostream& out) const;
};

#endif // GRAPH_H
