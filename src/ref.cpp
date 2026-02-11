#include <queue>
#include <unordered_set>

#include "ref.h"
#include "utils.h"

static Logger ref_logger("ref");

void RefIndex::addString(struct lps* l, const char* seq) {
    for (int i = 0; i < l->size; ++i) {
        core& c = l->cores[i];
        std::string region(seq + c.start, c.end - c.start);
        node_hash h = hash_string(region, 0);
        
        size_t idx = cores.size();
        cores.push_back({h, (uint64_t)c.start, (uint64_t)c.end});
        index[h].push_back(idx);
    }
    ref_logger.debug("Indexed " + std::to_string(cores.size()) + " reference cores");
}

int64_t RefIndex::anchor(node_hash h) const {
    auto it = index.find(h);
    if (it == index.end() || it->second.size() != 1) return -1;
    return cores[it->second[0]].start;
}

int64_t RefIndex::anchorBubble(const Bubble& bubble, const LCPGraph& graph) const {
    node_hash start = bubble.start_node;
    
    // Try start_node directly
    int64_t pos = anchor(start);
    if (pos >= 0) {
        auto it = index.find(start);
        size_t idx = it->second[0];
        return cores[idx].end;
    }
    
    // Walk backward to find unique anchor
    std::queue<node_hash> q;
    std::unordered_set<node_hash> seen;
    q.push(start);
    seen.insert(start);
    
    while (!q.empty()) {
        node_hash cur = q.front(); q.pop();
        
        auto prev_it = graph.prev_edges.find(cur);
        if (prev_it == graph.prev_edges.end()) continue;
        
        for (node_hash ph : prev_it->second) {
            if (seen.count(ph)) continue;
            seen.insert(ph);
            
            int64_t anchor_pos = anchor(ph);
            if (anchor_pos >= 0) {
                // Found unique anchor, now find start_node's position in reference
                auto ait = index.find(ph);
                size_t anchor_idx = ait->second[0];
                
                // Look forward in reference cores to find start_node
                for (size_t i = anchor_idx + 1; i < cores.size(); i++) {
                    if (cores[i].hash == start) return cores[i].end;
                }
                return -1;
            }
            q.push(ph);
        }
    }
    
    // Fallback: walk forward from end_node
    node_hash end = bubble.end_node;
    pos = anchor(end);
    if (pos >= 0) {
        auto it = index.find(end);
        size_t idx = it->second[0];
        return cores[idx].start - bubble.normal_seq.size();
    }
    
    return -1;
}
