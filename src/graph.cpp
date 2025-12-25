#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>

#include "graph.h"
#include "utils.h"

static Logger graph_logger("graph");

void LCPGraph::addString(struct lps* lcp_string, int color, const char* sequence, uint32_t seq_id) {
    if (lcp_string->size == 0) return;
    
    for (int i = 0; i < lcp_string->size; ++i) {
        node_hash hash = hash_core(&lcp_string->cores[i], sequence);
        
        Node& node = nodes[hash];
        if (node.hash == 0) node.hash = hash;
        node.count++;
        node.setColor(color);
        
        // Record occurrence
        node.occurrences.push_back({
            seq_id,
            lcp_string->cores[i].start,
            lcp_string->cores[i].end
        });

        if (i + 1 < lcp_string->size) {
            node_hash next_hash = hash_core(&lcp_string->cores[i + 1], sequence);
            next_edges[hash].insert(next_hash);
            prev_edges[next_hash].insert(hash);
        }
    }
}

BubbleResult LCPGraph::detectBubbles() const {
    BubbleResult result;
    
    for (const auto& [hash, node] : nodes) {
        // Check for bubble start: out-degree > 1 and outgoing edges go to different color paths
        auto next_it = next_edges.find(hash);
        if (next_it != next_edges.end() && next_it->second.size() > 1) {
            bool has_normal_path = false;
            bool has_tumor_path = false;
            
            for (node_hash next_h : next_it->second) {
                auto it = nodes.find(next_h);
                if (it != nodes.end()) {
                    uint32_t c = it->second.colors;
                    // Color 0 = Ref, Color 1 = Normal, Color 2 = Tumor
                    if ((c & 0x2) && !(c & 0x4)) has_normal_path = true;  // Normal only
                    if ((c & 0x4) && !(c & 0x2)) has_tumor_path = true;   // Tumor only
                }
            }
            
            if (has_normal_path && has_tumor_path) {
                result.bubble_starts.insert(hash);
            }
        }
        
        // Check for bubble end: in-degree > 1 and incoming edges come from different color paths
        auto prev_it = prev_edges.find(hash);
        if (prev_it != prev_edges.end() && prev_it->second.size() > 1) {
            bool has_normal_path = false;
            bool has_tumor_path = false;
            
            for (node_hash prev_h : prev_it->second) {
                auto it = nodes.find(prev_h);
                if (it != nodes.end()) {
                    uint32_t c = it->second.colors;
                    if ((c & 0x2) && !(c & 0x4)) has_normal_path = true;
                    if ((c & 0x4) && !(c & 0x2)) has_tumor_path = true;
                }
            }
            
            if (has_normal_path && has_tumor_path) {
                result.bubble_ends.insert(hash);
            }
        }
    }
    
    graph_logger.debug("Detected " + std::to_string(result.bubble_starts.size()) + " bubble starts, " 
                       + std::to_string(result.bubble_ends.size()) + " bubble ends");
    
    return result;
}

void LCPGraph::markBubbles(const BubbleResult& bubbles) {
    // Color 3 = Bubble Start, Color 4 = Bubble End
    for (node_hash h : bubbles.bubble_starts) {
        auto it = nodes.find(h);
        if (it != nodes.end()) {
            it->second.setColor(3);
        }
    }
    
    for (node_hash h : bubbles.bubble_ends) {
        auto it = nodes.find(h);
        if (it != nodes.end()) {
            it->second.setColor(4);
        }
    }
}

void LCPGraph::markDeletions(const std::unordered_set<node_hash>& deletion_nodes) {
    // Color 5 = Deletion
    for (node_hash h : deletion_nodes) {
        auto it = nodes.find(h);
        if (it != nodes.end()) {
            it->second.setColor(5);
        }
    }
}

void LCPGraph::exportGFA(const std::string& output_dir, const std::string& experiment_name, const std::map<uint32_t, uint32_t>* color_map) const {
    namespace fs = std::filesystem;
    fs::path base_dir = output_dir;
    fs::path gfa_path = base_dir / "gfa" / (experiment_name + ".gfa");
    fs::path csv_path = base_dir / "csv" / (experiment_name + ".csv");
    
    fs::create_directories(gfa_path.parent_path());
    fs::create_directories(csv_path.parent_path());
    
    std::ofstream gfa_out(gfa_path);
    if (!gfa_out) return;
    
    gfa_out << "H\tVN:Z:1.0\n";
    
    std::unordered_map<node_hash, size_t> hash_to_id;
    size_t node_id = 1;
    for (const auto& [hash, _] : nodes) {
        hash_to_id[hash] = node_id++;
    }
    
    std::map<std::string, std::vector<std::string>> csv_data;
    
    for (const auto& [hash, node] : nodes) {
        size_t id = hash_to_id[hash];
        
        uint64_t node_length = node.count;
        
        gfa_out << "S\t" << id << "\t*\tLN:i:" << node_length << "\tCN:i:" << node.count << "\tHASH:i:" << hash;
        for (int i = 0; i < 32; ++i) {
            if (node.hasColor(i)) {
                gfa_out << "\tC" << i << ":i:1";
            }
        }
        gfa_out << "\n";
        
        uint32_t c = 0;
        if (color_map) {
            auto it = color_map->find(node.colors);
            if (it != color_map->end()) {
                c = it->second;
            }
        }
        
        c &= 0xFFFFFF;
        char color_buf[9];
        snprintf(color_buf, sizeof(color_buf), "#%06X", c);
        csv_data[std::to_string(id)].push_back(color_buf);
    }

    for (const auto& [hash, next_set] : next_edges) {
        size_t from_id = hash_to_id[hash];
        for (node_hash next_hash : next_set) {
            if (hash_to_id.count(next_hash) && next_hash != hash) {
                size_t to_id = hash_to_id[next_hash];
                gfa_out << "L\t" << from_id << "\t+\t" << to_id << "\t+\t0M\n";
            }
        }
    }
    
    writeCSV(csv_data, csv_path.string(), {"Name", "Colour"}, ',');
}

void LCPGraph::compressGraph() {
    graph_logger.debug("Starting graph compression...");

    std::unordered_set<node_hash> claimed;
    std::vector<std::vector<node_hash>> chains;

    std::vector<node_hash> keys;
    keys.reserve(nodes.size());
    for (const auto& kv : nodes) keys.push_back(kv.first);

    for (node_hash h : keys) {
        if (!isChainStart(h, claimed)) continue;

        const uint32_t c = getNodeColor(h);

        std::vector<node_hash> chain;
        node_hash cur = h;

        while (true) {
            if (!nodeExists(cur) || claimed.count(cur)) break;

            chain.push_back(cur);
            claimed.insert(cur);

            const auto& curNext = getNextNodes(cur);
            if (curNext.size() != 1) break;

            node_hash nxt = *curNext.begin();
            if (!nodeExists(nxt) || claimed.count(nxt)) break;
            if (getNodeColor(nxt) != c) break;

            const auto& nxtPrev = getPrevNodes(nxt);
            if (nxtPrev.size() != 1) break;
            if (*nxtPrev.begin() != cur) break;

            cur = nxt;
        }

        if (chain.size() > 1) chains.push_back(std::move(chain));
    }

    graph_logger.debug("Found " + std::to_string(chains.size()) + " paths to collapse");

    for (const auto& chain : chains) {
        if (chain.size() < 2) continue;

        const node_hash first = chain.front();
        const node_hash last  = chain.back();

        if (!nodeExists(first) || !nodeExists(last)) continue;

        std::unordered_set<node_hash> chain_set(chain.begin(), chain.end());

        uint64_t merged_count = 0;
        uint32_t chain_color = 0;
        std::vector<Occurrence> merged_occurrences;
        
        for (node_hash h : chain) {
            if (auto nit = nodes.find(h); nit != nodes.end()) {
                merged_count += nit->second.count;
                if (chain_color == 0) {
                    chain_color = nit->second.colors;
                }
                // Merge occurrences
                for (const auto& occ : nit->second.occurrences) {
                    merged_occurrences.push_back(occ);
                }
            }
        }

        std::unordered_set<node_hash> external_preds;
        if (auto it = prev_edges.find(first); it != prev_edges.end()) {
            for (node_hash p : it->second) {
                if (chain_set.count(p) == 0) external_preds.insert(p);
            }
        }

        std::unordered_set<node_hash> external_succs;
        if (auto it = next_edges.find(last); it != next_edges.end()) {
            for (node_hash s : it->second) {
                if (chain_set.count(s) == 0) external_succs.insert(s);
            }
        }

        for (node_hash n : chain) {
            if (auto pit = prev_edges.find(n); pit != prev_edges.end()) {
                for (node_hash p : pit->second) {
                    if (auto nit = next_edges.find(p); nit != next_edges.end()) {
                        nit->second.erase(n);
                    }
                }
            }
            if (auto nit = next_edges.find(n); nit != next_edges.end()) {
                for (node_hash s : nit->second) {
                    if (auto pit = prev_edges.find(s); pit != prev_edges.end()) {
                        pit->second.erase(n);
                    }
                }
            }
        }

        for (size_t i = 1; i < chain.size(); ++i) {
            node_hash n = chain[i];
            nodes.erase(n);
            next_edges.erase(n);
            prev_edges.erase(n);
        }

        prev_edges[first].clear();
        next_edges[first].clear();
        if (auto fit = nodes.find(first); fit != nodes.end()) {
            fit->second.count = merged_count;
            fit->second.colors = chain_color;
            fit->second.occurrences = std::move(merged_occurrences);
        }

        for (node_hash p : external_preds) {
            if (!nodeExists(p)) continue;
            next_edges[p].insert(first);
            prev_edges[first].insert(p);
        }

        for (node_hash s : external_succs) {
            if (!nodeExists(s)) continue;
            prev_edges[s].insert(first);
            next_edges[first].insert(s);
        }
    }

    graph_logger.debug("Graph compression completed. Nodes: " + std::to_string(nodes.size()));
}

bool LCPGraph::nodeExists(node_hash h) const {
    return nodes.find(h) != nodes.end();
}

uint32_t LCPGraph::getNodeColor(node_hash h) const {
    auto it = nodes.find(h);
    return (it == nodes.end()) ? 0u : it->second.colors;
}

const std::unordered_set<node_hash>& LCPGraph::getNextNodes(node_hash h) const {
    static const std::unordered_set<node_hash> empty_set;
    auto it = next_edges.find(h);
    return (it != next_edges.end()) ? it->second : empty_set;
}

const std::unordered_set<node_hash>& LCPGraph::getPrevNodes(node_hash h) const {
    static const std::unordered_set<node_hash> empty_set;
    auto it = prev_edges.find(h);
    return (it != prev_edges.end()) ? it->second : empty_set;
}

bool LCPGraph::isChainStart(node_hash h, const std::unordered_set<node_hash>& claimed) const {
    if (claimed.count(h) || !nodeExists(h)) return false;

    const uint32_t c = getNodeColor(h);
    const auto& out = getNextNodes(h);
    if (out.size() != 1) return false;

    node_hash n = *out.begin();
    if (!nodeExists(n) || getNodeColor(n) != c) return false;

    const auto& nPrev = getPrevNodes(n);
    if (nPrev.size() != 1) return false;
    if (*nPrev.begin() != h) return false;

    const auto& in = getPrevNodes(h);
    if (in.size() != 1) return true;

    node_hash p = *in.begin();
    if (!nodeExists(p) || getNodeColor(p) != c) return true;

    const auto& pOut = getNextNodes(p);
    if (pOut.size() != 1) return true;
    if (*pOut.begin() != h) return true;

    return false;
}

// Helper to get occurrence for a specific seq_id
static std::pair<uint64_t, uint64_t> getOccurrenceCoords(const Node& node, uint32_t seq_id) {
    for (const auto& occ : node.occurrences) {
        if (occ.seq_id == seq_id) {
            return {occ.start, occ.end};
        }
    }
    return {UINT64_MAX, UINT64_MAX};
}

std::vector<Bubble> LCPGraph::extractBubbles(
    const BubbleResult& bubble_markers,
    uint32_t ref_seq_id,
    uint32_t normal_seq_id,
    uint32_t tumor_seq_id
) const {
    std::vector<Bubble> bubbles;
    
    // For each bubble start, try to find corresponding bubble end
    for (node_hash start_h : bubble_markers.bubble_starts) {
        auto start_it = nodes.find(start_h);
        if (start_it == nodes.end()) continue;
        
        const Node& start_node = start_it->second;
        
        // Get reference coordinates from start node
        auto [ref_start, ref_start_end] = getOccurrenceCoords(start_node, ref_seq_id);
        auto [tumor_start, tumor_start_end] = getOccurrenceCoords(start_node, tumor_seq_id);
        
        // If no ref occurrence, try normal (same as ref)
        if (ref_start == UINT64_MAX) {
            auto [ns, ne] = getOccurrenceCoords(start_node, normal_seq_id);
            ref_start = ns;
            ref_start_end = ne;
        }
        
        if (ref_start == UINT64_MAX) continue;  // Can't anchor this bubble
        
        // BFS to find bubble end following normal-only path
        std::vector<node_hash> normal_path_nodes;
        std::vector<node_hash> tumor_path_nodes;
        node_hash end_node_hash = 0;
        
        // Find normal-only path (nodes with color 1 set but not color 2)
        {
            std::queue<node_hash> q;
            std::unordered_set<node_hash> visited;
            
            for (node_hash next_h : getNextNodes(start_h)) {
                auto it = nodes.find(next_h);
                if (it == nodes.end()) continue;
                uint32_t c = it->second.colors;
                // Normal-only: has normal (bit 1), no tumor (bit 2)
                if ((c & 0x2) && !(c & 0x4)) {
                    q.push(next_h);
                    visited.insert(next_h);
                }
            }
            
            while (!q.empty()) {
                node_hash cur = q.front();
                q.pop();
                
                // Check if this is a bubble end
                if (bubble_markers.bubble_ends.count(cur)) {
                    end_node_hash = cur;
                    break;
                }
                
                normal_path_nodes.push_back(cur);
                
                // Continue along normal-only path
                for (node_hash next_h : getNextNodes(cur)) {
                    if (visited.count(next_h)) continue;
                    
                    auto it = nodes.find(next_h);
                    if (it == nodes.end()) continue;
                    
                    // If it's a bubble end, accept it regardless of color
                    if (bubble_markers.bubble_ends.count(next_h)) {
                        end_node_hash = next_h;
                        break;
                    }
                    
                    uint32_t c = it->second.colors;
                    // Continue if normal-only or shared (to reach bubble end)
                    if (c & 0x2) {  // Has normal
                        q.push(next_h);
                        visited.insert(next_h);
                    }
                }
                
                if (end_node_hash != 0) break;
            }
        }
        
        if (end_node_hash == 0) continue;  // No bubble end found
        
        // Find tumor-only path
        {
            std::queue<node_hash> q;
            std::unordered_set<node_hash> visited;
            
            for (node_hash next_h : getNextNodes(start_h)) {
                auto it = nodes.find(next_h);
                if (it == nodes.end()) continue;
                uint32_t c = it->second.colors;
                // Tumor-only: has tumor (bit 2), no normal (bit 1)
                if ((c & 0x4) && !(c & 0x2)) {
                    q.push(next_h);
                    visited.insert(next_h);
                }
            }
            
            while (!q.empty()) {
                node_hash cur = q.front();
                q.pop();
                
                if (cur == end_node_hash) break;
                
                tumor_path_nodes.push_back(cur);
                
                for (node_hash next_h : getNextNodes(cur)) {
                    if (visited.count(next_h)) continue;
                    if (next_h == end_node_hash) break;
                    
                    auto it = nodes.find(next_h);
                    if (it == nodes.end()) continue;
                    
                    uint32_t c = it->second.colors;
                    if (c & 0x4) {  // Has tumor
                        q.push(next_h);
                        visited.insert(next_h);
                    }
                }
            }
        }
        
        // Get end node coordinates
        auto end_it = nodes.find(end_node_hash);
        if (end_it == nodes.end()) continue;
        
        const Node& end_node = end_it->second;
        auto [ref_end_start, ref_end] = getOccurrenceCoords(end_node, ref_seq_id);
        auto [tumor_end_start, tumor_end] = getOccurrenceCoords(end_node, tumor_seq_id);
        
        if (ref_end_start == UINT64_MAX) {
            auto [ns, ne] = getOccurrenceCoords(end_node, normal_seq_id);
            ref_end_start = ns;
            ref_end = ne;
        }
        
        if (ref_end == UINT64_MAX) continue;
        
        // Create bubble
        Bubble bubble;
        bubble.start_node = start_h;
        bubble.end_node = end_node_hash;
        bubble.ref_start = ref_start_end;    // End of start node
        bubble.ref_end = ref_end_start;       // Start of end node
        bubble.tumor_start = (tumor_start != UINT64_MAX) ? tumor_start_end : 0;
        bubble.tumor_end = (tumor_end_start != UINT64_MAX) ? tumor_end_start : 0;
        
        bubble.normal_path.nodes = std::move(normal_path_nodes);
        bubble.normal_path.seq_id = normal_seq_id;
        
        bubble.tumor_path.nodes = std::move(tumor_path_nodes);
        bubble.tumor_path.seq_id = tumor_seq_id;
        
        bubbles.push_back(std::move(bubble));
    }
    
    graph_logger.debug("Extracted " + std::to_string(bubbles.size()) + " complete bubbles");
    
    return bubbles;
}
