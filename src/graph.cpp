#include "graph.h"
#include "sv.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <queue>
#include <set>

static Logger graph_logger("graph");

void LCPGraph::addString(struct lps* lcp_string, int color, const char* sequence, const std::string& chrom) {
    if (lcp_string->size == 0) return;
    
    // For k=1: each core becomes a single node
    for (int i = 0; i < lcp_string->size; ++i) {
        node_hash hash = hash_core(&lcp_string->cores[i], sequence);
        
        Node& node = nodes[hash];
        if (node.hash == 0) node.hash = hash;
        node.count++;
        node.setColor(color);
        
        // Store position with cell_id, chrom, start, end
        node_positions[hash].emplace_back(color, chrom, lcp_string->cores[i].start, lcp_string->cores[i].end);

        if (i + 1 < lcp_string->size) {
            node_hash next_hash = hash_core(&lcp_string->cores[i + 1], sequence);
            next_edges[hash].insert(next_hash);
            prev_edges[next_hash].insert(hash);
        }
    }
}

void LCPGraph::exportGFA(const std::string& output_dir, const std::string& experiment_name) const {
    namespace fs = std::filesystem;
    fs::path base_dir = output_dir;
    fs::path gfa_path = base_dir / "gfa" / (experiment_name + ".gfa");
    fs::path csv_path = base_dir / "csv" / (experiment_name + ".csv");
    
    fs::create_directories(gfa_path.parent_path());
    fs::create_directories(csv_path.parent_path());
    
    std::ofstream gfa_out(gfa_path);
    if (!gfa_out) return;
    
    std::ofstream csv_out(csv_path);
    if (!csv_out) return;
    
    gfa_out << "H\tVN:Z:1.0\n";
    csv_out << "Name,Colour\n";
    
    // Map hash to sequential ID
    std::unordered_map<node_hash, size_t> hash_to_id;
    size_t node_id = 1;
    for (const auto& [hash, _] : nodes) {
        hash_to_id[hash] = node_id++;
    }
    
    // Write nodes
    for (const auto& [hash, node] : nodes) {
        size_t id = hash_to_id[hash];
        
        // Calculate node length from positions (sum of all position lengths)
        uint64_t node_length = 0;
        auto pos_it = node_positions.find(hash);
        if (pos_it != node_positions.end()) {
            for (const auto& pos_info : pos_it->second) {
                if (pos_info.end > pos_info.start) {
                    node_length += (pos_info.end - pos_info.start);
                }
            }
        }
        // If no positions found, use count as fallback
        if (node_length == 0) {
            node_length = node.count;
        }
        
        // GFA: S line
        gfa_out << "S\t" << id << "\t*\tLN:i:" << node_length << "\tCN:i:" << node.count << "\tHASH:i:" << hash;
        // Add color info
        for (int i = 0; i < 32; ++i) {
            if (node.hasColor(i)) {
                gfa_out << "\tC" << i << ":i:1";
            }
        }
        gfa_out << "\n";
        
        // CSV: color for visualization
        // 7-color system: N=0, T=1, GT=2 (deletion from ground truth)
        // N (0) → Blue (0x0000FF)
        // T (1) → Red (0xFF0000)
        // GT (2) → Green (0x00FF00)
        // N+T (0+1) → Purple (0x800080)
        // N+GT (0+2) → Cyan (0x00FFFF) - T yok, sadece N ve GT
        // T+GT (1+2) → Yellow (0xFFFF00)
        // N+T+GT (0+1+2) → Orange (0xFF8000) - T dallanıyor (bubble'lar)
        
        bool has_normal = node.hasColor(0);
        bool has_tumor = false;
        // Check for tumor color (1) specifically, not all colors > 0
        if (node.hasColor(1)) {
            has_tumor = true;
        }
        bool has_deletion = node.hasColor(2);
        
        uint32_t c = 0;
        if (has_normal && has_tumor && has_deletion) {
            c = 0xFF8000;  // Orange: N+T+GT (T dallanıyor - bubble'lar)
        } else if (has_tumor && has_deletion) {
            c = 0xFFFF00;  // Yellow: T+GT
        } else if (has_normal && has_deletion) {
            c = 0x00FFFF;  // Cyan: N+GT (T yok, sadece N ve GT)
        } else if (has_deletion) {
            c = 0x00FF00;  // Green: GT only
        } else if (has_normal && has_tumor) {
            c = 0x800080;  // Purple: N+T
        } else if (has_normal) {
            c = 0x0000FF;  // Blue: N only
        } else if (has_tumor) {
            c = 0xFF0000;  // Red: T only
        }
        // else: c = 0 (black) - should not happen
        
        c &= 0xFFFFFF;
        char color_buf[9];
        snprintf(color_buf, sizeof(color_buf), "#%06X", c);
        csv_out << id << "," << color_buf << "\n";
    }
    
    // Write edges
    for (const auto& [hash, next_set] : next_edges) {
        size_t from_id = hash_to_id[hash];
        for (node_hash next_hash : next_set) {
            if (hash_to_id.count(next_hash) && next_hash != hash) {
                size_t to_id = hash_to_id[next_hash];
                gfa_out << "L\t" << from_id << "\t+\t" << to_id << "\t+\t0M\n";
            }
        }
    }
}

void LCPGraph::compressGraph(std::unordered_set<node_hash>* deletion_nodes) {
    graph_logger.debug("Starting graph compression...");
    
    // First, mark all deletion nodes with color 2 (deletion color)
    // Treat deletion as a separate category: N=0, T=1, D=2
    if (deletion_nodes) {
        for (node_hash del_hash : *deletion_nodes) {
            auto it = nodes.find(del_hash);
            if (it != nodes.end()) {
                it->second.setColor(2);  // Color 2 = deletion (GT deletion)
            }
        }
    }
    
    // Track which nodes are deletion nodes (will be updated as we collapse)
    std::unordered_set<node_hash> is_deletion;
    if (deletion_nodes) {
        is_deletion = *deletion_nodes;
    }

    // --- Phase 1: discover disjoint maximal chains to collapse ---
    std::unordered_set<node_hash> claimed; // ensures no overlapping chains
    std::vector<std::vector<node_hash>> chains;

    // Snapshot keys so we don't depend on map iteration while mutating later
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

            const auto& out = getNextNodes(cur);
            if (out.size() != 1) break;

            node_hash nxt = *out.begin();
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

    // --- Phase 2: collapse chains (disjoint by construction) ---
    for (const auto& chain : chains) {
        if (chain.size() < 2) continue;

        const node_hash first = chain.front();
        const node_hash last  = chain.back();

        if (!nodeExists(first) || !nodeExists(last)) continue;

        std::unordered_set<node_hash> chain_set(chain.begin(), chain.end());

        // Merge positions, counts, and colors
        std::vector<PositionInfo> merged_positions;
        uint64_t merged_count = 0;
        uint32_t merged_colors = 0;  // Merge all colors from chain
        {
            size_t total_pos = 0;
            for (node_hash h : chain) {
                auto pit = node_positions.find(h);
                if (pit != node_positions.end()) total_pos += pit->second.size();
            }
            merged_positions.reserve(total_pos);

            for (node_hash h : chain) {
                if (auto nit = nodes.find(h); nit != nodes.end()) {
                    merged_count += nit->second.count;
                    merged_colors |= nit->second.colors;  // Merge colors (OR operation)
                }
                auto pit = node_positions.find(h);
                if (pit != node_positions.end()) {
                    merged_positions.insert(
                        merged_positions.end(),
                        pit->second.begin(),
                        pit->second.end()
                    );
                }
            }
        }

        // Snapshot external neighbors we want to preserve
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

        // Detach all incident edges for every node in the chain
        // (Then we rebuild only what we want on the surviving node.)
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

        // Remove nodes except the first
        for (size_t i = 1; i < chain.size(); ++i) {
            node_hash n = chain[i];
            nodes.erase(n);
            node_positions.erase(n);
            next_edges.erase(n);
            prev_edges.erase(n);
            is_deletion.erase(n);  // Remove from deletion set
        }

        // Rebuild adjacency for the collapsed node (= first)
        prev_edges[first].clear();
        next_edges[first].clear();
        if (auto fit = nodes.find(first); fit != nodes.end()) {
            fit->second.count = merged_count;
            fit->second.colors = merged_colors;  // Set merged colors (includes color 2 if any node had it)
            // Track in deletion set if has color 2
            if (fit->second.hasColor(2)) {
                is_deletion.insert(first);
            }
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

        // Store merged positions
        node_positions[first] = std::move(merged_positions);
    }
    
    // Update deletion_nodes set if provided (after all collapses)
    if (deletion_nodes) {
        *deletion_nodes = std::move(is_deletion);
    }

    graph_logger.debug("Graph compression completed. Nodes: " + std::to_string(nodes.size()));
}

// Helper method implementations
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
    if (*nPrev.begin() != h) return false; // robust: ensure the unique prev is h

    const auto& in = getPrevNodes(h);
    if (in.size() != 1) return true; // no unique predecessor => start

    node_hash p = *in.begin();
    if (!nodeExists(p) || getNodeColor(p) != c) return true;

    const auto& pOut = getNextNodes(p);
    if (pOut.size() != 1) return true;
    if (*pOut.begin() != h) return true;

    // p -> h exists and p is also "chain-like", so h is not a start
    return false;
}

std::vector<PositionInfo> LCPGraph::getPositions(node_hash hash, int cell_id) const {
    auto it = node_positions.find(hash);
    if (it == node_positions.end()) {
        return {};
    }
    
    if (cell_id == -1) {
        // Return all positions
        return it->second;
    }
    
    // Filter by cell_id
    std::vector<PositionInfo> result;
    for (const auto& pos : it->second) {
        if (pos.cell_id == cell_id) {
            result.push_back(pos);
        }
    }
    return result;
}

std::string LCPGraph::getSequence(int cell_id, const std::string& chrom) const {
    auto cell_it = sequences.find(cell_id);
    if (cell_it == sequences.end()) {
        return "";
    }
    
    auto chrom_it = cell_it->second.find(chrom);
    if (chrom_it == cell_it->second.end()) {
        return "";
    }
    
    return chrom_it->second;
}

void LCPGraph::setSequence(int cell_id, const std::string& chrom, const std::string& sequence) {
    sequences[cell_id][chrom] = sequence;
}
