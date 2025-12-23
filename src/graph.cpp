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

void LCPGraph::addString(struct lps* lcp_string, int color, const char* sequence) {
    if (lcp_string->size == 0) return;
    
    // For k=1: each core becomes a single node
    for (int i = 0; i < lcp_string->size; ++i) {
        node_hash_t hash = hash_core(&lcp_string->cores[i], sequence);
        
        Node& node = nodes[hash];
        if (node.hash == 0) node.hash = hash;
        node.count++;
        node.setColor(color);
        
        node_positions[hash].push_back({lcp_string->cores[i].start, lcp_string->cores[i].end});

        if (i + 1 < lcp_string->size) {
            node_hash_t next_hash = hash_core(&lcp_string->cores[i + 1], sequence);
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
    std::unordered_map<node_hash_t, size_t> hash_to_id;
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
            for (const auto& [start, end] : pos_it->second) {
                if (end > start) {
                    node_length += (end - start);
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
        for (node_hash_t next_hash : next_set) {
            if (hash_to_id.count(next_hash) && next_hash != hash) {
                size_t to_id = hash_to_id[next_hash];
                gfa_out << "L\t" << from_id << "\t+\t" << to_id << "\t+\t0M\n";
            }
        }
    }
}

void LCPGraph::compressGraph(std::unordered_set<node_hash_t>* deletion_nodes) {
    graph_logger.debug("Starting graph compression...");
    
    // First, mark all deletion nodes with color 2 (deletion color)
    // Treat deletion as a separate category: N=0, T=1, D=2
    if (deletion_nodes) {
        for (node_hash_t del_hash : *deletion_nodes) {
            auto it = nodes.find(del_hash);
            if (it != nodes.end()) {
                it->second.setColor(2);  // Color 2 = deletion (GT deletion)
            }
        }
    }
    
    // Track which nodes are deletion nodes (will be updated as we collapse)
    std::unordered_set<node_hash_t> is_deletion;
    if (deletion_nodes) {
        is_deletion = *deletion_nodes;
    }

    // --- Phase 1: discover disjoint maximal chains to collapse ---
    std::unordered_set<node_hash_t> claimed; // ensures no overlapping chains
    std::vector<std::vector<node_hash_t>> chains;

    // Snapshot keys so we don't depend on map iteration while mutating later
    std::vector<node_hash_t> keys;
    keys.reserve(nodes.size());
    for (const auto& kv : nodes) keys.push_back(kv.first);

    auto exists = [&](node_hash_t h) -> bool {
        return nodes.find(h) != nodes.end();
    };

    auto colorOf = [&](node_hash_t h) -> uint32_t {
        auto it = nodes.find(h);
        return (it == nodes.end()) ? 0u : it->second.colors;
    };

    // A node is a chain start if it can extend forward, but cannot be extended backward
    auto isChainStart = [&](node_hash_t h) -> bool {
        if (claimed.count(h) || !exists(h)) return false;

        const uint32_t c = colorOf(h);
        const auto& out = getNext(h);
        if (out.size() != 1) return false;

        node_hash_t n = *out.begin();
        if (!exists(n) || colorOf(n) != c) return false;

        const auto& nPrev = getPrev(n);
        if (nPrev.size() != 1) return false;
        if (*nPrev.begin() != h) return false; // robust: ensure the unique prev is h

        const auto& in = getPrev(h);
        if (in.size() != 1) return true; // no unique predecessor => start

        node_hash_t p = *in.begin();
        if (!exists(p) || colorOf(p) != c) return true;

        const auto& pOut = getNext(p);
        if (pOut.size() != 1) return true;
        if (*pOut.begin() != h) return true;

        // p -> h exists and p is also "chain-like", so h is not a start
        return false;
    };

    for (node_hash_t h : keys) {
        if (!isChainStart(h)) continue;

        const uint32_t c = colorOf(h);

        std::vector<node_hash_t> chain;
        node_hash_t cur = h;

        while (true) {
            if (!exists(cur) || claimed.count(cur)) break;

            chain.push_back(cur);
            claimed.insert(cur);

            const auto& out = getNext(cur);
            if (out.size() != 1) break;

            node_hash_t nxt = *out.begin();
            if (!exists(nxt) || claimed.count(nxt)) break;
            if (colorOf(nxt) != c) break;

            const auto& nxtPrev = getPrev(nxt);
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

        const node_hash_t first = chain.front();
        const node_hash_t last  = chain.back();

        if (!exists(first) || !exists(last)) continue;

        std::unordered_set<node_hash_t> chain_set(chain.begin(), chain.end());

        // Merge positions, counts, and colors
        std::vector<std::pair<uint64_t, uint64_t>> merged_positions;
        uint64_t merged_count = 0;
        uint32_t merged_colors = 0;  // Merge all colors from chain
        {
            size_t total_pos = 0;
            for (node_hash_t h : chain) {
                auto pit = node_positions.find(h);
                if (pit != node_positions.end()) total_pos += pit->second.size();
            }
            merged_positions.reserve(total_pos);

            for (node_hash_t h : chain) {
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
        std::unordered_set<node_hash_t> external_preds;
        if (auto it = prev_edges.find(first); it != prev_edges.end()) {
            for (node_hash_t p : it->second) {
                if (chain_set.count(p) == 0) external_preds.insert(p);
            }
        }

        std::unordered_set<node_hash_t> external_succs;
        if (auto it = next_edges.find(last); it != next_edges.end()) {
            for (node_hash_t s : it->second) {
                if (chain_set.count(s) == 0) external_succs.insert(s);
            }
        }

        // Detach all incident edges for every node in the chain
        // (Then we rebuild only what we want on the surviving node.)
        for (node_hash_t n : chain) {
            if (auto pit = prev_edges.find(n); pit != prev_edges.end()) {
                for (node_hash_t p : pit->second) {
                    if (auto nit = next_edges.find(p); nit != next_edges.end()) {
                        nit->second.erase(n);
                    }
                }
            }
            if (auto nit = next_edges.find(n); nit != next_edges.end()) {
                for (node_hash_t s : nit->second) {
                    if (auto pit = prev_edges.find(s); pit != prev_edges.end()) {
                        pit->second.erase(n);
                    }
                }
            }
        }

        // Remove nodes except the first
        for (size_t i = 1; i < chain.size(); ++i) {
            node_hash_t n = chain[i];
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

        for (node_hash_t p : external_preds) {
            if (!exists(p)) continue;
            next_edges[p].insert(first);
            prev_edges[first].insert(p);
        }

        for (node_hash_t s : external_succs) {
            if (!exists(s)) continue;
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
