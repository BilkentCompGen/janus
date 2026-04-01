#include <fstream>
#include <iomanip>
#include <queue>
#include <sstream>

#include "graph.h"
#include "utils.h"

static Logger graph_logger("graph");

void LCPGraph::addString(struct lps* l, int color, const char* seq, uint32_t) {
    for (int i = 0; i < l->size; ++i) {
        core& c = l->cores[i];
        std::string region(seq + c.start, c.end - c.start);
        node_hash h = hash_string(region, 0);
        
        Node& n = nodes[h];
        if (n.hash == 0) n.hash = h, n.sequence = region;
        n.setColor(color);

        if (i + 1 < l->size) {
            core& nc = l->cores[i + 1];
            std::string next_region(seq + nc.start, nc.end - nc.start);
            node_hash nh = hash_string(next_region, 0);
            next_edges[h].insert(nh);
            prev_edges[nh].insert(h);
        }
    }
}

BubbleResult LCPGraph::detectBubbles() const {
    BubbleResult result;

    for (const auto& [hash, node] : nodes) {
        bool node_has_normal = (node.colors & 0x1) != 0;
        bool node_has_tumor = (node.colors & 0x2) != 0;
        if (!(node_has_normal && node_has_tumor)) continue;

        auto next_it = next_edges.find(hash);
        if (next_it != next_edges.end()) {
            for (node_hash nh : next_it->second) {
                auto it = nodes.find(nh);
                if (it == nodes.end()) continue;
                bool has_normal = (it->second.colors & 0x1) != 0;
                bool has_tumor = (it->second.colors & 0x2) != 0;
                if (!(has_normal && has_tumor)) {
                    result.bubble_starts.insert(hash);
                    break;
                }
            }
        }

        auto prev_it = prev_edges.find(hash);
        if (prev_it != prev_edges.end()) {
            for (node_hash ph : prev_it->second) {
                auto it = nodes.find(ph);
                if (it == nodes.end()) continue;
                bool has_normal = (it->second.colors & 0x1) != 0;
                bool has_tumor = (it->second.colors & 0x2) != 0;
                if (!(has_normal && has_tumor)) {
                    result.bubble_ends.insert(hash);
                    break;
                }
            }
        }
    }

    graph_logger.debug(
        "Detected "
        + std::to_string(result.bubble_starts.size()) + " bubble starts, "
        + std::to_string(result.bubble_ends.size()) + " bubble ends"
    );
    return result;
}

node_hash LCPGraph::findColoredPath(node_hash start, const BubbleResult& markers, int color_bit, Bubble& bubble) const {
    uint32_t target = 1 << color_bit, avoid = 1 << (1 - color_bit);
    std::queue<std::pair<node_hash, Node>> q;

    auto it_start = next_edges.find(start);
    if (it_start == next_edges.end()) return 0;

    const Node& start_node = nodes.at(start);
    size_t start_len = start_node.sequence.size();

    for (node_hash n : it_start->second) {
        if (markers.bubble_ends.count(n)) {
            (color_bit == 0 ? bubble.normal_seq : bubble.tumor_seq) = "";
            return n;
        }
        auto it_node = nodes.find(n);
        if (it_node != nodes.end() && (it_node->second.colors & target) && !(it_node->second.colors & avoid)) {
            Node acc = start_node;
            acc.merge(it_node->second);
            q.push({n, acc});
        }
    }

    std::unordered_set<node_hash> seen;
    while (!q.empty()) {
        auto [cur, acc] = q.front(); q.pop();
        if (!seen.insert(cur).second) continue;

        auto it_next = next_edges.find(cur);
        if (it_next != next_edges.end()) {
            for (node_hash n : it_next->second) {
                if (markers.bubble_ends.count(n)) {
                    size_t before = acc.sequence.size();
                    acc.merge(nodes.at(n));
                    size_t end_overlap = nodes.at(n).sequence.size() - (acc.sequence.size() - before);
                    if (before >= start_len + end_overlap) {
                        (color_bit == 0 ? bubble.normal_seq : bubble.tumor_seq) = acc.sequence.substr(start_len, before - start_len - end_overlap);
                    } else {
                        (color_bit == 0 ? bubble.normal_seq : bubble.tumor_seq) = "";
                    }
                    return n;
                }
                auto it_n = nodes.find(n);
                if (it_n != nodes.end() && (it_n->second.colors & target) && !(it_n->second.colors & avoid)) {
                    Node next_acc = acc;
                    next_acc.merge(it_n->second);
                    q.push({n, next_acc});
                }
            }
        }
    }
    return 0;
}

std::vector<Bubble> LCPGraph::extractBubbles(const BubbleResult& markers) const {
    std::vector<Bubble> bubbles;
    for (node_hash start : markers.bubble_starts) {
        Bubble bubble;
        bubble.start_node = start;
        node_hash end_normal = findColoredPath(start, markers, 0, bubble);
        node_hash end_tumor = findColoredPath(start, markers, 1, bubble);

        if (end_normal != 0 && end_normal == end_tumor) {
            bubble.end_node = end_normal;
            bubbles.push_back(bubble);
        }
    }
    graph_logger.debug("Extracted " + std::to_string(bubbles.size()) + " bubbles");
    return bubbles;
}

void LCPGraph::writeGFA(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) {
        graph_logger.debug("Failed to open file for GFA output: " + filename);
        return;
    }
    writeGFA(out);
    graph_logger.debug("Wrote GFA to " + filename);
}

void LCPGraph::writeGFA(std::ostream& out) const {
    // Build a stable mapping from node_hash -> sequential uint64_t ID (starting at 1)
    std::unordered_map<node_hash, uint64_t> id_map;
    uint64_t next_id = 1;
    for (const auto& [hash, node] : nodes) {
        id_map[hash] = next_id++;
    }

    auto hashToName = [&id_map](node_hash h) -> std::string {
        auto it = id_map.find(h);
        if (it != id_map.end()) return std::to_string(it->second);
        return "0";
    };

    out << "H\tVN:Z:1.1\n";
    
    out << "# LCPGraph export - contains " << nodes.size() << " segments, " << seq_names.size() << " source sequences\n";
    
    for (size_t i = 0; i < seq_names.size(); ++i) {
        out << "# seq[" << i << "]: " << seq_names[i] << "\n";
    }
    
    for (const auto& [hash, node] : nodes) {
        std::string name = hashToName(hash);
        
        out << "S\t" << name << "\t" << node.sequence;
        
        out << "\tHA:H:";
        auto flags = out.flags();
        auto fill = out.fill();
        out << std::hex << std::uppercase << std::setfill('0') << std::setw(16) << hash;
        out.flags(flags);
        out.fill(fill);
        
        out << "\tCL:i:" << node.colors;
        
        out << "\tLN:i:" << node.sequence.length();
        out << "\n";
    }
    
    for (const auto& [from_hash, to_set] : next_edges) {
        std::string from_name = hashToName(from_hash);
        for (node_hash to_hash : to_set) {
            std::string to_name = hashToName(to_hash);
            out << "L\t" << from_name << "\t+\t" << to_name << "\t+\t0M\n";
        }
    }
    
    for (size_t color = 0; color < seq_names.size() && color < 32; ++color) {
        uint32_t color_mask = 1U << color;
        
        std::vector<node_hash> color_nodes;
        for (const auto& [hash, node] : nodes) {
            if (node.colors & color_mask) {
                color_nodes.push_back(hash);
            }
        }
        
        if (color_nodes.empty()) continue;
        
        std::vector<node_hash> entry_nodes;
        for (node_hash h : color_nodes) {
            auto prev_it = prev_edges.find(h);
            bool has_colored_pred = false;
            if (prev_it != prev_edges.end()) {
                for (node_hash ph : prev_it->second) {
                    auto it = nodes.find(ph);
                    if (it != nodes.end() && (it->second.colors & color_mask)) {
                        has_colored_pred = true;
                        break;
                    }
                }
            }
            if (!has_colored_pred) {
                entry_nodes.push_back(h);
            }
        }
        
        std::unordered_set<node_hash> visited;
        int walk_idx = 0;
        for (node_hash start : entry_nodes) {
            if (visited.count(start)) continue;
            
            std::ostringstream walk;
            node_hash cur = start;
            uint64_t walk_len = 0;
            
            while (cur != 0 && !visited.count(cur)) {
                visited.insert(cur);
                walk << ">" << hashToName(cur);
                
                auto it = nodes.find(cur);
                if (it != nodes.end()) {
                    walk_len += it->second.sequence.length();
                }
                
                node_hash next = 0;
                auto next_it = next_edges.find(cur);
                if (next_it != next_edges.end()) {
                    for (node_hash nh : next_it->second) {
                        auto nit = nodes.find(nh);
                        if (nit != nodes.end() && (nit->second.colors & color_mask)) {
                            next = nh;
                            break;
                        }
                    }
                }
                cur = next;
            }
            
            std::string sample = color < seq_names.size() ? seq_names[color] : "sample" + std::to_string(color);
            std::string seq_id = "seq" + std::to_string(walk_idx++);
            
            out << "W\t" << sample << "\t0\t" << seq_id << "\t0\t" << walk_len << "\t" << walk.str() << "\n";
        }
    }
}