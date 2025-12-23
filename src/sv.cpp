#include "sv.h"
#include "utils.h"
#include "vcf.h"
#include <algorithm>
#include <set>
#include <fstream>
#include <sstream>
#include <tuple>
#include <filesystem>
#include <map>
#include <functional>

static Logger sv_logger("sv");

namespace {
    // Helper function to compare chromosomes (numeric if possible, otherwise lexicographic)
    int compareChromosomes(const std::string& chrom1, const std::string& chrom2) {
        try {
            int num1 = std::stoi(chrom1);
            int num2 = std::stoi(chrom2);
            if (num1 != num2) return num1 < num2 ? -1 : 1;
            return 0;
        } catch (...) {
            if (chrom1 < chrom2) return -1;
            if (chrom1 > chrom2) return 1;
            return 0;
        }
    }
}

void LinearReference::addCore(node_hash_t hash, const std::string& chrom, uint64_t start, uint64_t end) {
    hashes.push_back(hash);
    chroms.push_back(chrom);
    starts.push_back(start);
    ends.push_back(end);
}

void LinearReference::buildIndex() {
    hash_positions.clear();
    for (size_t i = 0; i < hashes.size(); ++i) {
        hash_positions[hashes[i]].push_back(i);
    }
}

std::vector<size_t> LinearReference::align(node_hash_t hash) const {
    auto it = hash_positions.find(hash);
    if (it != hash_positions.end()) {
        return it->second;
    }
    return {};
}

std::vector<SVCall> detectDeletionsFromBubbles(const std::vector<Bubble>& bubbles, const LCPGraph& graph, const LinearReference& reference) {
    std::vector<SVCall> sv_calls;
    
    for (const auto& bubble : bubbles) {
        // Use start and end nodes from bubble
        node_hash_t start_hash = bubble.start_node;
        node_hash_t end_hash = bubble.end_node;
        
        // Get positions for start node
        std::vector<size_t> start_indices = reference.align(start_hash);
        std::vector<size_t> end_indices = reference.align(end_hash);
        
        if (start_indices.empty() || end_indices.empty()) continue;
        
        // Create SV call from bubble's start and end nodes
        SVCall call;
        call.type = SVType::DELETION;
        call.node_hash = start_hash; // Use start node as primary hash
        
        // Collect all unique positions from start and end nodes
        std::set<std::pair<std::string, std::pair<uint64_t, uint64_t>>> seen;
        
        for (size_t idx : start_indices) {
            std::string chrom = reference.getChrom(idx);
            uint64_t start = reference.getStart(idx);
            uint64_t end = reference.getEnd(idx);
            auto key = std::make_pair(chrom, std::make_pair(start, end));
            if (seen.insert(key).second) {
                if (call.chrom.empty()) {
                    call.chrom = chrom;
                }
                call.genomic_starts.push_back(start);
                call.genomic_ends.push_back(end);
            }
        }
        
        for (size_t idx : end_indices) {
            std::string chrom = reference.getChrom(idx);
            uint64_t start = reference.getStart(idx);
            uint64_t end = reference.getEnd(idx);
            auto key = std::make_pair(chrom, std::make_pair(start, end));
            if (seen.insert(key).second) {
                if (call.chrom.empty()) {
                    call.chrom = chrom;
                }
                call.genomic_starts.push_back(start);
                call.genomic_ends.push_back(end);
            }
        }
        
        // Set count from graph nodes
        const auto& nodes = graph.getNodes();
        auto start_it = nodes.find(start_hash);
        if (start_it != nodes.end()) {
            call.count = start_it->second.count;
        }
        
        call.tumor_samples.push_back(1); // Default tumor sample
        
        if (!call.genomic_starts.empty()) {
            sv_calls.push_back(call);
        }
    }
    
    return sv_calls;
}

void exportTSV(const std::vector<SVCall>& sv_calls, const std::string& filename) {
    std::ofstream sv_file(filename);
    if (!sv_file) return;
    
    sv_file << "type\tchrom\tstart\tend\ttumor_samples\tcount\n";
    
    // Collect all output rows: (chrom, start, line_content)
    using RowData = std::tuple<std::string, uint64_t, std::string>;
    std::vector<RowData> rows;
    
    for (const auto& sv : sv_calls) {
        std::string type_str = (sv.type == SVType::DELETION) ? "DEL" : "UNK";
        
        // Build tumor samples string
        std::string samples_str;
        for (size_t i = 0; i < sv.tumor_samples.size(); ++i) {
            if (i > 0) samples_str += ",";
            samples_str += std::to_string(sv.tumor_samples[i]);
        }
        if (samples_str.empty()) samples_str = ".";
        
        std::string chrom_str = sv.chrom.empty() ? "." : sv.chrom;
        
        for (size_t i = 0; i < sv.genomic_starts.size(); ++i) {
            std::ostringstream oss;
            oss << type_str << "\t" 
                << chrom_str << "\t"
                << sv.genomic_starts[i] << "\t" 
                << sv.genomic_ends[i] << "\t"
                << samples_str << "\t"
                << sv.count;
            
            rows.emplace_back(chrom_str, sv.genomic_starts[i], oss.str());
        }
    }
    
    // Sort rows: first by chromosome, then by start position
    std::sort(rows.begin(), rows.end(), [](const RowData& a, const RowData& b) {
        const auto& chrom_a = std::get<0>(a);
        const auto& chrom_b = std::get<0>(b);
        const auto& start_a = std::get<1>(a);
        const auto& start_b = std::get<1>(b);
        
        int chrom_cmp = compareChromosomes(chrom_a, chrom_b);
        if (chrom_cmp != 0) return chrom_cmp < 0;
        return start_a < start_b;
    });
    
    // Write sorted rows
    for (const auto& row : rows) {
        sv_file << std::get<2>(row) << "\n";
    }
}

void exportBubblesTSV(const std::vector<Bubble>& bubbles, const std::string& output_dir, const std::string& prefix) {
    namespace fs = std::filesystem;
    fs::path tsv_dir = output_dir;
    fs::create_directories(tsv_dir);
    fs::path filename = tsv_dir / (prefix + "_bubbles.tsv");
    std::ofstream out(filename);
    if (!out) return;
    
    out << "type\tchrom\tstart\tend\tstart_node\tend_node\t"
        << "tumor_path_count\tnormal_path_count\t"
        << "tumor_paths_info\tnormal_paths_info\n";
    
    for (const auto& bubble : bubbles) {
        std::string type_str = (bubble.type == SVType::DELETION) ? "DELETION" : "INSERTION";
        
        // Format: path1_length:node_count:hashes;path2_length:node_count:hashes;...
        std::string tumor_paths_info;
        for (size_t i = 0; i < bubble.tumor_paths.size(); ++i) {
            if (i > 0) tumor_paths_info += ";";
            const auto& tpath = bubble.tumor_paths[i];
            tumor_paths_info += std::to_string(tpath.path_length) + ":" + std::to_string(tpath.node_count) + ":";
            for (size_t j = 0; j < tpath.internal_nodes.size(); ++j) {
                if (j > 0) tumor_paths_info += ",";
                tumor_paths_info += std::to_string(tpath.internal_nodes[j]);
            }
        }
        
        std::string normal_paths_info;
        for (size_t i = 0; i < bubble.normal_paths.size(); ++i) {
            if (i > 0) normal_paths_info += ";";
            const auto& npath = bubble.normal_paths[i];
            normal_paths_info += std::to_string(npath.path_length) + ":" + std::to_string(npath.node_count) + ":";
            for (size_t j = 0; j < npath.internal_nodes.size(); ++j) {
                if (j > 0) normal_paths_info += ",";
                normal_paths_info += std::to_string(npath.internal_nodes[j]);
            }
        }
        
        out << type_str << "\t" << bubble.chrom << "\t"
            << bubble.start_pos << "\t" << bubble.end_pos << "\t"
            << bubble.start_node << "\t" << bubble.end_node << "\t"
            << bubble.tumor_paths.size() << "\t" << bubble.normal_paths.size() << "\t"
            << tumor_paths_info << "\t" << normal_paths_info << "\n";
    }
}

void exportBubblesDetailedTSV(const std::vector<Bubble>& bubbles, const LCPGraph& graph, const std::string& output_dir, const std::string& prefix) {
    namespace fs = std::filesystem;
    fs::path tsv_dir = output_dir;
    fs::create_directories(tsv_dir);
    fs::path filename = tsv_dir / (prefix + "_bubbles_detailed.tsv");
    std::ofstream out(filename);
    if (!out) return;
    
    out << "bubble_id\ttype\tchrom\tstart_pos\tend_pos\tbubble_length\t"
        << "start_node\tend_node\t"
        << "tumor_path_count\tnormal_path_count\t"
        << "tumor_paths_detail\tnormal_paths_detail\t"
        << "start_node_positions\tend_node_positions\n";
    
    for (size_t idx = 0; idx < bubbles.size(); ++idx) {
        const auto& bubble = bubbles[idx];
        std::string type_str = (bubble.type == SVType::DELETION) ? "DELETION" : "INSERTION";
        
        // Format tumor paths: path_idx:length:node_count:chrom:start:end:hashes:positions;...
        std::string tumor_paths_detail;
        for (size_t i = 0; i < bubble.tumor_paths.size(); ++i) {
            if (i > 0) tumor_paths_detail += ";";
            const auto& tpath = bubble.tumor_paths[i];
            tumor_paths_detail += std::to_string(i) + ":" + std::to_string(tpath.path_length) + ":" 
                + std::to_string(tpath.node_count) + ":" + tpath.chrom + ":" 
                + std::to_string(tpath.start_pos) + ":" + std::to_string(tpath.end_pos) + ":";
            // Hashes
            for (size_t j = 0; j < tpath.internal_nodes.size(); ++j) {
                if (j > 0) tumor_paths_detail += ",";
                tumor_paths_detail += std::to_string(tpath.internal_nodes[j]);
            }
            tumor_paths_detail += ":";
            // Positions
            bool first_pos = true;
            for (node_hash_t h : tpath.internal_nodes) {
                const auto& positions = graph.getPositions(h);
                for (const auto& pos : positions) {
                    if (!first_pos) tumor_paths_detail += ",";
                    tumor_paths_detail += std::to_string(pos.first) + "-" + std::to_string(pos.second);
                    first_pos = false;
                }
            }
        }
        
        // Format normal paths similarly
        std::string normal_paths_detail;
        for (size_t i = 0; i < bubble.normal_paths.size(); ++i) {
            if (i > 0) normal_paths_detail += ";";
            const auto& npath = bubble.normal_paths[i];
            normal_paths_detail += std::to_string(i) + ":" + std::to_string(npath.path_length) + ":" 
                + std::to_string(npath.node_count) + ":" + npath.chrom + ":" 
                + std::to_string(npath.start_pos) + ":" + std::to_string(npath.end_pos) + ":";
            // Hashes
            for (size_t j = 0; j < npath.internal_nodes.size(); ++j) {
                if (j > 0) normal_paths_detail += ",";
                normal_paths_detail += std::to_string(npath.internal_nodes[j]);
            }
            normal_paths_detail += ":";
            // Positions
            bool first_pos = true;
            for (node_hash_t h : npath.internal_nodes) {
                const auto& positions = graph.getPositions(h);
                for (const auto& pos : positions) {
                    if (!first_pos) normal_paths_detail += ",";
                    normal_paths_detail += std::to_string(pos.first) + "-" + std::to_string(pos.second);
                    first_pos = false;
                }
            }
        }
        
        // Start node positions
        std::string start_positions;
        bool first_start = true;
        const auto& start_poss = graph.getPositions(bubble.start_node);
        for (const auto& pos : start_poss) {
            if (!first_start) start_positions += ";";
            start_positions += std::to_string(pos.first) + "-" + std::to_string(pos.second);
            first_start = false;
        }
        if (start_positions.empty()) start_positions = ".";
        
        // End node positions
        std::string end_positions;
        bool first_end = true;
        const auto& end_poss = graph.getPositions(bubble.end_node);
        for (const auto& pos : end_poss) {
            if (!first_end) end_positions += ";";
            end_positions += std::to_string(pos.first) + "-" + std::to_string(pos.second);
            first_end = false;
        }
        if (end_positions.empty()) end_positions = ".";
        
        uint64_t bubble_length = (bubble.end_pos > bubble.start_pos) ? (bubble.end_pos - bubble.start_pos) : 0;
        
        out << idx << "\t" << type_str << "\t" << bubble.chrom << "\t"
            << bubble.start_pos << "\t" << bubble.end_pos << "\t" << bubble_length << "\t"
            << bubble.start_node << "\t" << bubble.end_node << "\t"
            << bubble.tumor_paths.size() << "\t" << bubble.normal_paths.size() << "\t"
            << tumor_paths_detail << "\t" << normal_paths_detail << "\t"
            << start_positions << "\t" << end_positions << "\n";
    }
}

void exportBubblePathsTSV(const std::vector<BubblePathReport>& paths, const std::string& output_dir, const std::string& prefix) {
    namespace fs = std::filesystem;
    fs::path tsv_dir = output_dir;
    fs::create_directories(tsv_dir);
    fs::path filename = tsv_dir / (prefix + "_bubble_paths.tsv");
    std::ofstream out(filename);
    if (!out) return;

    out << "path_type\tchrom\tstart\tend\tpath_length\tnode_count\tnode_hashes\n";

    for (const auto& rep : paths) {
        std::string hashes;
        for (size_t i = 0; i < rep.internal_nodes.size(); ++i) {
            if (i > 0) hashes += ",";
            hashes += std::to_string(rep.internal_nodes[i]);
        }

        out << rep.path_type << "\t"
            << rep.chrom << "\t"
            << rep.start_pos << "\t"
            << rep.end_pos << "\t"
            << rep.path_length << "\t"
            << rep.node_count << "\t"
            << hashes << "\n";
    }
}

// Bubble detection - find alternative paths with different colors
namespace {
    uint64_t pathSpan(const LCPGraph& graph, const std::vector<node_hash_t>& nodes) {
        if (nodes.empty()) return 0;
        uint64_t mn = std::numeric_limits<uint64_t>::max();
        uint64_t mx = 0;
        bool has = false;
        for (node_hash_t h : nodes) {
            const auto& positions = graph.getPositions(h);
            if (positions.empty()) continue;
            for (const auto& pos : positions) {
                mn = std::min<uint64_t>(mn, pos.first);
                mx = std::max<uint64_t>(mx, pos.second);
                has = true;
            }
        }
        if (!has) return nodes.size(); // fallback
        return (mx > mn) ? (mx - mn) : 0;
    }

    struct AlignmentInfo {
        std::string chrom;
        uint64_t start = 0;
        uint64_t end = 0;
    };

    std::vector<AlignmentInfo> collectAlignments(node_hash_t hash, const LinearReference& reference) {
        std::vector<AlignmentInfo> res;
        auto idxs = reference.align(hash);
        res.reserve(idxs.size());
        for (size_t idx : idxs) {
            res.push_back(AlignmentInfo{
                reference.getChrom(idx),
                reference.getStart(idx),
                reference.getEnd(idx)
            });
        }
        // deterministic order: chrom then start
        std::sort(res.begin(), res.end(), [](const AlignmentInfo& a, const AlignmentInfo& b) {
            int cmp = compareChromosomes(a.chrom, b.chrom);
            if (cmp != 0) return cmp < 0;
            return a.start < b.start;
        });
        return res;
    }

    // Pick one alignment per chromosome: earliest start and latest end on that chrom
    std::vector<AlignmentInfo> compressAlignmentsByChrom(const std::vector<AlignmentInfo>& aligns) {
        std::map<std::string, std::pair<uint64_t,uint64_t>> by_chrom;
        for (const auto& a : aligns) {
            auto& p = by_chrom[a.chrom];
            if (p.first == 0 && p.second == 0) {
                p.first = a.start; p.second = a.end;
            } else {
                p.first = std::min<uint64_t>(p.first, a.start);
                p.second = std::max<uint64_t>(p.second, a.end);
            }
        }
        std::vector<AlignmentInfo> res;
        res.reserve(by_chrom.size());
        for (const auto& kv : by_chrom) {
            res.push_back(AlignmentInfo{kv.first, kv.second.first, kv.second.second});
        }
        return res;
    }
}

std::vector<Bubble> detectBubbles(
    const LCPGraph& graph,
    const LinearReference& reference,
    std::vector<BubblePathReport>* tumor_paths_out,
    std::vector<BubblePathReport>* normal_paths_out) {
    sv_logger.debug("Starting bubble detection...");
    std::vector<Bubble> bubbles;
    
    const auto& nodes = graph.getNodes();
    
    // Debug: count branch points and T+N nodes
    size_t branch_points = 0;
    size_t tn_nodes = 0;
    size_t tn_branch_points = 0;
    
    for (const auto& [hash, node] : nodes) {
        const auto& next_set = graph.getNext(hash);
        if (next_set.size() >= 2) {
            branch_points++;
        }
        
        bool has_normal = node.hasColor(0);
        bool has_tumor = false;
        for (int i = 1; i < 32; ++i) {
            if (node.hasColor(i)) {
                has_tumor = true;
                break;
            }
        }
        
        if (has_normal && has_tumor) {
            tn_nodes++;
            if (next_set.size() >= 2) {
                tn_branch_points++;
            }
        }
    }
    
    sv_logger.debug("Graph stats: " + std::to_string(nodes.size()) + " nodes, " + 
                    std::to_string(branch_points) + " branch points, " +
                    std::to_string(tn_nodes) + " T+N nodes, " +
                    std::to_string(tn_branch_points) + " T+N branch points");
    
    // Find branch points (nodes with multiple outgoing edges)
    for (const auto& [hash, node] : nodes) {
        const auto& next_set = graph.getNext(hash);
        if (next_set.size() < 2) continue;  // Not a branch point
        
        // Check if this node has both T and N colors (T+N)
        bool has_normal = node.hasColor(0);
        bool has_tumor = false;
        for (int i = 1; i < 32; ++i) {
            if (node.hasColor(i)) {
                has_tumor = true;
                break;
            }
        }
        
        if (!has_normal || !has_tumor) continue;  // Start node must be T+N
        
        // Trace all paths from this branch point
        std::vector<std::vector<node_hash_t>> all_paths;
        std::set<node_hash_t> visited;
        
        // DFS to find all paths from branch point
        std::function<void(node_hash_t, std::vector<node_hash_t>)> dfs = 
            [&](node_hash_t current, std::vector<node_hash_t> path) {
                if (visited.count(current)) return;
                visited.insert(current);
                
                path.push_back(current);
                
                const auto& next = graph.getNext(current);
                if (next.empty()) {
                    // End of path
                    all_paths.push_back(path);
                } else {
                    for (node_hash_t next_node : next) {
                        auto next_it = nodes.find(next_node);
                        if (next_it != nodes.end()) {
                            // Check if next node is T+N (convergence point)
                            bool next_has_normal = next_it->second.hasColor(0);
                            bool next_has_tumor = false;
                            for (int i = 1; i < 32; ++i) {
                                if (next_it->second.hasColor(i)) {
                                    next_has_tumor = true;
                                    break;
                                }
                            }
                            
                            if (next_has_normal && next_has_tumor) {
                                // Convergence point reached
                                path.push_back(next_node);
                                all_paths.push_back(path);
                            } else {
                                dfs(next_node, path);
                            }
                        }
                    }
                }
                
                visited.erase(current);
            };
        
        for (node_hash_t next_node : next_set) {
            dfs(next_node, {hash});
        }
        
        // Group paths by color composition
        std::vector<std::vector<node_hash_t>> tumor_paths;
        std::vector<std::vector<node_hash_t>> normal_paths;
        
        for (const auto& path : all_paths) {
            if (path.size() < 2) continue;
            
            // Check color composition of path (excluding start and end T+N nodes)
            bool is_tumor_only = true;
            bool is_normal_only = true;
            
            for (size_t i = 1; i + 1 < path.size(); ++i) {
                auto it = nodes.find(path[i]);
                if (it == nodes.end()) continue;
                
                bool has_n = it->second.hasColor(0);
                bool has_t = false;
                for (int j = 1; j < 32; ++j) {
                    if (it->second.hasColor(j)) {
                        has_t = true;
                        break;
                    }
                }
                
                if (has_n) is_tumor_only = false;
                if (has_t) is_normal_only = false;
            }
            
            if (is_tumor_only && !is_normal_only) {
                tumor_paths.push_back(path);
            } else if (is_normal_only && !is_tumor_only) {
                normal_paths.push_back(path);
            }
        }
        
        // Emit per-path reports (optional)
        auto emit_report = [&](const std::vector<node_hash_t>& path, const std::string& label, std::vector<BubblePathReport>* out_vec) {
            if (!out_vec) return;
            if (path.size() < 2) return;

            // internal nodes exclude start/end T+N
            std::vector<node_hash_t> internal;
            for (size_t i = 1; i + 1 < path.size(); ++i) {
                internal.push_back(path[i]);
            }
            uint64_t path_span = pathSpan(graph, internal);

            auto start_aligns_raw = collectAlignments(path.front(), reference);
            auto end_aligns_raw   = collectAlignments(path.back(), reference);
            auto start_aligns = compressAlignmentsByChrom(start_aligns_raw);
            auto end_aligns   = compressAlignmentsByChrom(end_aligns_raw);

            if (start_aligns.empty() || end_aligns.empty()) {
                BubblePathReport rep;
                rep.path_type = label;
                rep.start_node = path.front();
                rep.end_node = path.back();
                rep.internal_nodes = internal;
                rep.node_count = internal.size();
                rep.path_length = path_span;
                out_vec->push_back(std::move(rep));
                return;
            }

            // Emit for all same-chrom combinations
            for (const auto& sa : start_aligns) {
                for (const auto& ea : end_aligns) {
                    if (sa.chrom != ea.chrom) continue;
                    BubblePathReport rep;
                    rep.path_type = label;
                    rep.start_node = path.front();
                    rep.end_node = path.back();
                    rep.internal_nodes = internal;
                    rep.node_count = internal.size();
                    rep.path_length = path_span;

                    // span from anchors + internal nodes
                    uint64_t mn = sa.start, mx = ea.end;
                    if (!internal.empty()) {
                        uint64_t inner_mn = std::numeric_limits<uint64_t>::max();
                        uint64_t inner_mx = 0;
                        bool has_inner = false;
                        for (node_hash_t h : internal) {
                            for (const auto& pos : graph.getPositions(h)) {
                                inner_mn = std::min<uint64_t>(inner_mn, pos.first);
                                inner_mx = std::max<uint64_t>(inner_mx, pos.second);
                                has_inner = true;
                            }
                        }
                        if (has_inner) {
                            mn = std::min<uint64_t>(mn, inner_mn);
                            mx = std::max<uint64_t>(mx, inner_mx);
                        }
                    }
                    rep.chrom = sa.chrom;
                    rep.start_pos = mn;
                    rep.end_pos = mx;
                    out_vec->push_back(std::move(rep));
                }
            }
        };

        for (const auto& p : tumor_paths) emit_report(p, "TUMOR_ONLY", tumor_paths_out);
        for (const auto& p : normal_paths) emit_report(p, "NORMAL_ONLY", normal_paths_out);

        // Create single bubble per branch point with all T-only and N-only paths
        // Group paths by end node (must converge to same end)
        std::map<node_hash_t, std::pair<std::vector<std::vector<node_hash_t>>, std::vector<std::vector<node_hash_t>>>> paths_by_end;
        
        for (const auto& tumor_path : tumor_paths) {
            paths_by_end[tumor_path.back()].first.push_back(tumor_path);
        }
        for (const auto& normal_path : normal_paths) {
            paths_by_end[normal_path.back()].second.push_back(normal_path);
        }
        
        // Create one bubble per end node convergence
        for (const auto& [end_node, path_groups] : paths_by_end) {
            const auto& tumor_group = path_groups.first;
            const auto& normal_group = path_groups.second;
            
            if (tumor_group.empty() || normal_group.empty()) continue; // Need both types
            
            Bubble bubble;
            bubble.start_node = hash;
            bubble.end_node = end_node;
            
            // Helper to create PathInfo from a path
            auto createPathInfo = [&](const std::vector<node_hash_t>& path) -> PathInfo {
                PathInfo info;
                // Extract internal nodes (excluding start and end T+N nodes)
                for (size_t i = 1; i + 1 < path.size(); ++i) {
                    info.internal_nodes.push_back(path[i]);
                }
                info.path_length = pathSpan(graph, info.internal_nodes);
                info.node_count = info.internal_nodes.size();
                
                // Get positions from reference alignment
                auto start_aligns_raw = collectAlignments(path.front(), reference);
                auto end_aligns_raw = collectAlignments(path.back(), reference);
                auto start_aligns = compressAlignmentsByChrom(start_aligns_raw);
                auto end_aligns = compressAlignmentsByChrom(end_aligns_raw);
                
                if (!start_aligns.empty() && !end_aligns.empty()) {
                    // Use first same-chrom combination for this path
                    for (const auto& sa : start_aligns) {
                        for (const auto& ea : end_aligns) {
                            if (sa.chrom == ea.chrom) {
                                uint64_t mn = sa.start, mx = ea.end;
                                // Include internal node spans
                                if (!info.internal_nodes.empty()) {
                                    uint64_t inner_mn = std::numeric_limits<uint64_t>::max();
                                    uint64_t inner_mx = 0;
                                    bool has_inner = false;
                                    for (node_hash_t h : info.internal_nodes) {
                                        for (const auto& pos : graph.getPositions(h)) {
                                            inner_mn = std::min<uint64_t>(inner_mn, pos.first);
                                            inner_mx = std::max<uint64_t>(inner_mx, pos.second);
                                            has_inner = true;
                                        }
                                    }
                                    if (has_inner) {
                                        mn = std::min<uint64_t>(mn, inner_mn);
                                        mx = std::max<uint64_t>(mx, inner_mx);
                                    }
                                }
                                info.chrom = sa.chrom;
                                info.start_pos = mn;
                                info.end_pos = mx;
                                return info;
                            }
                        }
                    }
                }
                return info; // Return with empty chrom if no alignment
            };
            
            // Create PathInfo for all tumor paths
            for (const auto& tumor_path : tumor_group) {
                bubble.tumor_paths.push_back(createPathInfo(tumor_path));
            }
            
            // Create PathInfo for all normal paths
            for (const auto& normal_path : normal_group) {
                bubble.normal_paths.push_back(createPathInfo(normal_path));
            }
            
            // Calculate overall bubble span (min start, max end across all paths)
            uint64_t bubble_mn = std::numeric_limits<uint64_t>::max();
            uint64_t bubble_mx = 0;
            std::string bubble_chrom;
            bool has_pos = false;
            
            for (const auto& tpath : bubble.tumor_paths) {
                if (!tpath.chrom.empty() && tpath.end_pos > tpath.start_pos) {
                    if (bubble_chrom.empty()) bubble_chrom = tpath.chrom;
                    if (tpath.chrom == bubble_chrom) {
                        bubble_mn = std::min<uint64_t>(bubble_mn, tpath.start_pos);
                        bubble_mx = std::max<uint64_t>(bubble_mx, tpath.end_pos);
                        has_pos = true;
                    }
                }
            }
            for (const auto& npath : bubble.normal_paths) {
                if (!npath.chrom.empty() && npath.end_pos > npath.start_pos) {
                    if (bubble_chrom.empty()) bubble_chrom = npath.chrom;
                    if (npath.chrom == bubble_chrom) {
                        bubble_mn = std::min<uint64_t>(bubble_mn, npath.start_pos);
                        bubble_mx = std::max<uint64_t>(bubble_mx, npath.end_pos);
                        has_pos = true;
                    }
                }
            }
            
            if (has_pos) {
                bubble.chrom = bubble_chrom;
                bubble.start_pos = bubble_mn;
                bubble.end_pos = bubble_mx;
                
                // Determine SV type: if any normal path is longer than any tumor path, it's a deletion
                uint64_t max_normal_len = 0, max_tumor_len = 0;
                for (const auto& npath : bubble.normal_paths) {
                    max_normal_len = std::max(max_normal_len, npath.path_length);
                }
                for (const auto& tpath : bubble.tumor_paths) {
                    max_tumor_len = std::max(max_tumor_len, tpath.path_length);
                }
                bubble.type = (max_normal_len > max_tumor_len) ? SVType::DELETION : SVType::DELETION;
                
                bubbles.push_back(std::move(bubble));
            }
        }
    }
    
    sv_logger.debug("Found " + std::to_string(bubbles.size()) + " bubbles");
    return bubbles;
}
