#include "sv.h"
#include "utils.h"
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

void LinearReference::addCore(ublock hash, const std::string& chrom, uint64_t start, uint64_t end) {
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

std::vector<size_t> LinearReference::align(ublock hash) const {
    auto it = hash_positions.find(hash);
    if (it != hash_positions.end()) {
        return it->second;
    }
    return {};
}

std::vector<SVCall> detectDeletions(const LCPGraph& graph, const LinearReference& reference, int num_tumor_samples) {
    std::vector<SVCall> sv_calls;
    const auto& nodes = graph.getNodes();
    
    for (const auto& [hash, node] : nodes) {
        // DELETION: Present in normal (color 0), absent in tumor
        if (!node.hasColor(0)) continue;  // Must be in normal
        
        // Check each tumor sample
        for (int tumor_id = 1; tumor_id <= num_tumor_samples; ++tumor_id) {
            if (!node.hasColor(tumor_id)) {
                // Found a deletion!
                SVCall call;
                call.type = SVType::DELETION;
                call.node_hash = hash;
                call.count = node.count;
                call.tumor_samples.push_back(tumor_id);
                
                // Align to reference
                std::vector<size_t> indices = reference.align(hash);
                
                std::set<std::pair<uint64_t, uint64_t>> seen;
                for (size_t idx : indices) {
                    if (call.chrom.empty()) {
                        call.chrom = reference.getChrom(idx);
                    }
                    auto pos = std::make_pair(reference.getStart(idx), reference.getEnd(idx));
                    if (seen.insert(pos).second) {
                        call.genomic_starts.push_back(pos.first);
                        call.genomic_ends.push_back(pos.second);
                    }
                }
                
                if (!call.genomic_starts.empty()) {
                    sv_calls.push_back(call);
                }
            }
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
    
    out << "type\tchrom\tstart\tend\ttumor_path_length\tnormal_path_length\t"
        << "tumor_node_count\tnormal_node_count\ttumor_node_hashes\tnormal_node_hashes\n";
    
    for (const auto& bubble : bubbles) {
        std::string type_str = (bubble.type == SVType::DELETION) ? "DELETION" : "INSERTION";
        
        std::string tumor_hashes;
        for (size_t i = 0; i < bubble.tumor_path.size(); ++i) {
            if (i > 0) tumor_hashes += ",";
            tumor_hashes += std::to_string(bubble.tumor_path[i]);
        }
        
        std::string normal_hashes;
        for (size_t i = 0; i < bubble.normal_path.size(); ++i) {
            if (i > 0) normal_hashes += ",";
            normal_hashes += std::to_string(bubble.normal_path[i]);
        }
        
        out << type_str << "\t" << bubble.chrom << "\t"
            << bubble.start_pos << "\t" << bubble.end_pos << "\t"
            << bubble.tumor_path_length << "\t" << bubble.normal_path_length << "\t"
            << bubble.tumor_path.size() << "\t" << bubble.normal_path.size() << "\t"
            << tumor_hashes << "\t" << normal_hashes << "\n";
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
        << "start_node\tedge_node\t"
        << "tumor_path_length\tnormal_path_length\t"
        << "tumor_node_count\tnormal_node_count\t"
        << "tumor_node_hashes\tnormal_node_hashes\t"
        << "tumor_node_positions\tnormal_node_positions\t"
        << "start_node_positions\tend_node_positions\n";
    
    for (size_t idx = 0; idx < bubbles.size(); ++idx) {
        const auto& bubble = bubbles[idx];
        std::string type_str = (bubble.type == SVType::DELETION) ? "DELETION" : "INSERTION";
        
        // Tumor path hashes
        std::string tumor_hashes;
        for (size_t i = 0; i < bubble.tumor_path.size(); ++i) {
            if (i > 0) tumor_hashes += ",";
            tumor_hashes += std::to_string(bubble.tumor_path[i]);
        }
        
        // Normal path hashes
        std::string normal_hashes;
        for (size_t i = 0; i < bubble.normal_path.size(); ++i) {
            if (i > 0) normal_hashes += ",";
            normal_hashes += std::to_string(bubble.normal_path[i]);
        }
        
        // Tumor path positions (all positions for all nodes)
        std::string tumor_positions;
        bool first_tumor = true;
        for (node_hash_t h : bubble.tumor_path) {
            const auto& positions = graph.getPositions(h);
            for (const auto& pos : positions) {
                if (!first_tumor) tumor_positions += ";";
                tumor_positions += std::to_string(pos.first) + "-" + std::to_string(pos.second);
                first_tumor = false;
            }
        }
        if (tumor_positions.empty()) tumor_positions = ".";
        
        // Normal path positions
        std::string normal_positions;
        bool first_normal = true;
        for (node_hash_t h : bubble.normal_path) {
            const auto& positions = graph.getPositions(h);
            for (const auto& pos : positions) {
                if (!first_normal) normal_positions += ";";
                normal_positions += std::to_string(pos.first) + "-" + std::to_string(pos.second);
                first_normal = false;
            }
        }
        if (normal_positions.empty()) normal_positions = ".";
        
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
            << bubble.tumor_path_length << "\t" << bubble.normal_path_length << "\t"
            << bubble.tumor_path.size() << "\t" << bubble.normal_path.size() << "\t"
            << tumor_hashes << "\t" << normal_hashes << "\t"
            << tumor_positions << "\t" << normal_positions << "\t"
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

    std::vector<AlignmentInfo> collectAlignments(ublock hash, const LinearReference& reference) {
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

        // Create bubbles: pair each tumor path with each normal path
        for (const auto& tumor_path : tumor_paths) {
            for (const auto& normal_path : normal_paths) {
                if (tumor_path.back() != normal_path.back()) continue;  // Must converge to same end node
                
                Bubble bubble;
                bubble.start_node = hash;
                bubble.end_node = tumor_path.back();
                
                // Extract internal nodes (excluding start and end T+N nodes)
                for (size_t i = 1; i + 1 < tumor_path.size(); ++i) {
                    bubble.tumor_path.push_back(tumor_path[i]);
                }
                for (size_t i = 1; i + 1 < normal_path.size(); ++i) {
                    bubble.normal_path.push_back(normal_path[i]);
                }
                
                // Calculate path spans (min/max) instead of sum of occurrences
                bubble.tumor_path_length = pathSpan(graph, bubble.tumor_path);
                bubble.normal_path_length = pathSpan(graph, bubble.normal_path);
                
                // Get positions from reference alignment (all combinations, same chrom)
                auto start_aligns_raw = collectAlignments(bubble.start_node, reference);
                auto end_aligns_raw   = collectAlignments(bubble.end_node, reference);
                auto start_aligns = compressAlignmentsByChrom(start_aligns_raw);
                auto end_aligns   = compressAlignmentsByChrom(end_aligns_raw);

                if (start_aligns.empty() || end_aligns.empty()) {
                    continue; // cannot place on reference
                }

                for (const auto& sa : start_aligns) {
                    for (const auto& ea : end_aligns) {
                        if (sa.chrom != ea.chrom) continue;

                        Bubble b = bubble; // copy base
                        b.chrom = sa.chrom;

                        uint64_t mn = sa.start, mx = ea.end;
                        // include internal node spans
                        if (!b.tumor_path.empty() || !b.normal_path.empty()) {
                            // aggregate all internal nodes of both paths for span
                            uint64_t inner_mn = std::numeric_limits<uint64_t>::max();
                            uint64_t inner_mx = 0;
                            bool has_inner = false;
                            for (node_hash_t h : b.tumor_path) {
                                for (const auto& pos : graph.getPositions(h)) {
                                    inner_mn = std::min<uint64_t>(inner_mn, pos.first);
                                    inner_mx = std::max<uint64_t>(inner_mx, pos.second);
                                    has_inner = true;
                                }
                            }
                            for (node_hash_t h : b.normal_path) {
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
                        b.start_pos = mn;
                        b.end_pos = mx;

                        // Determine SV type (still only DELETION for now)
                        b.type = (b.normal_path_length > b.tumor_path_length) ? 
                            SVType::DELETION : SVType::DELETION;

                        bubbles.push_back(std::move(b));
                    }
                }
            }
        }
    }
    
    sv_logger.debug("Found " + std::to_string(bubbles.size()) + " bubbles");
    return bubbles;
}
