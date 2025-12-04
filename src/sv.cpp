#include "sv.h"
#include <algorithm>
#include <set>

bool isInNormal(const CompressedNode& node) {
    return node.colors[0] > 0;
}

bool isInTumor(const CompressedNode& node) {
    for (size_t i = 1; i < MAX_COLORS; ++i) {
        if (node.colors[i] > 0) {
            return true;
        }
    }
    return false;
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

std::vector<size_t> LinearReference::align(const std::vector<ublock>& path) const {
    std::vector<size_t> positions;
    
    if (path.empty() || hashes.empty()) {
        return positions;
    }
    
    ublock first_hash = path[0];
    auto it = hash_positions.find(first_hash);
    
    if (it == hash_positions.end()) {
        for (size_t i = 0; i < path.size(); ++i) {
            auto pos_it = hash_positions.find(path[i]);
            if (pos_it != hash_positions.end()) {
                for (size_t pos : pos_it->second) {
                    if (pos >= i) {
                        positions.push_back(pos - i);
                    }
                }
            }
        }
        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
        return positions;
    }
    
    for (size_t start_pos : it->second) {
        size_t match_count = 0;
        
        for (size_t i = 0; i < path.size() && start_pos + i < hashes.size(); ++i) {
            if (hashes[start_pos + i] == path[i]) {
                match_count++;
            } else {
                break;
            }
        }
        
        if (match_count > 0) {
            positions.push_back(start_pos);
        }
    }
    
    return positions;
}

// Check if an edge exists only in tumor samples (not in normal)
// This is detected by looking at edges between nodes with different color patterns
bool isTumorOnlyEdge(const CompressedNode& from_node, const CompressedNode& to_node) {
    // Edge exists in tumor if both nodes exist in tumor
    bool in_tumor = false;
    for (size_t i = 1; i < MAX_COLORS; ++i) {
        if (from_node.colors[i] > 0 && to_node.colors[i] > 0) {
            in_tumor = true;
            break;
        }
    }
    // Edge is tumor-only if not present in normal
    bool in_normal = (from_node.colors[0] > 0 && to_node.colors[0] > 0);
    return in_tumor && !in_normal;
}

std::vector<SVCall> detectPresenceAbsenceSVs(const LCPGraph& graph, const LinearReference& reference, int num_tumor_samples) {
    std::vector<SVCall> sv_calls;
    
    if (!graph.isCompressed()) {
        return sv_calls;
    }
    
    const auto& compressed_nodes = graph.getCompressedNodes();
    
    // Build reverse edges map
    std::unordered_map<ublock, std::vector<ublock>> reverse_edges;
    for (const auto& [hash, node] : compressed_nodes) {
        for (ublock next_hash : node.next) {
            reverse_edges[next_hash].push_back(hash);
        }
    }
    
    // Track which nodes are part of back-edges (for cycle detection)
    std::unordered_set<ublock> back_edge_targets;
    
    for (const auto& [hash, node] : compressed_nodes) {
        bool normal_present = isInNormal(node);
        unsigned int normal_count = node.colors[0];
        
        // DELETION: Present in normal, absent in specific tumor(s)
        if (normal_present) {
            for (int tumor_id = 1; tumor_id <= num_tumor_samples; ++tumor_id) {
                if (node.colors[tumor_id] == 0) {
                    SVCall call;
                    call.type = SVType::DELETION;
                    call.node_hash = hash;
                    call.path = node.path;
                    call.tumor_samples.push_back(tumor_id);
                    call.normal_count = normal_count;
                    call.tumor_count = 0;
                    call.reference_indices = reference.align(node.path);
                    
                    std::set<std::pair<uint64_t, uint64_t>> seen;
                    for (size_t idx : call.reference_indices) {
                        if (call.chrom.empty()) {
                            call.chrom = reference.getChrom(idx);
                        }
                        auto pos = std::make_pair(reference.getStart(idx), reference.getEnd(idx));
                        if (seen.insert(pos).second) {
                            call.genomic_starts.push_back(pos.first);
                            call.genomic_ends.push_back(pos.second);
                        }
                    }
                    sv_calls.push_back(call);
                }
            }
        }
        
        // DUPLICATION (count-based): Present in normal AND tumor has higher count
        if (normal_present && normal_count > 0) {
            for (int tumor_id = 1; tumor_id <= num_tumor_samples; ++tumor_id) {
                unsigned int tumor_count = node.colors[tumor_id];
                if (tumor_count > normal_count) {
                    SVCall call;
                    call.type = SVType::DUPLICATION;
                    call.node_hash = hash;
                    call.path = node.path;
                    call.tumor_samples.push_back(tumor_id);
                    call.normal_count = normal_count;
                    call.tumor_count = tumor_count;
                    call.reference_indices = reference.align(node.path);
                    
                    std::set<std::pair<uint64_t, uint64_t>> seen;
                    for (size_t idx : call.reference_indices) {
                        if (call.chrom.empty()) {
                            call.chrom = reference.getChrom(idx);
                        }
                        auto pos = std::make_pair(reference.getStart(idx), reference.getEnd(idx));
                        if (seen.insert(pos).second) {
                            call.genomic_starts.push_back(pos.first);
                            call.genomic_ends.push_back(pos.second);
                        }
                    }
                    sv_calls.push_back(call);
                }
            }
        }
        
        // DUPLICATION (cycle-based): Look for tumor-only back-edges
        // A tandem duplication ABC->ABCABC creates edge C->A only in tumor
        if (normal_present) {
            for (ublock next_hash : node.next) {
                auto next_it = compressed_nodes.find(next_hash);
                if (next_it == compressed_nodes.end()) continue;
                
                const auto& next_node = next_it->second;
                
                // Check if this edge only exists in tumor (potential back-edge from duplication)
                // We look for cases where normal has A->B but tumor also has A->B->...->A
                if (isInNormal(next_node)) {
                    // Both nodes exist in normal, check if tumor has extra copies
                    // This would show as tumor_count > normal_count (handled above)
                } else {
                    // next_node only in tumor - this edge is tumor-only
                    // Could be part of a duplication if it eventually loops back
                    // Track for potential cycle detection
                    back_edge_targets.insert(next_hash);
                }
            }
        }
        
        // INSERTION: Absent in normal, present in specific tumor(s)
        if (!normal_present) {
            for (int tumor_id = 1; tumor_id <= num_tumor_samples; ++tumor_id) {
                if (node.colors[tumor_id] > 0) {
                    SVCall call;
                    call.type = SVType::INSERTION;
                    call.node_hash = hash;
                    call.path = node.path;
                    call.tumor_samples.push_back(tumor_id);
                    call.normal_count = 0;
                    call.tumor_count = node.colors[tumor_id];
                    
                    // Find flanking nodes in normal
                    auto rev_it = reverse_edges.find(hash);
                    if (rev_it != reverse_edges.end()) {
                        for (ublock prev_hash : rev_it->second) {
                            auto prev_it = compressed_nodes.find(prev_hash);
                            if (prev_it != compressed_nodes.end() && isInNormal(prev_it->second)) {
                                for (ublock h : prev_it->second.path) {
                                    call.flanking_before.push_back(h);
                                }
                                break;
                            }
                        }
                    }
                    
                    for (ublock next_hash : node.next) {
                        auto next_it = compressed_nodes.find(next_hash);
                        if (next_it != compressed_nodes.end() && isInNormal(next_it->second)) {
                            for (ublock h : next_it->second.path) {
                                call.flanking_after.push_back(h);
                            }
                            break;
                        }
                    }

                    if (!call.flanking_before.empty()) {
                        call.reference_indices = reference.align(call.flanking_before);
                        for (auto& idx : call.reference_indices) {
                            idx += call.flanking_before.size();
                        }
                    } else if (!call.flanking_after.empty()) {
                        call.reference_indices = reference.align(call.flanking_after);
                    }
                    
                    std::set<std::pair<uint64_t, uint64_t>> seen;
                    for (size_t idx : call.reference_indices) {
                        if (call.chrom.empty()) {
                            call.chrom = reference.getChrom(idx);
                        }
                        auto pos = std::make_pair(reference.getStart(idx), reference.getEnd(idx));
                        if (seen.insert(pos).second) {
                            call.genomic_starts.push_back(pos.first);
                            call.genomic_ends.push_back(pos.second);
                        }
                    }
                    sv_calls.push_back(call);
                }
            }
        }
    }
    
    return sv_calls;
}
