#include "sv.h"
#include "utils.h"
#include <algorithm>
#include <set>
#include <fstream>

bool isInNormal(const CompressedNode& node) {
    return node.colors[0] > 0;
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
                    if (pos >= i) positions.push_back(pos - i);
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

// std::vector<SVCall> detectDeletions(const LCPGraph& graph, const LinearReference& reference, int num_tumor_samples) {
//     std::vector<SVCall> sv_calls;
    
//     if (!graph.isCompressed()) {
//         return sv_calls;
//     }
    
//     const auto& compressed_nodes = graph.getCompressedNodes();
    
//     for (const auto& [hash, node] : compressed_nodes) {
//         bool normal_present = isInNormal(node);
//         unsigned int normal_count = node.colors[0];
        
//         // DELETION: Present in normal, absent in specific tumor(s)
//         if (normal_present) {
//             for (int tumor_id = 1; tumor_id <= num_tumor_samples; ++tumor_id) {
//                 if (node.colors[tumor_id] == 0) {
//                     SVCall call;
//                     call.type = SVType::DELETION;
//                     call.node_hash = hash;
//                     call.path = node.path;
//                     call.tumor_samples.push_back(tumor_id);
//                     call.normal_count = normal_count;
//                     call.tumor_count = 0;
//                     call.reference_indices = reference.align(node.path);
                    
//                     std::set<std::pair<uint64_t, uint64_t>> seen;
//                     for (size_t idx : call.reference_indices) {
//                         if (call.chrom.empty()) {
//                             call.chrom = reference.getChrom(idx);
//                         }
//                         auto pos = std::make_pair(reference.getStart(idx), reference.getEnd(idx));
//                         if (seen.insert(pos).second) {
//                             call.genomic_starts.push_back(pos.first);
//                             call.genomic_ends.push_back(pos.second);
//                         }
//                     }
//                     sv_calls.push_back(call);
//                 }
//             }
//         }
//     }
    
//     return sv_calls;
// }

// Lightweight logger for debugging SV detection.
static Logger sv_logger("sv_logger");

std::vector<SVCall> detectDeletions(const LCPGraph& graph, const LinearReference& reference, int num_tumor_samples) {
    std::vector<SVCall> sv_calls;
    
    if (!graph.isCompressed()) {
        return sv_calls;
    }
    
    // Use original nodes instead of compressed nodes to avoid color summing issues
    const auto& node_dict = graph.getNodeDict();
    
    // Debug counters to understand why calls may be missing.
    std::vector<size_t> candidate_per_tumor(num_tumor_samples + 1, 0); // index by tumor_id
    std::vector<size_t> with_align_per_tumor(num_tumor_samples + 1, 0);
    std::vector<size_t> no_align_per_tumor(num_tumor_samples + 1, 0);
    
    for (const auto& [hash, data] : node_dict) {
        bool normal_present = (data.colors[0] > 0);
        unsigned int normal_count = data.colors[0];
        
        // DELETION: Present in normal, absent in specific tumor(s)
        if (normal_present) {
            for (int tumor_id = 1; tumor_id <= num_tumor_samples; ++tumor_id) {
                if (data.colors[tumor_id] == 0) {
                    candidate_per_tumor[tumor_id]++;
                    SVCall call;
                    call.type = SVType::DELETION;
                    call.node_hash = hash;
                    // For original nodes, path is just the single node
                    call.path = {hash};
                    call.tumor_samples.push_back(tumor_id);
                    call.normal_count = normal_count;
                    call.tumor_count = 0;
                    call.reference_indices = reference.align(call.path);

                    if (call.reference_indices.empty()) {
                        no_align_per_tumor[tumor_id]++;
                    } else {
                        with_align_per_tumor[tumor_id]++;
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

    // Emit a concise debug summary to help diagnose missing calls.
    for (int tumor_id = 1; tumor_id <= num_tumor_samples; ++tumor_id) {
        sv_logger.log_impl(
            "tumor_id=" + std::to_string(tumor_id) +
            " candidates=" + std::to_string(candidate_per_tumor[tumor_id]) +
            " align_ok=" + std::to_string(with_align_per_tumor[tumor_id]) +
            " align_empty=" + std::to_string(no_align_per_tumor[tumor_id]),
            __FILE__, __LINE__);
    }
    
    return sv_calls;
}

void exportTSV(const std::vector<SVCall>& sv_calls, const std::string& filename) {
    std::ofstream sv_file(filename);
    if (!sv_file) return;
    
    sv_file << "type\tchrom\tstart\tend\ttumor_samples\tnormal_count\ttumor_count\n";
    for (const auto& sv : sv_calls) {
        std::string type_str;
        switch (sv.type) {
            case SVType::DELETION: type_str = "DEL"; break;
        }
        
        // Build tumor samples string
        std::string samples_str;
        for (size_t i = 0; i < sv.tumor_samples.size(); ++i) {
            if (i > 0) samples_str += ",";
            samples_str += std::to_string(sv.tumor_samples[i]);
        }
        if (samples_str.empty()) samples_str = ".";
        
        for (size_t i = 0; i < sv.genomic_starts.size(); ++i) {
            sv_file << type_str << "\t" 
                    << (sv.chrom.empty() ? "." : sv.chrom) << "\t"
                    << sv.genomic_starts[i] << "\t" 
                    << sv.genomic_ends[i] << "\t"
                    << samples_str << "\t"
                    << sv.normal_count << "\t"
                    << sv.tumor_count << "\n";
        }
    }
    sv_file.close();
}
