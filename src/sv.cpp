#include "sv.h"
#include "utils.h"
#include <algorithm>
#include <set>
#include <fstream>
#include <sstream>
#include <tuple>

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
