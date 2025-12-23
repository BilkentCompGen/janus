#include "vcf.h"
#include "sv.h"
#include <fstream>
#include <sstream>
#include <regex>

std::vector<DeletionRegion> parseVCFDeletions(const std::string& vcf_path) {
    std::vector<DeletionRegion> deletions;
    std::ifstream vcf_file(vcf_path);
    
    if (!vcf_file.is_open()) {
        return deletions;
    }
    
    std::string line;
    std::regex end_regex(R"(END=(\d+))");
    
    while (std::getline(vcf_file, line)) {
        // Skip header lines
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // Parse VCF line: CHROM POS ID REF ALT QUAL FILTER INFO ...
        std::istringstream iss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        
        if (!(iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            continue;
        }
        
        // Only process deletions
        if (alt != "<DEL>") {
            continue;
        }
        
        // Extract END from INFO field
        std::smatch match;
        if (!std::regex_search(info, match, end_regex)) {
            continue;
        }
        
        uint64_t start_pos = std::stoull(pos);
        uint64_t end_pos = std::stoull(match[1].str());
        
        if (end_pos <= start_pos) {
            continue;
        }
        
        // Remove 'chr' prefix if present
        if (chrom.size() > 3 && chrom.substr(0, 3) == "chr") {
            chrom = chrom.substr(3);
        }
        
        DeletionRegion del;
        del.chrom = chrom;
        del.start = start_pos;
        del.end = end_pos;
        deletions.push_back(del);
    }
    
    return deletions;
}

std::unordered_set<node_hash_t> findDeletionNodeHashes(
    const std::vector<DeletionRegion>& deletions,
    const LinearReference& reference) {
    std::unordered_set<node_hash_t> deletion_hashes;
    
    for (const auto& del : deletions) {
        // Find all core hashes that overlap with this deletion region
        for (size_t i = 0; i < reference.size(); ++i) {
            std::string core_chrom = reference.getChrom(i);
            uint64_t core_start = reference.getStart(i);
            uint64_t core_end = reference.getEnd(i);
            
            // Check if chromosome matches
            std::string del_chrom = del.chrom;
            if (core_chrom.size() > 3 && core_chrom.substr(0, 3) == "chr") {
                core_chrom = core_chrom.substr(3);
            }
            if (del_chrom.size() > 3 && del_chrom.substr(0, 3) == "chr") {
                del_chrom = del_chrom.substr(3);
            }
            
            if (core_chrom != del_chrom) {
                continue;
            }
            
            // Check if core overlaps with deletion region
            // Overlap: max(0, min(core_end, del.end) - max(core_start, del.start))
            uint64_t overlap_start = std::max(core_start, del.start);
            uint64_t overlap_end = std::min(core_end, del.end);
            
            if (overlap_end > overlap_start) {
                // Core overlaps with deletion, add its hash
                node_hash_t hash = reference.hashes[i];
                deletion_hashes.insert(hash);
            }
        }
    }
    
    return deletion_hashes;
}

