#include "analysis.h"
#include "utils.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>

static Logger analysis_logger("analysis");

void CollisionTracker::addHash(ublock hash, const std::string& chrom, uint64_t start, uint64_t end) {
    total_hashes++;
    hash_positions[hash].emplace_back(chrom, start, end);
    if (hash_positions[hash].size() == 1) {
        unique_hashes++;
    }
}

void CollisionTracker::generateTSVReport(const std::filesystem::path& output_dir, const std::string& output_prefix) const {
    std::filesystem::path report_file = output_dir / (output_prefix + "_hash_collisions.tsv");
    std::ofstream report(report_file);
    
    if (!report) {
        std::cerr << "Error: Could not open collision TSV file: " << report_file << std::endl;
        return;
    }
    
    // Write TSV header
    report << "hash\thash_decimal\tcollision_count\tchrom\tstart\tend\n";
    
    // Write all hash positions with their collision counts
    for (const auto& [hash, positions] : hash_positions) {
        size_t collision_count = positions.size();
        
        for (const auto& [chrom, start, end] : positions) {
            report << "0x" << std::hex << std::uppercase << hash << std::dec 
                   << "\t" << hash 
                   << "\t" << collision_count
                   << "\t" << chrom
                   << "\t" << start
                   << "\t" << end
                   << "\n";
        }
    }
    
    // Calculate and log statistics
    size_t collision_hashes = 0;
    for (const auto& [hash, positions] : hash_positions) {
        if (positions.size() > 1) {
            collision_hashes++;
        }
    }
    
    double collision_rate = total_hashes > 0 ? (100.0 * collision_hashes / unique_hashes) : 0.0;
    
    analysis_logger.debug("Hash collision TSV report written to " + report_file.string());
    analysis_logger.debug("Found " + std::to_string(collision_hashes) + " hashes with collisions out of " + std::to_string(unique_hashes) + " unique hashes (" + std::to_string(collision_rate) + "%)");
}
