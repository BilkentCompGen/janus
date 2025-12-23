#include "analysis.h"
#include "utils.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>

static Logger analysis_logger("analysis");

void CollisionTracker::setSequence(const std::string& chrom, const std::string& sequence) {
    chrom_sequences[chrom] = sequence;
}

void CollisionTracker::addHash(node_hash_t hash, const std::string& chrom, uint64_t start, uint64_t end) {
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
    size_t collisions_with_same_length = 0;
    size_t true_duplicates = 0;  // Same hash AND same sequence
    size_t hash_collisions = 0;  // Same hash BUT different sequences
    size_t total_collision_groups = 0;
    
    for (const auto& [hash, positions] : hash_positions) {
        if (positions.size() > 1) {
            collision_hashes++;
            total_collision_groups++;
            
            // Check if all positions in this collision have the same length
            bool all_same_length = true;
            uint64_t first_length = 0;
            bool first_set = false;
            
            // Check if all sequences are actually the same (true duplicates vs hash collisions)
            bool all_same_sequence = true;
            std::string first_sequence;
            bool first_seq_set = false;
            
            for (const auto& [chrom, start, end] : positions) {
                uint64_t length = end > start ? (end - start) : 0;
                if (!first_set) {
                    first_length = length;
                    first_set = true;
                } else if (length != first_length) {
                    all_same_length = false;
                }
                
                // Get sequence substring from global sequence (only when needed for collision check)
                auto seq_it = chrom_sequences.find(chrom);
                if (seq_it != chrom_sequences.end() && end > start && end <= seq_it->second.size()) {
                    std::string seq_substring = seq_it->second.substr(start, end - start);
                    
                    // Check sequence content
                    if (!first_seq_set) {
                        first_sequence = seq_substring;
                        first_seq_set = true;
                    } else if (seq_substring != first_sequence) {
                        all_same_sequence = false;
                    }
                } else {
                    // Can't verify sequence (chrom not found or invalid range)
                    all_same_sequence = false;
                }
            }
            
            if (all_same_length) {
                collisions_with_same_length++;
            }
            
            if (all_same_sequence) {
                true_duplicates++;
            } else {
                hash_collisions++;  // Same hash but different sequences = real hash collision!
            }
        }
    }
    
    double collision_rate = total_hashes > 0 ? (100.0 * collision_hashes / unique_hashes) : 0.0;
    double same_length_rate = total_collision_groups > 0 ? (100.0 * collisions_with_same_length / total_collision_groups) : 0.0;
    double true_duplicate_rate = total_collision_groups > 0 ? (100.0 * true_duplicates / total_collision_groups) : 0.0;
    double hash_collision_rate = total_collision_groups > 0 ? (100.0 * hash_collisions / total_collision_groups) : 0.0;
    
    // Calculate total duplicate instances (not just groups)
    size_t total_duplicate_instances = 0;
    size_t total_hash_collision_instances = 0;
    for (const auto& [hash, positions] : hash_positions) {
        if (positions.size() > 1) {
            // Check if this is a true duplicate or hash collision
            bool is_duplicate = true;
            std::string first_sequence;
            bool first_seq_set = false;
            
            for (const auto& [chrom, start, end] : positions) {
                auto seq_it = chrom_sequences.find(chrom);
                if (seq_it != chrom_sequences.end() && end > start && end <= seq_it->second.size()) {
                    std::string seq_substring = seq_it->second.substr(start, end - start);
                    if (!first_seq_set) {
                        first_sequence = seq_substring;
                        first_seq_set = true;
                    } else if (seq_substring != first_sequence) {
                        is_duplicate = false;
                        break;
                    }
                }
            }
            
            if (is_duplicate) {
                total_duplicate_instances += (positions.size() - 1);  // All except first are duplicates
            } else {
                total_hash_collision_instances += (positions.size() - 1);  // All except first are collisions
            }
        }
    }
    
    analysis_logger.debug("Hash collision TSV report written to " + report_file.string());
    analysis_logger.debug("Found " + std::to_string(collision_hashes) + " hashes with collisions out of " + std::to_string(unique_hashes) + " unique hashes (" + std::to_string(collision_rate) + "%)");
    analysis_logger.debug("Collisions with same length: " + std::to_string(collisions_with_same_length) + " out of " + std::to_string(total_collision_groups) + " collision groups (" + std::to_string(same_length_rate) + "%)");
    analysis_logger.debug("TRUE DUPLICATES: " + std::to_string(true_duplicates) + " groups (" + std::to_string(true_duplicate_rate) + "%), " + std::to_string(total_duplicate_instances) + " total duplicate instances");
    analysis_logger.debug("HASH COLLISIONS: " + std::to_string(hash_collisions) + " groups (" + std::to_string(hash_collision_rate) + "%), " + std::to_string(total_hash_collision_instances) + " total collision instances");
}
