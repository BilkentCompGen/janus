#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "graph.h"
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <filesystem>

// Track hash collisions during reference building
struct CollisionTracker {
    std::map<node_hash_t, std::vector<std::tuple<std::string, uint64_t, uint64_t>>> hash_positions;  // chrom, start, end
    std::map<std::string, std::string> chrom_sequences;  // chrom -> full sequence
    size_t total_hashes = 0;
    size_t unique_hashes = 0;
    
    void setSequence(const std::string& chrom, const std::string& sequence);
    void addHash(node_hash_t hash, const std::string& chrom, uint64_t start, uint64_t end);
    void generateTSVReport(const std::filesystem::path& output_dir, const std::string& output_prefix) const;
};

#endif // ANALYSIS_H
