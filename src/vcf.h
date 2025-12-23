#ifndef VCF_H
#define VCF_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_set>

struct DeletionRegion {
    std::string chrom;
    uint64_t start;   // 1-based inclusive
    uint64_t end;     // 1-based inclusive
};

typedef uint64_t node_hash;

class LCPGraph;

std::vector<DeletionRegion>
parseVCFDeletions(const std::string& vcf_path);

std::unordered_set<node_hash>
findDeletionNodeHashes(const std::vector<DeletionRegion>& deletions, const LCPGraph& graph);

#endif // VCF_H