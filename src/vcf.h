#ifndef VCF_H
#define VCF_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_set>

struct DeletionRegion {
    std::string chrom;
    uint64_t start;
    uint64_t end;
};

// Forward declarations
struct LinearReference;
typedef uint64_t node_hash_t;

std::vector<DeletionRegion> parseVCFDeletions(const std::string& vcf_path);
std::unordered_set<node_hash_t> findDeletionNodeHashes(
    const std::vector<DeletionRegion>& deletions,
    const LinearReference& reference);

#endif // VCF_H

