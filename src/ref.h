#ifndef REF_H
#define REF_H

#include <cstdint>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

typedef uint64_t node_hash;

class LCPGraph;

struct DeletionRegion {
    std::string chrom;
    uint64_t start;   // 1-based inclusive
    uint64_t end;     // 1-based inclusive
};

struct Region {
    uint64_t start;
    uint64_t end;
};

struct VariantRegion {
    std::string chrom;
    uint64_t ref_start;
    uint64_t ref_end;
    std::string svtype;  // DEL, INS
};

// VCF parsing
std::vector<DeletionRegion> parseVCFDeletions(const std::string& vcf_path);

// Deletion node marking
std::unordered_set<node_hash> findDeletionNodeHashes(
    const std::vector<DeletionRegion>& deletions, 
    const LCPGraph& graph, 
    uint32_t ref_seq_id
);

// Interval merging
std::vector<Region> mergeOverlappingRegions(std::vector<Region>& regions);

// Variant collection
std::vector<VariantRegion> collectVariantRegions(
    const LCPGraph& graph, 
    uint32_t ref_seq_id,
    const std::string& chrom
);

// Anchor-based coordinate mapping
std::pair<uint64_t, uint64_t> findAnchorWindow(
    node_hash target, 
    const LCPGraph& graph,
    uint32_t ref_seq_id
);

// VCF export (old - region based)
void exportVCF(
    const std::string& vcf_path,
    const std::vector<VariantRegion>& regions
);

// Forward declaration
struct VariantCall;

// VCF export (new - alignment based)
void exportAlignedVCF(
    const std::string& vcf_path,
    const std::vector<VariantCall>& variants
);

#endif // REF_H
