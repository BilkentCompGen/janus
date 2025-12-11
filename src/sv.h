#ifndef SV_H
#define SV_H

#include "graph.h"
#include <vector>
#include <unordered_map>

enum class SVType { DELETION };

struct LinearReference {
    std::vector<ublock> hashes;
    std::vector<std::string> chroms;
    std::vector<uint64_t> starts;
    std::vector<uint64_t> ends;
    std::unordered_map<ublock, std::vector<size_t>> hash_positions;
    
    void addCore(ublock hash, const std::string& chrom, uint64_t start, uint64_t end);
    
    void buildIndex();
    
    std::vector<size_t> align(const std::vector<ublock>& path) const;
    
    std::string getChrom(size_t idx) const { return idx < chroms.size() ? chroms[idx] : ""; }
    uint64_t getStart(size_t idx) const { return idx < starts.size() ? starts[idx] : 0; }
    uint64_t getEnd(size_t idx) const { return idx < ends.size() ? ends[idx] : 0; }
    
    size_t size() const { return hashes.size(); }
};

struct SVCall {
    SVType type;
    ublock node_hash;
    std::vector<ublock> path;
    std::vector<ublock> flanking_before;
    std::vector<ublock> flanking_after;
    std::vector<size_t> reference_indices;      // Indices in linear reference
    std::string chrom;                          // Chromosome
    std::vector<uint64_t> genomic_starts;       // Real genomic start positions
    std::vector<uint64_t> genomic_ends;         // Real genomic end positions
    std::vector<int> tumor_samples;             // Which tumor samples have this SV
    unsigned int normal_count = 0;              // Count in normal sample
    unsigned int tumor_count = 0;               // Count in tumor sample (for duplications)
};

std::vector<SVCall> detectDeletions(const LCPGraph& graph, const LinearReference& reference, int num_tumor_samples);

void exportTSV(const std::vector<SVCall>& sv_calls, const std::string& filename);

bool isInNormal(const CompressedNode& node);

#endif // SV_H
