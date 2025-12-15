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
    std::vector<size_t> align(ublock hash) const;
    
    std::string getChrom(size_t idx) const { return idx < chroms.size() ? chroms[idx] : ""; }
    uint64_t getStart(size_t idx) const { return idx < starts.size() ? starts[idx] : 0; }
    uint64_t getEnd(size_t idx) const { return idx < ends.size() ? ends[idx] : 0; }
    size_t size() const { return hashes.size(); }
};

struct SVCall {
    SVType type;
    ublock node_hash;
    std::string chrom;
    std::vector<uint64_t> genomic_starts;
    std::vector<uint64_t> genomic_ends;
    std::vector<int> tumor_samples;
    unsigned int count = 0;  // Global count
};

std::vector<SVCall> detectDeletions(const LCPGraph& graph, const LinearReference& reference, int num_tumor_samples);
void exportTSV(const std::vector<SVCall>& sv_calls, const std::string& filename);

#endif // SV_H
