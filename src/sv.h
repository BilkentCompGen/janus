#ifndef SV_H
#define SV_H

#include "graph.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>

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

struct Bubble {
    node_hash_t start_node;  // T+N node
    node_hash_t end_node;    // T+N node
    std::vector<node_hash_t> tumor_path;
    std::vector<node_hash_t> normal_path;
    uint64_t tumor_path_length;
    uint64_t normal_path_length;
    std::string chrom;
    uint64_t start_pos;
    uint64_t end_pos;
    SVType type;
};

struct BubblePathReport {
    std::string path_type; // "TUMOR_ONLY" or "NORMAL_ONLY"
    node_hash_t start_node;
    node_hash_t end_node;
    std::vector<node_hash_t> internal_nodes; // excludes start/end T+N nodes
    uint64_t path_length;
    size_t node_count;
    std::string chrom;
    uint64_t start_pos;
    uint64_t end_pos;
};

std::vector<SVCall> detectDeletions(const LCPGraph& graph, const LinearReference& reference, int num_tumor_samples);
void exportTSV(const std::vector<SVCall>& sv_calls, const std::string& filename);

void exportBubblesTSV(const std::vector<Bubble>& bubbles, const std::string& output_dir, const std::string& prefix);
void exportBubblePathsTSV(const std::vector<BubblePathReport>& paths, const std::string& output_dir, const std::string& prefix);
void exportBubblesDetailedTSV(const std::vector<Bubble>& bubbles, const LCPGraph& graph, const std::string& output_dir, const std::string& prefix);
std::vector<Bubble> detectBubbles(
    const LCPGraph& graph,
    const LinearReference& reference,
    std::vector<BubblePathReport>* tumor_paths_out = nullptr,
    std::vector<BubblePathReport>* normal_paths_out = nullptr);

#endif // SV_H
