#ifndef SV_H
#define SV_H

#include "graph.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>

enum class SVType { DELETION };

struct LinearReference {
    std::vector<node_hash_t> hashes;
    std::vector<std::string> chroms;
    std::vector<uint64_t> starts;
    std::vector<uint64_t> ends;
    std::unordered_map<node_hash_t, std::vector<size_t>> hash_positions;
    
    void addCore(node_hash_t hash, const std::string& chrom, uint64_t start, uint64_t end);
    void buildIndex();
    std::vector<size_t> align(node_hash_t hash) const;
    
    std::string getChrom(size_t idx) const { return idx < chroms.size() ? chroms[idx] : ""; }
    uint64_t getStart(size_t idx) const { return idx < starts.size() ? starts[idx] : 0; }
    uint64_t getEnd(size_t idx) const { return idx < ends.size() ? ends[idx] : 0; }
    size_t size() const { return hashes.size(); }
};

struct SVCall {
    SVType type;
    node_hash_t node_hash;
    std::string chrom;
    std::vector<uint64_t> genomic_starts;
    std::vector<uint64_t> genomic_ends;
    std::vector<int> tumor_samples;
    unsigned int count = 0;  // Global count
};

// Simple path info - used for both bubble paths and standalone path reports
struct PathInfo {
    std::vector<node_hash_t> internal_nodes; // excludes start/end T+N nodes
    uint64_t path_length;
    size_t node_count;
    std::string chrom;
    uint64_t start_pos;
    uint64_t end_pos;
};

struct Bubble {
    node_hash_t start_node;  // T+N node
    node_hash_t end_node;    // T+N node
    std::vector<PathInfo> tumor_paths;   // All T-only paths: T1, T2, T3...
    std::vector<PathInfo> normal_paths; // All N-only paths: N1, N2...
    std::string chrom;
    uint64_t start_pos;  // Overall bubble span (min of all paths)
    uint64_t end_pos;    // Overall bubble span (max of all paths)
    SVType type;
};

// Legacy struct for backward compatibility (can be removed later)
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

std::vector<SVCall> detectDeletionsFromBubbles(const std::vector<Bubble>& bubbles, const LCPGraph& graph, const LinearReference& reference);
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
