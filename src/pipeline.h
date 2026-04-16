#ifndef PIPELINE_H
#define PIPELINE_H

#include <cstdint>
#include <string>
#include <unordered_set>

#include "graph.h"
#include "ref.h"

void processGenome(const std::string& file, int color, uint32_t& sid,
                   std::string& seq, LCPGraph& graph,
                   const std::string& chr_str, int depth);

std::unordered_set<uint64_t> buildRefHashes(
    const std::string& ref_file, const std::string& chr_str, int depth,
    RefIndex& ref_index, std::string& ref_sequence);

std::unordered_set<uint64_t> buildRcrefHashes(
    const std::string& ref_sequence, int depth);

std::unordered_set<uint64_t> buildTumorHashes(
    const std::string& tumor_sequence, int depth);

#endif
