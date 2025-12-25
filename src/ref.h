#ifndef REF_H
#define REF_H

#include <string>
#include <vector>
#include <unordered_set>
#include "graph.h"
#include "align.h"

struct DeletionRegion { std::string chrom; uint64_t start, end; };

std::vector<DeletionRegion> parseVCFDeletions(const std::string& path);
std::unordered_set<node_hash> findDeletionNodeHashes(const std::vector<DeletionRegion>& dels, const LCPGraph& g, uint32_t rsid);
void exportAlignedVCF(const std::string& path, std::vector<VariantCall>& vars);

#endif
