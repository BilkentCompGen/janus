#ifndef REF_H
#define REF_H

#include <string>
#include <vector>
#include "graph.h"
#include "align.h"

struct VariantRegion { std::string chrom, type; uint64_t start, end; };

std::vector<VariantRegion> parseVCF(const std::string& path);
void exportAlignedVCF(const std::string& path, std::vector<VariantCall>& vars);

#endif
