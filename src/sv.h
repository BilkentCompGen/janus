#ifndef SV_H
#define SV_H

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include "graph.h"

class RefIndex;

struct SVRegion { int64_t start; int64_t length; std::string call;};

class SVCall {
public:
    std::string normal_seq;
    std::string tumor_seq;
    std::string svtype;
    int64_t start;
    int64_t length;
    int64_t ref_pos;

    SVCall(const std::string& normal_seq, const std::string& tumor_seq, const std::string& svtype, int64_t start, int64_t length, int64_t ref_pos);
    SVCall() : start(0), length(0), ref_pos(0) {}
};

SVRegion detectSV(const std::string& sequence, const std::string& query);
std::vector<SVCall> callSV(const Bubble& bubble, int64_t ref_offset);
std::vector<SVCall> extractAndCallSVs(const LCPGraph& graph, const RefIndex& ref_index, const BubbleResult& markers,
    const std::vector<SVCall>& large_inv_regions = {});
std::vector<SVCall> detectLargeINVs(const std::string& tumor_seq, int depth,
    const std::unordered_set<uint64_t>& ref_hashes,
    const std::unordered_set<uint64_t>& rcref_hashes,
    const RefIndex& ref_index);

std::vector<SVCall> detectLargeDELs(const RefIndex& ref_index,
    const std::unordered_set<uint64_t>& tumor_hashes,
    const std::vector<SVCall>& large_inv_regions);

std::vector<SVCall> detectLargeINS(const std::string& tumor_seq, int depth,
    const std::unordered_set<uint64_t>& ref_hashes,
    const std::unordered_set<uint64_t>& rcref_hashes,
    const RefIndex& ref_index);

#endif