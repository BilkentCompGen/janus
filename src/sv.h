#ifndef SV_H
#define SV_H

#include <cstdint>
#include <string>
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
};

SVRegion detectSV(const std::string& sequence, const std::string& query);
std::vector<SVCall> callSV(const Bubble& bubble, int64_t ref_offset);
std::vector<SVCall> extractAndCallSVs(const LCPGraph& graph, const RefIndex& ref_index, const BubbleResult& markers);

#endif