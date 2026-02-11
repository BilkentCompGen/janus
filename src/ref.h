#ifndef REF_H
#define REF_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "core.h"
#include "lps.h"
#include "graph.h"

struct RefCore { node_hash hash; uint64_t start, end; };

class RefIndex {
public:
    void addString(struct lps* l, const char* seq);
    int64_t anchor(node_hash h) const;
    int64_t anchorBubble(const Bubble& bubble, const LCPGraph& graph) const;
private:
    std::vector<RefCore> cores;
    std::unordered_map<node_hash, std::vector<size_t>> index;
};

#endif
