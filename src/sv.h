#ifndef SV_H
#define SV_H

#include "graph.h"
#include <vector>

enum class SVType {DELETION, INSERTION};

struct SVCall {
    SVType type;
    ublock node_hash;
};

std::vector<SVCall> detectPresenceAbsenceSVs(const CompressedGraph& compressed_graph);

#endif // SV_H

