#include "sv.h"

bool isInNormal(const CompressedNode& node) {
    return node.colors[0] > 0;
}

bool isInTumor(const CompressedNode& node) {
    for (size_t i = 1; i < MAX_COLORS; ++i) {
        if (node.colors[i] > 0) {
            return true;
        }
    }
    return false;
}

std::vector<SVCall> detectPresenceAbsenceSVs(const CompressedGraph& compressed_graph) {
    std::vector<SVCall> sv_calls;
    const auto& nodes = compressed_graph.getNodes();
    
    for (const auto& [hash, node] : nodes) {
        bool normal_present = isInNormal(node);
        bool tumor_present = isInTumor(node);
        
        if (normal_present && !tumor_present) {
            SVCall call;
            call.type = SVType::DELETION;
            call.node_hash = hash;
            sv_calls.push_back(call);
        }
        
        if (!normal_present && tumor_present) {
            SVCall call;
            call.type = SVType::INSERTION;
            call.node_hash = hash;
            sv_calls.push_back(call);
        }
    }
    
    return sv_calls;
}

