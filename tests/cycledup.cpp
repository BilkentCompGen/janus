#include <cassert>
#include <iostream>
#include <string>

#include "../src/graph.h"
#include "../src/utils.h"

using namespace std;

static node_hash addNode(LCPGraph& g, const string& seq, uint32_t colors) {
    Node n;
    n.hash = hash_string(seq, 0);
    n.sequence = seq;
    n.colors = colors;
    g.nodes[n.hash] = n;
    return n.hash;
}

static void addEdge(LCPGraph& g, node_hash from, node_hash to) {
    g.next_edges[from].insert(to);
    g.prev_edges[to].insert(from);
}

void testLooseBoundaryDetection() {
    cout << "Test: loose boundary detection\n";
    LCPGraph g;

    node_hash hStart = addNode(g, "START_SHARED", 3);
    node_hash hEnd = addNode(g, "END_SHARED", 3);
    node_hash hNormalOnly = addNode(g, "N_ONLY", 1);
    node_hash hTumorOnly = addNode(g, "T_ONLY", 2);

    addEdge(g, hStart, hNormalOnly);
    addEdge(g, hStart, hTumorOnly);
    addEdge(g, hNormalOnly, hEnd);
    addEdge(g, hTumorOnly, hEnd);

    BubbleResult markers = g.detectBubbles();
    assert(markers.bubble_starts.count(hStart));
    assert(markers.bubble_ends.count(hEnd));
    cout << "PASS\n\n";
}

void testNoBoundaryOnFullySharedNeighbors() {
    cout << "Test: no loose boundary on fully shared neighbors\n";
    LCPGraph g;

    node_hash hA = addNode(g, "A_SHARED", 3);
    node_hash hB = addNode(g, "B_SHARED", 3);
    node_hash hC = addNode(g, "C_SHARED", 3);

    addEdge(g, hA, hB);
    addEdge(g, hB, hC);

    BubbleResult markers = g.detectBubbles();
    assert(markers.bubble_starts.empty());
    assert(markers.bubble_ends.empty());
    cout << "PASS\n\n";
}

void testExtractBubblesStillWorksWithLooseMarkers() {
    cout << "Test: extractBubbles with loose markers\n";
    LCPGraph g;

    node_hash hS = addNode(g, "ACGTACGT", 3);
    node_hash hN = addNode(g, "CGTGCGTG", 1);
    node_hash hT = addNode(g, "CGTTCGTG", 2);
    node_hash hE = addNode(g, "GTGATCGA", 3);

    addEdge(g, hS, hN);
    addEdge(g, hS, hT);
    addEdge(g, hN, hE);
    addEdge(g, hT, hE);

    BubbleResult markers = g.detectBubbles();
    auto bubbles = g.extractBubbles(markers);

    assert(!bubbles.empty());
    bool has_pair = false;
    for (const auto& b : bubbles) {
        if (b.start_node == hS && b.end_node == hE) {
            has_pair = true;
            break;
        }
    }
    assert(has_pair);
    cout << "PASS\n\n";
}

int main() {
    cout << "=== Loose Boundary Tests ===\n\n";
    testLooseBoundaryDetection();
    testNoBoundaryOnFullySharedNeighbors();
    testExtractBubblesStillWorksWithLooseMarkers();
    cout << "All loose boundary tests passed!\n";
    return 0;
}
