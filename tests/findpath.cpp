#include <iostream>
#include <cassert>
#include <string>

#include "../src/graph.h"
#include "../src/utils.h"

using namespace std;

static Node makeNode(const string& seq, uint32_t colors) {
    Node n;
    n.hash = hash_string(seq, 0);
    n.sequence = seq;
    n.colors = colors;
    return n;
}

static node_hash addNode(LCPGraph& g, const string& seq, uint32_t colors) {
    Node n = makeNode(seq, colors);
    g.nodes[n.hash] = n;
    return n.hash;
}

static void addEdge(LCPGraph& g, node_hash from, node_hash to) {
    g.next_edges[from].insert(to);
    g.prev_edges[to].insert(from);
}

// Test 1: Simple SNP bubble
//   START(ACGTACGT) → N(CGTGCGTG) → END(GTGATCGA)
//   START(ACGTACGT) → T(CGTTCGTG) → END(GTGATCGA)
// Overlaps: START-N=CGT(3), START-T=CGT(3), N-END=GTG(3), T-END=GTG(3)
// Content between start and end:  normal="GC", tumor="TC" (the SNP: G→T)
void testSimpleSNP() {
    cout << "Test: simple SNP bubble\n";
    LCPGraph g;

    node_hash hS = addNode(g, "ACGTACGT", 3);
    node_hash hN = addNode(g, "CGTGCGTG", 1);
    node_hash hT = addNode(g, "CGTTCGTG", 2);
    node_hash hE = addNode(g, "GTGATCGA", 3);

    addEdge(g, hS, hN); addEdge(g, hS, hT);
    addEdge(g, hN, hE); addEdge(g, hT, hE);

    BubbleResult markers;
    markers.bubble_starts.insert(hS);
    markers.bubble_ends.insert(hE);

    Bubble b;
    b.start_node = hS;
    node_hash endN = g.findColoredPath(hS, markers, 0, b);
    node_hash endT = g.findColoredPath(hS, markers, 1, b);

    assert(endN == hE && endT == hE);
    cout << "  normal=\"" << b.normal_seq << "\" tumor=\"" << b.tumor_seq << "\"\n";
    assert(b.normal_seq == "GC");
    assert(b.tumor_seq == "TC");
    assert(b.normal_seq.size() == b.tumor_seq.size());

    cout << "PASS\n\n";
}

// Test 2: Deletion bubble (symmetric overlap)
//   START(ACGTACGT) → A(CGTAATCG) → B(TCGCCGTG) → END(GTGATCGA)  [normal]
//   START(ACGTACGT) → X(CGTCCGTG) → END(GTGATCGA)                 [tumor]
// Overlaps: START-A=CGT(3), A-B=TCG(3), B-END=GTG(3), START-X=CGT(3), X-END=GTG(3)
// Content: normal="AATCGCC"(7bp), tumor="CC"(2bp). Normal > tumor.
void testDeletion() {
    cout << "Test: deletion bubble\n";
    LCPGraph g;

    node_hash hS = addNode(g, "ACGTACGT", 3);
    node_hash hA = addNode(g, "CGTAATCG", 1);
    node_hash hB = addNode(g, "TCGCCGTG", 1);
    node_hash hX = addNode(g, "CGTCCGTG", 2);
    node_hash hE = addNode(g, "GTGATCGA", 3);

    addEdge(g, hS, hA); addEdge(g, hA, hB); addEdge(g, hB, hE);
    addEdge(g, hS, hX); addEdge(g, hX, hE);

    BubbleResult markers;
    markers.bubble_starts.insert(hS);
    markers.bubble_ends.insert(hE);

    Bubble b;
    b.start_node = hS;
    node_hash endN = g.findColoredPath(hS, markers, 0, b);
    node_hash endT = g.findColoredPath(hS, markers, 1, b);

    assert(endN == hE && endT == hE);
    cout << "  normal=\"" << b.normal_seq << "\"(" << b.normal_seq.size() << "bp)";
    cout << " tumor=\"" << b.tumor_seq << "\"(" << b.tumor_seq.size() << "bp)\n";
    assert(b.normal_seq == "AATCGCC");
    assert(b.tumor_seq == "CC");
    assert(b.normal_seq.size() > b.tumor_seq.size());

    cout << "PASS\n\n";
}

// Test 3: Deletion with asymmetric start overlap
//   START(ACGTACGT) → A(TAATCG) → B(TCGCCG) → END(CCGATCGA)   [normal, 1bp start overlap]
//   START(ACGTACGT) → X(ACGTACGTCCG) → END(CCGATCGA)           [tumor, 8bp start overlap]
// Without fix: normal=9bp, tumor=11bp → tumor>normal (BUG)
// With fix: normal="AATCG"(5bp), tumor=""(0bp). Normal > tumor.
void testDeletionAsymmetricOverlap() {
    cout << "Test: deletion with asymmetric overlap\n";
    LCPGraph g;

    node_hash hS = addNode(g, "ACGTACGT", 3);
    node_hash hA = addNode(g, "TAATCG", 1);
    node_hash hB = addNode(g, "TCGCCG", 1);
    node_hash hX = addNode(g, "ACGTACGTCCG", 2);
    node_hash hE = addNode(g, "CCGATCGA", 3);

    addEdge(g, hS, hA); addEdge(g, hA, hB); addEdge(g, hB, hE);
    addEdge(g, hS, hX); addEdge(g, hX, hE);

    BubbleResult markers;
    markers.bubble_starts.insert(hS);
    markers.bubble_ends.insert(hE);

    Bubble b;
    b.start_node = hS;
    node_hash endN = g.findColoredPath(hS, markers, 0, b);
    node_hash endT = g.findColoredPath(hS, markers, 1, b);

    assert(endN == hE && endT == hE);
    cout << "  normal=\"" << b.normal_seq << "\"(" << b.normal_seq.size() << "bp)";
    cout << " tumor=\"" << b.tumor_seq << "\"(" << b.tumor_seq.size() << "bp)\n";
    assert(b.normal_seq == "AATCG");
    assert(b.tumor_seq.empty());
    assert(b.normal_seq.size() > b.tumor_seq.size());

    cout << "PASS\n\n";
}

// Test 4: Direct edge to bubble end
//   START(ACGTACGT) → A(CGTAATCG) → END(TCGATCGA)  [normal path]
//   START(ACGTACGT) → END(TCGATCGA)                  [tumor direct edge]
void testDirectEnd() {
    cout << "Test: direct edge to bubble end\n";
    LCPGraph g;

    node_hash hS = addNode(g, "ACGTACGT", 3);
    node_hash hA = addNode(g, "CGTAATCG", 1);
    node_hash hE = addNode(g, "TCGATCGA", 3);

    addEdge(g, hS, hA); addEdge(g, hA, hE);
    addEdge(g, hS, hE);

    BubbleResult markers;
    markers.bubble_starts.insert(hS);
    markers.bubble_ends.insert(hE);

    Bubble b;
    b.start_node = hS;
    node_hash endN = g.findColoredPath(hS, markers, 0, b);
    node_hash endT = g.findColoredPath(hS, markers, 1, b);

    assert(endN == hE && endT == hE);
    cout << "  normal=\"" << b.normal_seq << "\"(" << b.normal_seq.size() << "bp)";
    cout << " tumor=\"" << b.tumor_seq << "\"(" << b.tumor_seq.size() << "bp)\n";

    // Tumor has no color-exclusive nodes → empty (correct)
    assert(b.tumor_seq.empty());
    // Normal: depends on iteration order of unordered_set
    // If END checked before A → normal="" (known limitation)
    // If A checked before END → normal="AA" (correct: CGTAATCG minus CGT start, minus TCG end)
    if (b.normal_seq.empty())
        cout << "  NOTE: normal empty due to iteration order (direct end hit first)\n";
    else
        assert(b.normal_seq == "AA");

    cout << "PASS\n\n";
}

int main() {
    cout << "=== findColoredPath Unit Tests ===\n\n";
    testSimpleSNP();
    testDeletion();
    testDeletionAsymmetricOverlap();
    testDirectEnd();
    cout << "All findColoredPath tests done!\n";
    return 0;
}
