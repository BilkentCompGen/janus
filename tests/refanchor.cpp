#include <iostream>
#include <cassert>
#include <string>

#include "../lib/lcptools/encoding.h"
#include "../lib/lcptools/lps.h"
#include "../src/ref.h"
#include "../src/graph.h"
#include "../src/sv.h"
#include "../src/utils.h"

using namespace std;

// Helper: parse into RefIndex
void parseRef(RefIndex& ri, const string& seq, int depth) {
    struct lps l; init_lps(&l, seq.c_str(), seq.size());
    lps_deepen(&l, depth);
    ri.addString(&l, seq.c_str());
    free_lps(&l);
}

// Helper: parse into LCPGraph with color
void parseGraph(LCPGraph& g, const string& seq, int color, uint32_t sid) {
    struct lps l; init_lps(&l, seq.c_str(), seq.size());
    lps_deepen(&l, 2);
    g.addString(&l, color, seq.c_str(), sid);
    free_lps(&l);
}

void testRefIndexBasic() {
    cout << "Testing RefIndex basic operations\n";
    LCP_INIT();

    // Diverse 100bp reference — no long single-char runs, no trivial repeats
    string ref_seq = "TACGATCGAACTGTCAGGCTAAGCTTGCATGCTAGCAGTTACGCATGCTTAACGTCGATGACTACGCAGTTAGCGCATCGATCAGTGCTAGCTTACGCTG";

    RefIndex ri;
    parseRef(ri, ref_seq, 2);

    // Parse same seq into a graph just to enumerate core hashes
    LCPGraph tmp;
    parseGraph(tmp, ref_seq, 0, 0);

    int unique = 0, repeat = 0;
    for (const auto& [h, node] : tmp.nodes) {
        int64_t pos = ri.anchor(h);
        pos >= 0 ? unique++ : repeat++;
    }
    cout << "  Cores in graph: " << tmp.nodeCount() << ", unique anchors: " << unique << ", repeated: " << repeat << "\n";
    assert(unique > 0);

    cout << "PASS: RefIndex basic\n\n";
}

void testEndToEnd() {
    cout << "Testing end-to-end: ref parsing → graph → bubbles → anchor → SV call\n";
    LCP_INIT();

    // 100bp reference — diverse, no single-char runs
    //       0         1         2         3         4         5         6         7         8         9
    //       0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    string ref_seq = "TACGATCGAACTGTCAGGCTAAGCTTGCATGCTAGCAGTTACGCATGCTTAACGTCGATGACTACGCAGTTAGCGCATCGATCAGTGCTAGCTTACGCTG";
    // normal = ref (identical)
    string norm_seq = ref_seq;
    // tumor: 3 somatic variants
    //   SNP at pos 25: T→A  (isolated, ~40bp from the close pair)
    //   DEL at pos 65-67: GCA deleted (close to next, ~7bp gap)
    //   SNP at pos 74: G→T  (close to previous)
    string tumor_seq = ref_seq.substr(0, 25) + "A" + ref_seq.substr(26, 39) + ref_seq.substr(68, 6) + "T" + ref_seq.substr(75);

    cout << "  ref:   " << ref_seq.size() << "bp\n";
    cout << "  norm:  " << norm_seq.size() << "bp (= ref)\n";
    cout << "  tumor: " << tumor_seq.size() << "bp (SNP@25, DEL@65-67, SNP@74)\n";

    // 1. Parse reference into RefIndex
    RefIndex ri;
    parseRef(ri, ref_seq, 2);

    // 2. Parse normal and tumor into graph
    LCPGraph graph;
    parseGraph(graph, norm_seq, 0, 0);
    parseGraph(graph, tumor_seq, 1, 1);

    int shared = 0, norm_only = 0, tumor_only = 0;
    for (const auto& [h, n] : graph.nodes) {
        if ((n.colors & 0x3) == 0x3) shared++;
        else if (n.colors & 0x1) norm_only++;
        else if (n.colors & 0x2) tumor_only++;
    }
    cout << "  Graph nodes: " << graph.nodeCount() << " (shared=" << shared << " norm_only=" << norm_only << " tumor_only=" << tumor_only << ")\n";

    // 3. Detect and extract bubbles
    BubbleResult bres = graph.detectBubbles();

    // Debug: show bubble start/end successor info
    for (node_hash s : bres.bubble_starts) {
        auto& succs = graph.next_edges.at(s);
        cout << "  BubbleStart successors:";
        for (node_hash n : succs) {
            auto it = graph.nodes.find(n);
            int c = it != graph.nodes.end() ? it->second.colors : -1;
            bool is_end = bres.bubble_ends.count(n);
            cout << " [col=" << c << (is_end ? " END" : "") << " seq=" << (it != graph.nodes.end() ? it->second.sequence.substr(0,15) : "?") << "]";
        }
        cout << "\n";
    }

    auto bubbles = graph.extractBubbles(bres);
    cout << "  Bubbles: " << bubbles.size() << "\n";
    assert(bubbles.size() >= 1);

    // 4. Anchor each bubble and call SVs
    int anchored = 0, total_calls = 0;
    for (const auto& b : bubbles) {
        int64_t ref_offset = ri.anchorBubble(b, graph);
        cout << "  Bubble: normal=" << b.normal_seq.size() << "bp tumor=" << b.tumor_seq.size() << "bp anchor=" << ref_offset << "\n";

        if (ref_offset < 0) continue;
        anchored++;

        auto calls = callSV(b, ref_offset);
        total_calls += calls.size();
        for (const auto& c : calls) {
            cout << "    " << c.svtype << " ref_pos=" << c.ref_pos << " len=" << c.length
                 << " ref=" << c.normal_seq << " alt=" << c.tumor_seq << "\n";
            assert(c.ref_pos > 0 && c.ref_pos < (int64_t)ref_seq.size());
        }
    }

    cout << "  Anchored: " << anchored << "/" << bubbles.size() << "\n";
    cout << "  Total SV calls: " << total_calls << "\n";
    assert(anchored > 0);
    assert(total_calls >= 1);

    cout << "PASS: end-to-end\n\n";
}

int main() {
    cout << "=== RefIndex and Anchoring Tests ===\n\n";
    testRefIndexBasic();
    testEndToEnd();
    cout << "All RefIndex tests passed!\n";
    return 0;
}
