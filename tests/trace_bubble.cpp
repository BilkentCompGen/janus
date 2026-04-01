#include <any>
#include <cstdint>
#include <iostream>
#include <string>
#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"
#include "src/graph.h"
#include "src/ref.h"
#include "src/sv.h"
#include "src/utils.h"
using namespace std;

string colorStr(uint32_t c) {
    bool n = (c & 0x1), t = (c & 0x2);
    if (n && t) return "N+T";
    if (n) return "N";
    if (t) return "T";
    return "?";
}

int main() {
    LCP_INIT();
    string ref_file = "simulation/ref.fasta";
    string tumor_file = "simulation/duplication_simulation.fa";
    int depth = 5;
    string chr_str = "19";

    RefIndex ref_index;
    LCPGraph graph;
    string ref_seq, normal_seq, tumor_seq;

    int64_t target = 40745261;
    int64_t target_end = 40745315; // +54

    // Build ref index
    FastaStreamer rs(ref_file, false);
    for (const FastaPair& p : rs) {
        string c = p.first; size_t hp = c.find("_hap"); if (hp != string::npos) c = c.substr(0, hp);
        if (c != chr_str) continue;
        ref_seq = p.second;
        struct lps l; init_lps(&l, ref_seq.c_str(), ref_seq.size()); lps_deepen(&l, depth);
        ref_index.addString(&l, ref_seq.c_str()); free_lps(&l);
    }

    // Dump ref cores near target
    {
        struct lps l; init_lps(&l, ref_seq.c_str(), ref_seq.size()); lps_deepen(&l, depth);
        cout << "=== REF CORES near DUP (pos=" << target << " end=" << target_end << " len=54) ===" << endl;
        for (int i = 0; i < l.size; i++) {
            int64_t cs = l.cores[i].start, ce = l.cores[i].end;
            if (ce < target - 400 || cs > target_end + 400) continue;
            string region(ref_seq.c_str() + cs, ce - cs);
            uint64_t h = hash_string(region, 0);
            string marker = "";
            if (cs <= target && ce >= target) marker += " <<<SPANS DUP START";
            if (cs <= target_end && ce >= target_end) marker += " <<<SPANS DUP END";
            if (cs >= target && ce <= target_end) marker += " <<<INSIDE DUP";
            cout << "  [" << i << "] " << cs << "-" << ce << " len=" << (ce-cs)
                 << " hash=" << h << marker << endl;
            cout << "       seq=" << region.substr(0, 80) << (region.size() > 80 ? "..." : "") << endl;
        }
        free_lps(&l);
    }

    // Build normal graph
    { FastaStreamer s(ref_file, false);
      for (const FastaPair& p : s) {
        string c = p.first; size_t hp = c.find("_hap"); if (hp != string::npos) c = c.substr(0, hp);
        if (c != chr_str) continue;
        uint32_t sid = graph.seq_names.size(); graph.seq_names.push_back(p.first);
        normal_seq = p.second;
        struct lps l; init_lps(&l, normal_seq.c_str(), normal_seq.size()); lps_deepen(&l, depth);
        graph.addString(&l, 0, normal_seq.c_str(), sid); free_lps(&l);
    }}

    // Build tumor graph and dump cores
    { FastaStreamer s(tumor_file, false);
      for (const FastaPair& p : s) {
        string c = p.first; size_t hp = c.find("_hap"); if (hp != string::npos) c = c.substr(0, hp);
        if (c != chr_str) continue;
        uint32_t sid = graph.seq_names.size(); graph.seq_names.push_back(p.first);
        tumor_seq = p.second;
        struct lps l; init_lps(&l, tumor_seq.c_str(), tumor_seq.size()); lps_deepen(&l, depth);

        // Find the tumor offset: DUPs before target shift coordinates
        // Search for the DUP region sequence in tumor
        string dup_region = ref_seq.substr(target, 54);
        size_t tumor_pos = tumor_seq.find(dup_region);
        cout << "\n=== TUMOR INFO ===" << endl;
        cout << "DUP region sequence: " << dup_region << endl;
        if (tumor_pos != string::npos) {
            // Check if it appears twice (tandem dup)
            size_t second = tumor_seq.find(dup_region, tumor_pos + 1);
            cout << "Found in tumor at: " << tumor_pos << endl;
            if (second != string::npos) {
                cout << "Second copy at: " << second << endl;
                cout << "Gap between copies: " << (second - tumor_pos) << "bp" << endl;
            }

            // Show tumor cores in this region
            cout << "\n=== TUMOR CORES near DUP region (tumor coords) ===" << endl;
            for (int i = 0; i < l.size; i++) {
                int64_t cs = l.cores[i].start, ce = l.cores[i].end;
                if (ce < (int64_t)tumor_pos - 400 || cs > (int64_t)tumor_pos + 400) continue;
                string region(tumor_seq.c_str() + cs, ce - cs);
                uint64_t h = hash_string(region, 0);
                auto nit = graph.nodes.find(h);
                // Color will be set after addString, but we can check if hash matches a ref core
                cout << "  [" << i << "] " << cs << "-" << ce << " len=" << (ce-cs)
                     << " hash=" << h << endl;
                cout << "       seq=" << region.substr(0, 80) << (region.size() > 80 ? "..." : "") << endl;
            }
        }

        graph.addString(&l, 1, tumor_seq.c_str(), sid); free_lps(&l);
    }}

    // Now show node colors for all cores near the region
    cout << "\n=== NODE COLORS (after both genomes added) ===" << endl;
    {
        struct lps l; init_lps(&l, ref_seq.c_str(), ref_seq.size()); lps_deepen(&l, depth);
        for (int i = 0; i < l.size; i++) {
            int64_t cs = l.cores[i].start, ce = l.cores[i].end;
            if (ce < target - 400 || cs > target_end + 400) continue;
            string region(ref_seq.c_str() + cs, ce - cs);
            uint64_t h = hash_string(region, 0);
            auto nit = graph.nodes.find(h);
            string col = nit != graph.nodes.end() ? colorStr(nit->second.colors) : "NOT_IN_GRAPH";
            cout << "  ref[" << cs << "-" << ce << "] hash=" << h << " color=" << col << endl;
        }
        free_lps(&l);
    }

    // Detect bubbles and trace the one near target
    BubbleResult markers = graph.detectBubbles();
    auto bubbles = graph.extractBubbles(markers);

    for (const auto& b : bubbles) {
        int64_t ref_offset = ref_index.anchorBubble(b, graph);
        if (ref_offset < 0) continue;
        if (abs(ref_offset - target) > 500) continue;

        cout << "\n=== EXTRACTED BUBBLE (ref_offset=" << ref_offset << ") ===" << endl;

        // Start node
        const Node& sn = graph.nodes.at(b.start_node);
        int64_t spos = ref_index.anchor(b.start_node);
        cout << "\nSTART NODE:" << endl;
        cout << "  hash=" << b.start_node << " color=" << colorStr(sn.colors)
             << " ref_pos=" << spos << " seq_len=" << sn.sequence.size() << endl;
        cout << "  seq=" << sn.sequence.substr(0, 100) << (sn.sequence.size() > 100 ? "..." : "") << endl;

        // Successors of start
        auto nxt = graph.next_edges.find(b.start_node);
        if (nxt != graph.next_edges.end()) {
            cout << "  SUCCESSORS:" << endl;
            for (node_hash nh : nxt->second) {
                const Node& nn = graph.nodes.at(nh);
                int64_t np = ref_index.anchor(nh);
                cout << "    hash=" << nh << " color=" << colorStr(nn.colors)
                     << " ref_pos=" << np << " seq_len=" << nn.sequence.size() << endl;
                cout << "    seq=" << nn.sequence.substr(0, 100) << (nn.sequence.size() > 100 ? "..." : "") << endl;
            }
        }

        // End node
        const Node& en = graph.nodes.at(b.end_node);
        int64_t epos = ref_index.anchor(b.end_node);
        cout << "\nEND NODE:" << endl;
        cout << "  hash=" << b.end_node << " color=" << colorStr(en.colors)
             << " ref_pos=" << epos << " seq_len=" << en.sequence.size() << endl;
        cout << "  seq=" << en.sequence.substr(0, 100) << (en.sequence.size() > 100 ? "..." : "") << endl;

        // Predecessors of end
        auto prv = graph.prev_edges.find(b.end_node);
        if (prv != graph.prev_edges.end()) {
            cout << "  PREDECESSORS:" << endl;
            for (node_hash ph : prv->second) {
                const Node& pn = graph.nodes.at(ph);
                int64_t pp = ref_index.anchor(ph);
                cout << "    hash=" << ph << " color=" << colorStr(pn.colors)
                     << " ref_pos=" << pp << " seq_len=" << pn.sequence.size() << endl;
                cout << "    seq=" << pn.sequence.substr(0, 100) << (pn.sequence.size() > 100 ? "..." : "") << endl;
            }
        }

        // Extracted sequences
        cout << "\nEXTRACTED SEQUENCES:" << endl;
        cout << "  normal_seq (" << b.normal_seq.size() << "bp): " << b.normal_seq << endl;
        cout << "  tumor_seq  (" << b.tumor_seq.size() << "bp): " << b.tumor_seq << endl;

        // Classification
        cout << "\nCLASSIFICATION:" << endl;
        cout << "  normal(" << b.normal_seq.size() << ") vs tumor(" << b.tumor_seq.size() << ")" << endl;
        if (b.normal_seq.size() > b.tumor_seq.size())
            cout << "  normal > tumor → treated as DEL (normal has 'extra' sequence)" << endl;
        else if (b.tumor_seq.size() > b.normal_seq.size())
            cout << "  tumor > normal → treated as INS" << endl;
        else
            cout << "  same size → SNP/multi" << endl;

        // Reference context
        cout << "\nREFERENCE CONTEXT:" << endl;
        cout << "  DUP region ref[" << target << ".." << target_end << "]: " << ref_seq.substr(target, 54) << endl;
    }

    return 0;
}

// APPEND: Check edges between the two tumor-only cores
// hash 7856366017321890231 and 15783397869307142715
