#include <any>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

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
    string ref_sequence, normal_sequence, tumor_sequence;

    // Parse reference
    FastaStreamer ref_streamer(ref_file, false);
    for (const FastaPair& pair : ref_streamer) {
        string chrom = pair.first;
        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
        if (chrom != chr_str) continue;
        ref_sequence = pair.second;
        struct lps lcp_string;
        init_lps(&lcp_string, ref_sequence.c_str(), ref_sequence.size());
        lps_deepen(&lcp_string, depth);
        ref_index.addString(&lcp_string, ref_sequence.c_str());
        free_lps(&lcp_string);
    }

    // Parse normal (= ref)
    {
        FastaStreamer streamer(ref_file, false);
        for (const FastaPair& pair : streamer) {
            string chrom = pair.first;
            size_t hap_pos = chrom.find("_hap");
            if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
            if (chrom != chr_str) continue;
            uint32_t sid = graph.seq_names.size();
            graph.seq_names.push_back(pair.first);
            normal_sequence = pair.second;
            struct lps lcp_string;
            init_lps(&lcp_string, normal_sequence.c_str(), normal_sequence.size());
            lps_deepen(&lcp_string, depth);
            graph.addString(&lcp_string, 0, normal_sequence.c_str(), sid);
            free_lps(&lcp_string);
        }
    }

    // Parse tumor
    {
        FastaStreamer streamer(tumor_file, false);
        for (const FastaPair& pair : streamer) {
            string chrom = pair.first;
            size_t hap_pos = chrom.find("_hap");
            if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
            if (chrom != chr_str) continue;
            uint32_t sid = graph.seq_names.size();
            graph.seq_names.push_back(pair.first);
            tumor_sequence = pair.second;
            struct lps lcp_string;
            init_lps(&lcp_string, tumor_sequence.c_str(), tumor_sequence.size());
            lps_deepen(&lcp_string, depth);
            graph.addString(&lcp_string, 1, tumor_sequence.c_str(), sid);
            free_lps(&lcp_string);
        }
    }

    cout << "Graph: " << graph.nodes.size() << " nodes" << endl;

    BubbleResult boundaries = graph.detectBubbles();
    cout << "Bubble starts: " << boundaries.bubble_starts.size() << endl;
    cout << "Bubble ends: " << boundaries.bubble_ends.size() << endl;

    // For each bubble start, dump topology
    cout << "\n=== BUBBLE STARTS ===" << endl;
    int idx = 0;
    for (node_hash start : boundaries.bubble_starts) {
        const Node& sn = graph.nodes.at(start);
        int64_t ref_pos = ref_index.anchor(start);

        // Check how many times this hash appears in ref
        // (anchor returns -1 if not unique or not found)
        bool unique_in_ref = (ref_pos >= 0);

        cout << "\n--- Bubble Start #" << idx++ << " ---" << endl;
        cout << "  hash=" << start << " color=" << colorStr(sn.colors)
             << " seq_len=" << sn.sequence.size()
             << " ref_pos=" << ref_pos
             << " unique_in_ref=" << (unique_in_ref ? "yes" : "no") << endl;

        // Predecessors
        auto prev_it = graph.prev_edges.find(start);
        if (prev_it != graph.prev_edges.end()) {
            cout << "  predecessors (" << prev_it->second.size() << "):" << endl;
            for (node_hash ph : prev_it->second) {
                auto it = graph.nodes.find(ph);
                if (it == graph.nodes.end()) continue;
                int64_t pref = ref_index.anchor(ph);
                cout << "    " << colorStr(it->second.colors)
                     << " seq_len=" << it->second.sequence.size()
                     << " ref_pos=" << pref
                     << " unique=" << (pref >= 0 ? "yes" : "no") << endl;
            }
        }

        // Successors (these diverge - that's why it's a bubble start)
        auto next_it = graph.next_edges.find(start);
        if (next_it != graph.next_edges.end()) {
            cout << "  successors (" << next_it->second.size() << "):" << endl;
            int normal_only = 0, tumor_only = 0, shared = 0;
            for (node_hash nh : next_it->second) {
                auto it = graph.nodes.find(nh);
                if (it == graph.nodes.end()) continue;
                bool n = (it->second.colors & 0x1), t = (it->second.colors & 0x2);
                int64_t nref = ref_index.anchor(nh);
                cout << "    " << colorStr(it->second.colors)
                     << " seq_len=" << it->second.sequence.size()
                     << " ref_pos=" << nref
                     << " unique=" << (nref >= 0 ? "yes" : "no") << endl;
                if (n && !t) normal_only++;
                else if (t && !n) tumor_only++;
                else if (n && t) shared++;
            }
            cout << "  successor pattern: N-only=" << normal_only
                 << " T-only=" << tumor_only
                 << " shared=" << shared << endl;
        }
    }

    // For each bubble end, dump topology
    cout << "\n=== BUBBLE ENDS ===" << endl;
    idx = 0;
    for (node_hash end : boundaries.bubble_ends) {
        const Node& en = graph.nodes.at(end);
        int64_t ref_pos = ref_index.anchor(end);
        bool unique_in_ref = (ref_pos >= 0);

        cout << "\n--- Bubble End #" << idx++ << " ---" << endl;
        cout << "  hash=" << end << " color=" << colorStr(en.colors)
             << " seq_len=" << en.sequence.size()
             << " ref_pos=" << ref_pos
             << " unique_in_ref=" << (unique_in_ref ? "yes" : "no") << endl;

        // Predecessors (these converge - that's why it's a bubble end)
        auto prev_it = graph.prev_edges.find(end);
        if (prev_it != graph.prev_edges.end()) {
            cout << "  predecessors (" << prev_it->second.size() << "):" << endl;
            int normal_only = 0, tumor_only = 0, shared = 0;
            for (node_hash ph : prev_it->second) {
                auto it = graph.nodes.find(ph);
                if (it == graph.nodes.end()) continue;
                bool n = (it->second.colors & 0x1), t = (it->second.colors & 0x2);
                int64_t pref = ref_index.anchor(ph);
                cout << "    " << colorStr(it->second.colors)
                     << " seq_len=" << it->second.sequence.size()
                     << " ref_pos=" << pref
                     << " unique=" << (pref >= 0 ? "yes" : "no") << endl;
                if (n && !t) normal_only++;
                else if (t && !n) tumor_only++;
                else if (n && t) shared++;
            }
            cout << "  predecessor pattern: N-only=" << normal_only
                 << " T-only=" << tumor_only
                 << " shared=" << shared << endl;
        }

        // Successors
        auto next_it = graph.next_edges.find(end);
        if (next_it != graph.next_edges.end()) {
            cout << "  successors (" << next_it->second.size() << "):" << endl;
            for (node_hash nh : next_it->second) {
                auto it = graph.nodes.find(nh);
                if (it == graph.nodes.end()) continue;
                int64_t nref = ref_index.anchor(nh);
                cout << "    " << colorStr(it->second.colors)
                     << " seq_len=" << it->second.sequence.size()
                     << " ref_pos=" << nref
                     << " unique=" << (nref >= 0 ? "yes" : "no") << endl;
            }
        }
    }

    // Summary: count successor/predecessor patterns
    cout << "\n=== PATTERN SUMMARY ===" << endl;
    map<string, int> start_patterns, end_patterns;
    for (node_hash start : boundaries.bubble_starts) {
        auto next_it = graph.next_edges.find(start);
        if (next_it == graph.next_edges.end()) continue;
        int n_only = 0, t_only = 0, sh = 0;
        for (node_hash nh : next_it->second) {
            auto it = graph.nodes.find(nh);
            if (it == graph.nodes.end()) continue;
            bool n = (it->second.colors & 0x1), t = (it->second.colors & 0x2);
            if (n && !t) n_only++;
            else if (t && !n) t_only++;
            else if (n && t) sh++;
        }
        string pat = "N=" + to_string(n_only) + " T=" + to_string(t_only) + " S=" + to_string(sh);
        start_patterns[pat]++;
    }
    for (node_hash end : boundaries.bubble_ends) {
        auto prev_it = graph.prev_edges.find(end);
        if (prev_it == graph.prev_edges.end()) continue;
        int n_only = 0, t_only = 0, sh = 0;
        for (node_hash ph : prev_it->second) {
            auto it = graph.nodes.find(ph);
            if (it == graph.nodes.end()) continue;
            bool n = (it->second.colors & 0x1), t = (it->second.colors & 0x2);
            if (n && !t) n_only++;
            else if (t && !n) t_only++;
            else if (n && t) sh++;
        }
        string pat = "N=" + to_string(n_only) + " T=" + to_string(t_only) + " S=" + to_string(sh);
        end_patterns[pat]++;
    }

    cout << "\nBubble start successor patterns:" << endl;
    for (auto& [pat, cnt] : start_patterns) cout << "  " << pat << " : " << cnt << endl;

    cout << "\nBubble end predecessor patterns:" << endl;
    for (auto& [pat, cnt] : end_patterns) cout << "  " << pat << " : " << cnt << endl;

    // How many bubble starts/ends have non-unique ref positions?
    int starts_non_unique = 0, ends_non_unique = 0;
    for (node_hash h : boundaries.bubble_starts)
        if (ref_index.anchor(h) < 0) starts_non_unique++;
    for (node_hash h : boundaries.bubble_ends)
        if (ref_index.anchor(h) < 0) ends_non_unique++;

    cout << "\nRef uniqueness:" << endl;
    cout << "  Bubble starts non-unique in ref: " << starts_non_unique << "/" << boundaries.bubble_starts.size() << endl;
    cout << "  Bubble ends non-unique in ref: " << ends_non_unique << "/" << boundaries.bubble_ends.size() << endl;

    // Extraction stats
    auto bubbles = graph.extractBubbles(boundaries);
    cout << "\nExtracted " << bubbles.size() << " bubbles from "
         << boundaries.bubble_starts.size() << " starts" << endl;

    int starts_that_extracted = 0;
    for (node_hash start : boundaries.bubble_starts) {
        bool found = false;
        for (auto& b : bubbles) {
            if (b.start_node == start) { found = true; break; }
        }
        if (found) starts_that_extracted++;
    }
    cout << "Starts that produced a bubble: " << starts_that_extracted << endl;
    cout << "Starts that failed extraction: " << (boundaries.bubble_starts.size() - starts_that_extracted) << endl;

    return 0;
}
