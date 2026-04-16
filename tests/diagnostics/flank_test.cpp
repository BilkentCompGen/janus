#include <parasail.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"
#include "src/graph.h"
#include "src/ref.h"
#include "src/sv.h"
#include "src/utils.h"
using namespace std;

struct InsCall {
    int64_t pos;
    int64_t svlen;
    string insseq;
};

struct GtDup {
    int64_t pos;
    int64_t svlen;
};

int main() {
    LCP_INIT();

    // Build graph and run pipeline to get INS calls
    string ref_file = "simulation/ref.fasta";
    string tumor_file = "simulation/duplication_simulation.fa";
    int depth = 5;
    string chr_str = "19";

    RefIndex ref_index;
    LCPGraph graph;
    string ref_seq;

    // Load ref
    FastaStreamer rs(ref_file, false);
    for (const FastaPair& p : rs) {
        string c = p.first; size_t hp = c.find("_hap"); if (hp != string::npos) c = c.substr(0, hp);
        if (c != chr_str) continue;
        ref_seq = p.second;
        struct lps l; init_lps(&l, ref_seq.c_str(), ref_seq.size()); lps_deepen(&l, depth);
        ref_index.addString(&l, ref_seq.c_str()); free_lps(&l);
    }
    // Normal
    { FastaStreamer s(ref_file, false);
      for (const FastaPair& p : s) {
        string c = p.first; size_t hp = c.find("_hap"); if (hp != string::npos) c = c.substr(0, hp);
        if (c != chr_str) continue;
        uint32_t sid = graph.seq_names.size(); graph.seq_names.push_back(p.first);
        string seq = p.second;
        struct lps l; init_lps(&l, seq.c_str(), seq.size()); lps_deepen(&l, depth);
        graph.addString(&l, 0, seq.c_str(), sid); free_lps(&l);
    }}
    // Tumor
    { FastaStreamer s(tumor_file, false);
      for (const FastaPair& p : s) {
        string c = p.first; size_t hp = c.find("_hap"); if (hp != string::npos) c = c.substr(0, hp);
        if (c != chr_str) continue;
        uint32_t sid = graph.seq_names.size(); graph.seq_names.push_back(p.first);
        string seq = p.second;
        struct lps l; init_lps(&l, seq.c_str(), seq.size()); lps_deepen(&l, depth);
        graph.addString(&l, 1, seq.c_str(), sid); free_lps(&l);
    }}

    BubbleResult markers = graph.detectBubbles();
    auto calls = extractAndCallSVs(graph, ref_index, markers);

    // Collect INS calls with insseq
    vector<InsCall> ins_calls;
    for (auto& c : calls) {
        if (c.svtype == "INS" && !c.tumor_seq.empty() && c.length > 0 && c.length <= 2000) {
            ins_calls.push_back({c.ref_pos, c.length, c.tumor_seq});
        }
    }
    cout << "INS calls to test: " << ins_calls.size() << endl;

    // Load GT DUPs for comparison
    vector<GtDup> gt_dups;
    ifstream gt("simulation/duplication_simulation.vcf");
    string line;
    while (getline(gt, line)) {
        if (line[0] == '#') continue;
        // Parse chr
        size_t tab1 = line.find('\t');
        string chrom = line.substr(0, tab1);
        if (chrom.find("19") == string::npos) continue;
        // Parse pos
        size_t tab2 = line.find('\t', tab1+1);
        int64_t pos = stoll(line.substr(tab1+1, tab2-tab1-1));
        // Parse SVLEN
        size_t svlen_pos = line.find("SVLEN=");
        if (svlen_pos == string::npos) continue;
        int64_t svlen = abs(stoll(line.substr(svlen_pos+6)));
        // Check SVTYPE=DUP
        if (line.find("SVTYPE=DUP") == string::npos) continue;
        gt_dups.push_back({pos, svlen});
    }

    // Parasail setup
    int match = 2, mismatch = -4, gap_open = 5, gap_extend = 1;
    parasail_matrix_t* matrix = parasail_matrix_create("ACGT", match, mismatch);

    // Test: for each INS call, align tumor_seq to nearby reference window
    int old_method_pass = 0, new_method_pass = 0;
    int both_pass = 0, new_only = 0, old_only = 0, neither = 0;

    cout << "\n=== TESTING EACH INS CALL ===" << endl;
    cout << "pos       | insseq_len | window | old_sim | new_score | new_norm | old_pass | new_pass | gt_match" << endl;

    for (auto& ins : ins_calls) {
        int64_t k = min(ins.svlen, (int64_t)ins.insseq.size());
        
        // Old method: character-by-character flank comparison
        double old_best = 0.0;
        if (k > 0 && ins.pos >= k && ins.pos <= (int64_t)ref_seq.size()) {
            string prev_k = ref_seq.substr(ins.pos - k, k);
            size_t len = min(ins.insseq.size(), prev_k.size());
            size_t matches = 0;
            for (size_t i = 0; i < len; i++) {
                if (toupper(ins.insseq[i]) == toupper(prev_k[i])) matches++;
            }
            old_best = max(old_best, (double)matches / len);
        }
        if (k > 0 && ins.pos >= 0 && ins.pos + k <= (int64_t)ref_seq.size()) {
            string next_k = ref_seq.substr(ins.pos, k);
            size_t len = min(ins.insseq.size(), next_k.size());
            size_t matches = 0;
            for (size_t i = 0; i < len; i++) {
                if (toupper(ins.insseq[i]) == toupper(next_k[i])) matches++;
            }
            old_best = max(old_best, (double)matches / len);
        }
        bool old_pass = old_best >= 0.90;

        // New method: Smith-Waterman local alignment against reference window
        int64_t window_pad = max((int64_t)200, 2 * k);
        int64_t win_start = max((int64_t)0, ins.pos - window_pad);
        int64_t win_end = min((int64_t)ref_seq.size(), ins.pos + window_pad);
        string window = ref_seq.substr(win_start, win_end - win_start);

        parasail_result_t* result = parasail_sw_trace(
            ins.insseq.c_str(), (int)ins.insseq.size(),
            window.c_str(), (int)window.size(),
            gap_open, gap_extend, matrix);

        int sw_score = result->score;
        int perfect_score = (int)ins.insseq.size() * match;
        double normalized = (perfect_score > 0) ? (double)sw_score / perfect_score : 0.0;
        parasail_result_free(result);

        bool new_pass = normalized >= 0.70;

        // Check GT match
        bool gt_match = false;
        for (auto& g : gt_dups) {
            if (abs(ins.pos - g.pos) <= 200) { gt_match = true; break; }
        }

        if (old_pass) old_method_pass++;
        if (new_pass) new_method_pass++;
        if (old_pass && new_pass) both_pass++;
        else if (new_pass && !old_pass) new_only++;
        else if (old_pass && !new_pass) old_only++;
        else neither++;

        printf("%9ld | %10ld | %6ld | %7.3f | %9d | %8.3f | %8s | %8s | %s\n",
            ins.pos, (long)ins.insseq.size(), win_end - win_start,
            old_best, sw_score, normalized,
            old_pass ? "PASS" : "fail",
            new_pass ? "PASS" : "fail",
            gt_match ? "YES" : "no");
    }

    cout << "\n=== SUMMARY ===" << endl;
    cout << "Total INS calls tested: " << ins_calls.size() << endl;
    cout << "Old method (flank sim >= 0.90): " << old_method_pass << endl;
    cout << "New method (SW norm >= 0.70):   " << new_method_pass << endl;
    cout << "Both pass:   " << both_pass << endl;
    cout << "New only:    " << new_only << endl;
    cout << "Old only:    " << old_only << endl;
    cout << "Neither:     " << neither << endl;

    // Try different thresholds
    cout << "\n=== THRESHOLD SWEEP ===" << endl;
    for (double thresh : {0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90}) {
        int pass = 0;
        for (auto& ins : ins_calls) {
            int64_t window_pad = max((int64_t)200, 2 * (int64_t)ins.insseq.size());
            int64_t win_start = max((int64_t)0, ins.pos - window_pad);
            int64_t win_end = min((int64_t)ref_seq.size(), ins.pos + window_pad);
            string window = ref_seq.substr(win_start, win_end - win_start);

            parasail_result_t* result = parasail_sw_trace(
                ins.insseq.c_str(), (int)ins.insseq.size(),
                window.c_str(), (int)window.size(),
                gap_open, gap_extend, matrix);
            int sw_score = result->score;
            int perfect = (int)ins.insseq.size() * match;
            double norm = (perfect > 0) ? (double)sw_score / perfect : 0.0;
            parasail_result_free(result);

            if (norm >= thresh) pass++;
        }
        printf("  threshold=%.2f: %d/%zu pass\n", thresh, pass, ins_calls.size());
    }

    parasail_matrix_free(matrix);
    return 0;
}
