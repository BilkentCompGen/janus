#include <any>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <unordered_set>

#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"

#include "src/align.h"
#include "src/graph.h"
#include "src/ref.h"
#include "src/utils.h"

using namespace std;

Logger main_logger("main");

int main(int argc, char** argv) {
    LCP_INIT();
    ArgParser parser("Build LCP graph from reference, normal, and tumor FASTA files.");
    parser.add("ref", "--ref", "-r", "Reference genome FASTA file", "string", true);
    parser.add("normal", "--normal", "-n", "Normal genome FASTA file", "string", true);
    parser.add("tumor", "--tumor", "-t", "Tumor genome FASTA file", "string", true);
    parser.add("vcf", "--vcf", "-v", "VCF file with ground truth deletions", "string", false);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("output", "--output", "-o", "Output prefix (without extension)", "string", true);

    map<string, any> args = parser.parse(argc, argv);
    string ref_file = any_cast<string>(args["ref"]), normal_file = any_cast<string>(args["normal"]), tumor_file = any_cast<string>(args["tumor"]), output_prefix = any_cast<string>(args["output"]);
    int depth = any_cast<int>(args["depth"]);
    string vcf_file = args.count("vcf") ? any_cast<string>(args["vcf"]) : "";

    namespace fs = std::filesystem;
    fs::path output_dir = "outputs"; fs::create_directories(output_dir);

    LCPGraph graph;
    uint32_t ref_seq_id = 0, normal_seq_id = 0, tumor_seq_id = 0;
    string ref_sequence, normal_sequence, tumor_sequence;

    auto process = [&](const string& file, int color, uint32_t& sid, string& seq) {
        FastaStreamer streamer(file, false);
        for (const FastaPair& pair : streamer) {
            string chrom = pair.first; size_t hap_pos = chrom.find("_hap");
            if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
            if (chrom != "19") continue;
            sid = graph.seq_names.size(); graph.seq_names.push_back(pair.first);
            seq = pair.second;
            struct lps lcp_string; init_lps(&lcp_string, seq.c_str(), seq.size()); lps_deepen(&lcp_string, depth);
            graph.addString(&lcp_string, color, seq.c_str(), sid); free_lps(&lcp_string);
        }
    };

    process(ref_file, 0, ref_seq_id, ref_sequence);
    process(normal_file, 1, normal_seq_id, normal_sequence);
    process(tumor_file, 2, tumor_seq_id, tumor_sequence);

    BubbleResult bubbles = graph.detectBubbles();
    graph.markBubbles(bubbles);
    auto extracted_bubbles = graph.extractBubbles(bubbles, ref_seq_id, normal_seq_id, tumor_seq_id);

    if (!vcf_file.empty()) {
        auto gt_vars = parseVCF(vcf_file);
        for (const auto& b : extracted_bubbles) {
            for (const auto& gt : gt_vars) {
                if (std::abs((long long)b.ref_s - (long long)gt.start) < 50) {
                    int bit = (gt.type == "DEL") ? 5 : (gt.type == "INS" ? 6 : 7);
                    bool mark_normal = (gt.type == "DEL");
                    uint32_t target_bit = mark_normal ? 0x2 : 0x4;
                    uint32_t avoid_bit = mark_normal ? 0x4 : 0x1;

                    std::queue<node_hash> q; std::unordered_set<node_hash> vis;
                    for (node_hash nh : graph.next_edges.at(b.start_node)) {
                        auto& n = graph.nodes.at(nh);
                        if ((n.colors & target_bit) && !(n.colors & avoid_bit)) { q.push(nh); vis.insert(nh); }
                    }
                    while (!q.empty()) {
                        node_hash cur = q.front(); q.pop();
                        if (cur == b.end_node) continue;
                        graph.nodes.at(cur).setColor(bit);
                        for (node_hash nh : graph.next_edges.at(cur)) {
                            if (!vis.count(nh) && nh != b.end_node && (graph.nodes.at(nh).colors & target_bit)) {
                                q.push(nh); vis.insert(nh);
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<VariantCall> all_variants;
    size_t total_bubbles = extracted_bubbles.size();
    size_t processed = 0;
    for (const auto& b : extracted_bubbles) {
        ++processed;
        if (b.ref_s == UINT64_MAX || b.norm_s == UINT64_MAX || b.tum_s == UINT64_MAX) {
            main_logger.debug("Bubble " + to_string(processed) + "/" + to_string(total_bubbles) + ": Invalid coordinates, skipping");
            continue;
        }
        if (b.norm_s >= b.norm_e || b.norm_e > normal_sequence.size() || b.tum_s >= b.tum_e || b.tum_e > tumor_sequence.size()) {
            main_logger.debug("Bubble " + to_string(processed) + "/" + to_string(total_bubbles) + ": Out of bounds, skipping");
            continue;
        }
        size_t n_len = b.norm_e - b.norm_s, t_len = b.tum_e - b.tum_s;
        main_logger.debug("Bubble " + to_string(processed) + "/" + to_string(total_bubbles) + ": n_len=" + to_string(n_len) + ", t_len=" + to_string(t_len));
        string n_seq = normal_sequence.substr(b.norm_s, n_len), t_seq = tumor_sequence.substr(b.tum_s, t_len);
        main_logger.debug("Bubble " + to_string(processed) + ": Starting alignment...");
        auto vars = alignAndCall(n_seq, t_seq, b.ref_s + 1, "19");
        main_logger.debug("Bubble " + to_string(processed) + ": Alignment done, found " + to_string(vars.size()) + " variants");
        for (auto& v : vars) all_variants.push_back(std::move(v));
    }
    exportAlignedVCF((output_dir / (output_prefix + "_variants.vcf")).string(), all_variants);

    map<uint32_t, uint32_t> color_map;
    for (uint32_t i = 0; i < 256; ++i) {
        if (i & (1 << 7)) color_map[i] = 0xFF00FF;
        else if (i & (1 << 6)) color_map[i] = 0x00FFFF;
        else if (i & (1 << 5)) color_map[i] = 0xFFFF00;
        else if (i & (1 << 3)) color_map[i] = 0x00FF00;
        else if (i & (1 << 4)) color_map[i] = 0xFF8C00;
        else if ((i & 0x07) == 0x07) color_map[i] = 0xD3D3D3;
        else if (i & 0x02) color_map[i] = 0x0000FF;
        else if (i & 0x04) color_map[i] = 0xFF0000;
        else if (i & 0x01) color_map[i] = 0x808080;
    }

    graph.exportGFA(output_dir.string(), output_prefix + "_raw", &color_map);
    graph.compressGraph();
    graph.exportGFA(output_dir.string(), output_prefix + "_compressed", &color_map);

    main_logger.debug("Done!");
    return 0;
}
