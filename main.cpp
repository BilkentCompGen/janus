#include <any>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <queue>
#include <string>
#include <unordered_set>

#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"

#include "src/align.h"
#include "src/debug.h"
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
    parser.add("debug", "--debug", "-D", "Enable verbose debug output with coloring", "bool", false);
    parser.add("tags", "--tags", "-T", "Tag file with positions to highlight (format: chrom:pos or chrom\tpos)", "string", false);

    map<string, any> args = parser.parse(argc, argv);
    string ref_file = any_cast<string>(args["ref"]), normal_file = any_cast<string>(args["normal"]), tumor_file = any_cast<string>(args["tumor"]), output_prefix = any_cast<string>(args["output"]);
    int depth = any_cast<int>(args["depth"]);
    string vcf_file = args.count("vcf") ? any_cast<string>(args["vcf"]) : "";
    string tag_file = args.count("tags") ? any_cast<string>(args["tags"]) : "";
    bool debug_mode = args.count("debug") && any_cast<bool>(args["debug"]);

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
                    uint32_t target_bit = mark_normal ? 0x2 : 0x4, avoid_bit = mark_normal ? 0x4 : 0x1;

                    std::queue<node_hash> q; std::unordered_set<node_hash> vis;
                    auto start_it = graph.next_edges.find(b.start_node);
                    if (start_it != graph.next_edges.end()) {
                        for (node_hash nh : start_it->second) {
                            auto nit = graph.nodes.find(nh);
                            if (nit != graph.nodes.end() && (nit->second.colors & target_bit) && !(nit->second.colors & avoid_bit)) { q.push(nh); vis.insert(nh); }
                        }
                    }
                    while (!q.empty()) {
                        node_hash cur = q.front(); q.pop();
                        if (cur == b.end_node) continue;
                        auto cur_it = graph.nodes.find(cur);
                        if (cur_it != graph.nodes.end()) cur_it->second.setColor(bit);
                        auto next_it = graph.next_edges.find(cur);
                        if (next_it != graph.next_edges.end()) {
                            for (node_hash nh : next_it->second) {
                                auto nit = graph.nodes.find(nh);
                                if (!vis.count(nh) && nh != b.end_node && nit != graph.nodes.end() && (nit->second.colors & target_bit) && !(nit->second.colors & avoid_bit)) { q.push(nh); vis.insert(nh); }
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<VariantCall> all_variants;
    size_t total_bubbles = extracted_bubbles.size(), processed = 0;
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

        // Only skip full alignment if one of the lengths is > 30k bp
        // Precise boundary scanning is O(N) and better for clean indel bubbles
        if (n_len > 30000 || t_len > 30000) {
            string n_seq = normal_sequence.substr(b.norm_s, n_len);
            string t_seq = tumor_sequence.substr(b.tum_s, t_len);
            
            // 1. Find longest common prefix
            size_t prefix_len = 0;
            while (prefix_len < n_len && prefix_len < t_len && n_seq[prefix_len] == t_seq[prefix_len]) {
                prefix_len++;
            }
            
            // 2. Find longest common suffix (stopping at prefix)
            size_t suffix_len = 0;
            while (suffix_len < n_len - prefix_len && suffix_len < t_len - prefix_len && 
                   n_seq[n_len - 1 - suffix_len] == t_seq[t_len - 1 - suffix_len]) {
                suffix_len++;
            }
            
            // 3. Calculate exact event size and position
            int64_t exact_n_diff = n_len - prefix_len - suffix_len;
            int64_t exact_t_diff = t_len - prefix_len - suffix_len;
            int64_t size_diff = exact_n_diff - exact_t_diff;

            // Capture even small events if they fall into this path
            if (std::abs(size_diff) >= 1) {
                if (size_diff > 0) all_variants.push_back({"19", "N", "<DEL>", "DEL", b.ref_s + 1 + prefix_len, -size_diff});
                else all_variants.push_back({"19", "N", "<INS>", "INS", b.ref_s + 1 + prefix_len, -size_diff});
                main_logger.debug("Bubble " + to_string(processed) + ": Large bubble, precise scan found prefix=" + to_string(prefix_len) + ", suffix=" + to_string(suffix_len) + ", SV size=" + to_string(std::abs(size_diff)));
            }
        } else {
            string n_seq = normal_sequence.substr(b.norm_s, n_len), t_seq = tumor_sequence.substr(b.tum_s, t_len);
            main_logger.debug("Bubble " + to_string(processed) + ": Starting full alignment...");
            auto vars = alignAndCall(n_seq, t_seq, b.ref_s + 1, "19");
            main_logger.debug("Bubble " + to_string(processed) + ": Alignment done, found " + to_string(vars.size()) + " variants");
            for (auto& v : vars) all_variants.push_back(std::move(v));
        }
    }
    fs::path vcf_dir = output_dir / "vcf"; fs::create_directories(vcf_dir);
    fs::path report_dir = output_dir / "report"; fs::create_directories(report_dir);
    exportAlignedVCF((vcf_dir / (output_prefix + "_variants.vcf")).string(), all_variants);

    std::map<uint32_t, uint32_t> color_map;
    std::map<node_hash, uint32_t> node_color_map;
    bool use_tagged_colors = !tag_file.empty();

    if (use_tagged_colors) {
        auto tags = debug::parseTagFile(tag_file);
        auto tagged_nodes = debug::findNodesNearPositions(graph, extracted_bubbles, tags, 100);
        node_color_map = debug::buildTaggedColorMap(graph, tagged_nodes);
    } else {
        color_map = debug::buildColorMap();
    }

    if (debug_mode) {
        string report_path = (report_dir / (output_prefix + "_report.txt")).string();
        // Clear previous report if exists
        std::ofstream(report_path, std::ios::trunc);
        
        auto bubble_stats = debug::analyzeBubbles(extracted_bubbles, normal_sequence, tumor_sequence);
        debug::printBubbleStats(bubble_stats, report_path);
        auto color_stats = debug::collectColorStats(graph);
        debug::printColorStats(color_stats, report_path);
        if (!vcf_file.empty()) debug::runFullAnalysis(vcf_file, (vcf_dir / (output_prefix + "_variants.vcf")).string(), "19", report_path);
    }

    if (use_tagged_colors) graph.exportGFA(output_dir.string(), output_prefix + "_raw", node_color_map);
    else graph.exportGFA(output_dir.string(), output_prefix + "_raw", &color_map);

    graph.compressGraph();

    if (use_tagged_colors) graph.exportGFA(output_dir.string(), output_prefix + "_compressed", node_color_map);
    else graph.exportGFA(output_dir.string(), output_prefix + "_compressed", &color_map);

    main_logger.debug("Done!");
    return 0;
}
