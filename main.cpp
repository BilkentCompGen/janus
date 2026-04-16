#include <any>
#include <filesystem>
#include <map>
#include <string>

#include "lib/lcptools/lps.h"

#include "src/diag.h"
#include "src/graph.h"
#include "src/pipeline.h"
#include "src/ref.h"
#include "src/sv.h"
#include "src/utils.h"
#include "src/vcf.h"

using namespace std;

int main(int argc, char** argv) {
    LCP_INIT();
    ArgParser parser("Build LCP graph from reference, normal, and tumor FASTA files.");
    parser.add("ref", "--ref", "-r", "Reference genome FASTA file", "string", true);
    parser.add("normal", "--normal", "-n", "Normal genome FASTA file", "string", true);
    parser.add("tumor", "--tumor", "-t", "Tumor genome FASTA file", "string", true);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("chr", "--chr", "-c", "Chromosome to process", "int", false);
    parser.add("output", "--output", "-o", "Output prefix (without extension)", "string", true);
    parser.add("debug", "--debug", "-D", "Enable diagnostic output (writes TSV files)", "bool", false);

    map<string, any> args = parser.parse(argc, argv);
    string ref_file = any_cast<string>(args["ref"]);
    string normal_file = any_cast<string>(args["normal"]);
    string tumor_file = any_cast<string>(args["tumor"]);
    string output_prefix = any_cast<string>(args["output"]);
    int depth = any_cast<int>(args["depth"]);
    int chr = args.count("chr") ? any_cast<int>(args["chr"]) : 19;
    string chr_str = to_string(chr);
    bool debug_mode = args.count("debug") > 0;

    namespace fs = std::filesystem;
    fs::path output_dir = "outputs";
    fs::path gfa_dir = output_dir / "gfa";
    fs::path vcf_dir = output_dir / "vcf";
    fs::create_directories(gfa_dir);
    fs::create_directories(vcf_dir);

    if (debug_mode) {
        string diag_prefix = (vcf_dir / output_prefix).string();
        DiagLog::instance().enable(diag_prefix);
    }

    RefIndex ref_index;
    LCPGraph graph;
    uint32_t normal_seq_id = 0, tumor_seq_id = 0;
    string ref_sequence, normal_sequence, tumor_sequence;

    auto ref_hashes = buildRefHashes(ref_file, chr_str, depth, ref_index, ref_sequence);
    auto rcref_hashes = buildRcrefHashes(ref_sequence, depth);

    processGenome(normal_file, 0, normal_seq_id, normal_sequence, graph, chr_str, depth);
    processGenome(tumor_file, 1, tumor_seq_id, tumor_sequence, graph, chr_str, depth);

    auto large_inv_calls = detectLargeINVs(tumor_sequence, depth, ref_hashes, rcref_hashes, ref_index);

    auto tumor_hashes = buildTumorHashes(tumor_sequence, depth);
    auto large_del_calls = detectLargeDELs(ref_index, tumor_hashes, large_inv_calls);
    auto large_ins_calls = detectLargeINS(tumor_sequence, depth, ref_hashes, rcref_hashes, ref_index);

    BubbleResult boundaries = graph.detectBubbles();
    auto calls = extractAndCallSVs(graph, ref_index, boundaries, large_inv_calls);
    calls.insert(calls.end(), large_inv_calls.begin(), large_inv_calls.end());

    // Deduplicate core-index DEL/INS calls against bubble calls:
    // skip if a bubble call of the same type overlaps within 1kb
    auto dedup = [&calls](const vector<SVCall>& new_calls) {
        int added = 0;
        for (const auto& nc : new_calls) {
            bool duplicate = false;
            for (const auto& ec : calls) {
                if (ec.svtype != nc.svtype) continue;
                // Check if positions overlap within tolerance
                int64_t ec_end = ec.ref_pos + ec.length;
                int64_t nc_end = nc.ref_pos + nc.length;
                if (nc.ref_pos < ec_end + 1000 && nc_end > ec.ref_pos - 1000) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                calls.push_back(nc);
                added++;
            }
        }
        return added;
    };
    int added_dels = dedup(large_del_calls);
    int added_ins = dedup(large_ins_calls);

    if (debug_mode) {
        DiagLog::instance().addStat("large_del_after_dedup", to_string(added_dels));
        DiagLog::instance().addStat("large_ins_after_dedup", to_string(added_ins));
    }

    if (debug_mode) {
        auto& diag = DiagLog::instance();
        diag.addStat("graph_nodes", to_string(graph.nodeCount()));
        diag.addStat("graph_edges", to_string(graph.next_edges.size()));
        diag.addStat("ref_cores", to_string(ref_index.cores.size()));
        diag.addStat("bubble_starts", to_string(boundaries.bubble_starts.size()));
        diag.addStat("bubble_ends", to_string(boundaries.bubble_ends.size()));

        // Count node colors
        int shared = 0, normal_only = 0, tumor_only = 0;
        for (const auto& [h, n] : graph.nodes) {
            bool has_n = (n.colors & 0x1) != 0;
            bool has_t = (n.colors & 0x2) != 0;
            if (has_n && has_t) shared++;
            else if (has_n) normal_only++;
            else if (has_t) tumor_only++;
        }
        diag.addStat("nodes_shared", to_string(shared));
        diag.addStat("nodes_normal_only", to_string(normal_only));
        diag.addStat("nodes_tumor_only", to_string(tumor_only));

        // Count unique vs non-unique ref cores
        int unique_cores = 0, ambiguous_cores = 0;
        for (const auto& [h, indices] : ref_index.index) {
            if (indices.size() == 1) unique_cores++;
            else ambiguous_cores++;
        }
        diag.addStat("ref_unique_hashes", to_string(unique_cores));
        diag.addStat("ref_ambiguous_hashes", to_string(ambiguous_cores));
    }

    string vcf_file = (vcf_dir / (output_prefix + ".vcf")).string();
    writeVCF(vcf_file, calls, ref_sequence, chr_str);

    string gfa_file = (gfa_dir / (output_prefix + "_graph.gfa")).string();
    graph.writeGFA(gfa_file);

    if (debug_mode) {
        DiagLog::instance().flush();
    }

    return 0;
}
