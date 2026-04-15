#include <parasail.h>
#include <any>
#include <cstdint>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <unordered_set>

#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"

#include "src/graph.h"
#include "src/ref.h"
#include "src/sv.h"
#include "src/utils.h"

using namespace std;

Logger main_logger("main");

void processGenome(const string& file, int color, uint32_t& sid, string& seq, LCPGraph& graph, const string& chr_str, int depth) {
    FastaStreamer streamer(file, false);
    for (const FastaPair& pair : streamer) {
        string chrom = pair.first;
        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
        if (chrom != chr_str) continue;

        sid = graph.seq_names.size();
        graph.seq_names.push_back(pair.first);
        seq = pair.second;

        struct lps lcp_string;
        init_lps(&lcp_string, seq.c_str(), seq.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, color, seq.c_str(), sid);
        free_lps(&lcp_string);
    }
}

int main(int argc, char** argv) {
    LCP_INIT();
    ArgParser parser("Build LCP graph from reference, normal, and tumor FASTA files.");
    parser.add("ref", "--ref", "-r", "Reference genome FASTA file", "string", true);
    parser.add("normal", "--normal", "-n", "Normal genome FASTA file", "string", true);
    parser.add("tumor", "--tumor", "-t", "Tumor genome FASTA file", "string", true);
    parser.add("depth", "--depth", "-d", "Depth of LCP", "int", true);
    parser.add("chr", "--chr", "-c", "Chromosome to process", "int", false);
    parser.add("output", "--output", "-o", "Output prefix (without extension)", "string", true);

    map<string, any> args = parser.parse(argc, argv);
    string ref_file = any_cast<string>(args["ref"]), normal_file = any_cast<string>(args["normal"]), tumor_file = any_cast<string>(args["tumor"]), output_prefix = any_cast<string>(args["output"]);
    int depth = any_cast<int>(args["depth"]);
    int chr = args.count("chr") ? any_cast<int>(args["chr"]) : 19;
    string chr_str = to_string(chr);

    namespace fs = std::filesystem;
    fs::path output_dir = "outputs";
    fs::path gfa_dir = output_dir / "gfa";
    fs::path vcf_dir = output_dir / "vcf";
    fs::create_directories(output_dir);
    fs::create_directories(gfa_dir);
    fs::create_directories(vcf_dir);

    RefIndex ref_index;
    LCPGraph graph;
    uint32_t normal_seq_id = 0, tumor_seq_id = 0;
    string ref_sequence, normal_sequence, tumor_sequence;

    // Parse reference into RefIndex + build ref hash set
    unordered_set<uint64_t> ref_hashes;
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
        for (int i = 0; i < lcp_string.size; ++i) {
            string region(ref_sequence.c_str() + lcp_string.cores[i].start,
                          lcp_string.cores[i].end - lcp_string.cores[i].start);
            ref_hashes.insert(hash_string(region, 0));
        }
        free_lps(&lcp_string);
    }

    // Build revcomp-ref hash set for large INV detection
    unordered_set<uint64_t> rcref_hashes;
    {
        string rcref = reverse_complement(ref_sequence);
        struct lps rcref_lcp;
        init_lps(&rcref_lcp, rcref.c_str(), rcref.size());
        lps_deepen(&rcref_lcp, depth);
        for (int i = 0; i < rcref_lcp.size; ++i) {
            string region(rcref.c_str() + rcref_lcp.cores[i].start,
                          rcref_lcp.cores[i].end - rcref_lcp.cores[i].start);
            rcref_hashes.insert(hash_string(region, 0));
        }
        free_lps(&rcref_lcp);
        main_logger.debug("Built rcref hash set: " + to_string(rcref_hashes.size()) + " unique hashes");
    }

    processGenome(normal_file, 0, normal_seq_id, normal_sequence, graph, chr_str, depth);  // Normal = bit 0 (0x1)
    processGenome(tumor_file, 1, tumor_seq_id, tumor_sequence, graph, chr_str, depth);     // Tumor = bit 1 (0x2)

    // Detect large INVs via revcomp-ref core index
    auto large_inv_calls = detectLargeINVs(tumor_sequence, depth, ref_hashes, rcref_hashes, ref_index);

    BubbleResult boundaries = graph.detectBubbles();
    auto calls = extractAndCallSVs(graph, ref_index, boundaries, large_inv_calls);

    // Add core-index INV calls (bubble-based calls inside these regions already skipped)
    calls.insert(calls.end(), large_inv_calls.begin(), large_inv_calls.end());

    // Call SVs and write VCF
    string vcf_file = (vcf_dir / (output_prefix + ".vcf")).string();
    ofstream vcf(vcf_file);
    vcf << "##fileformat=VCFv4.2\n"
        << "##source=svp\n"
        << "##contig=<ID=" << chr_str << ",length=" << ref_sequence.size() << ">\n"
        << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n"
        << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
        << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of variant\">\n"
        << "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Inserted sequence\">\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    main_logger.debug("Writing " + to_string(calls.size()) + " SV calls");
    constexpr int64_t INS_DUP_MAX_BP = 2000;
    constexpr double INS_DUP_SW_THRESHOLD = 0.70;
    constexpr int64_t INS_DUP_WINDOW_PAD = 300;

    int sw_match = 2, sw_mismatch = -4, sw_gap_open = 5, sw_gap_extend = 1;
    parasail_matrix_t* sw_matrix = parasail_matrix_create("ACGT", sw_match, sw_mismatch);

    int sv_count = 0;
    int ins_retyped_as_dup = 0;
    for (auto c : calls) {
        if (c.svtype == "INS" && c.length > 0 && c.length < INS_DUP_MAX_BP && !c.tumor_seq.empty()) {
            int64_t pad = std::max(INS_DUP_WINDOW_PAD, 2 * static_cast<int64_t>(c.tumor_seq.size()));
            int64_t win_start = std::max<int64_t>(0, c.ref_pos - pad);
            int64_t win_end = std::min<int64_t>(static_cast<int64_t>(ref_sequence.size()), c.ref_pos + pad);
            string window = ref_sequence.substr(static_cast<size_t>(win_start), static_cast<size_t>(win_end - win_start));

            parasail_result_t* sw_result = parasail_sw(
                c.tumor_seq.c_str(), static_cast<int>(c.tumor_seq.size()),
                window.c_str(), static_cast<int>(window.size()),
                sw_gap_open, sw_gap_extend, sw_matrix);

            int perfect_score = static_cast<int>(c.tumor_seq.size()) * sw_match;
            double normalized = (perfect_score > 0) ? static_cast<double>(sw_result->score) / perfect_score : 0.0;
            parasail_result_free(sw_result);

            if (normalized >= INS_DUP_SW_THRESHOLD) {
                c.svtype = "DUP";
                c.length = static_cast<int64_t>(c.tumor_seq.size());
                ins_retyped_as_dup++;
            }
        }

        sv_count++;
        string ref_allele, alt_allele, info;
        if (c.svtype == "DEL") {
            ref_allele = "N";
            alt_allele = "<DEL>";
            int64_t end_pos = c.ref_pos + c.length;
            info = "END=" + to_string(end_pos) + ";SVTYPE=DEL;SVLEN=" + to_string(c.length);
        } else if (c.svtype == "INS") {
            ref_allele = "N";
            alt_allele = "<INS>";
            info = "END=" + to_string(c.ref_pos) + ";SVTYPE=INS;SVLEN=" + to_string(c.length);
            if (!c.tumor_seq.empty()) info += ";INSSEQ=" + c.tumor_seq;
        } else if (c.svtype == "DUP") {
            ref_allele = "N";
            alt_allele = "<DUP>";
            int64_t interval_len = std::max<int64_t>(1, c.length >= 0 ? c.length : -c.length);
            int64_t end_pos = c.ref_pos + interval_len;
            info = "END=" + to_string(end_pos) + ";SVTYPE=DUP;SVLEN=" + to_string(c.length);
        } else if (c.svtype == "INV") {
            ref_allele = "N";
            alt_allele = "<INV>";
            int64_t end_pos = c.ref_pos + c.length;
            info = "END=" + to_string(end_pos) + ";SVTYPE=INV;SVLEN=" + to_string(c.length);
        } else {
            ref_allele = c.normal_seq.empty() ? "." : c.normal_seq;
            alt_allele = c.tumor_seq.empty() ? "." : c.tumor_seq;
            info = "SVTYPE=" + c.svtype + ";SVLEN=" + to_string(c.length);
        }
        vcf << chr_str << "\t" << c.ref_pos << "\t.\t" << ref_allele << "\t" << alt_allele << "\t.\tPASS\t" << info << "\n";
    }
    vcf.close();
    parasail_matrix_free(sw_matrix);
    main_logger.debug("Wrote " + to_string(sv_count) + " SVs (" + to_string(ins_retyped_as_dup) + " INS→DUP by local alignment) → " + vcf_file);

    string gfa_file = (gfa_dir / (output_prefix + "_graph.gfa")).string();
    main_logger.debug("Writing GFA to " + gfa_file);
    graph.writeGFA(gfa_file);

    main_logger.debug("Done!");
    return 0;
}
