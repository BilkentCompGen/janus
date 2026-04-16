#include <parasail.h>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <string>

#include "vcf.h"
#include "utils.h"

using namespace std;

static Logger vcf_logger("vcf");

void writeVCF(const string& vcf_path, vector<SVCall>& calls,
              const string& ref_sequence, const string& chr_str) {
    ofstream vcf(vcf_path);
    vcf << "##fileformat=VCFv4.2\n"
        << "##source=svp\n"
        << "##contig=<ID=" << chr_str << ",length=" << ref_sequence.size() << ">\n"
        << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n"
        << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
        << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of variant\">\n"
        << "##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Inserted sequence\">\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    vcf_logger.debug("Writing " + to_string(calls.size()) + " SV calls");

    constexpr int64_t INS_DUP_MAX_BP = 2000;
    constexpr double INS_DUP_SW_THRESHOLD = 0.70;
    constexpr int64_t INS_DUP_WINDOW_PAD = 300;

    int sw_match = 2, sw_mismatch = -4, sw_gap_open = 5, sw_gap_extend = 1;
    parasail_matrix_t* sw_matrix = parasail_matrix_create("ACGT", sw_match, sw_mismatch);

    int sv_count = 0;
    int ins_retyped_as_dup = 0;
    for (auto& c : calls) {
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

    vcf_logger.debug("Wrote " + to_string(sv_count) + " SVs (" + to_string(ins_retyped_as_dup) + " INS->DUP by local alignment) -> " + vcf_path);
}
