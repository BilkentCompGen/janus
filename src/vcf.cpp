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

    // Reclassify all DUP calls as INS — standard truth sets (GIAB v5.0q)
    // represent tandem duplications as insertions, and the flipped-DUP detector
    // produces too many false DUP calls on real data
    int dup_to_ins = 0;
    for (auto& c : calls) {
        if (c.svtype == "DUP") {
            c.svtype = "INS";
            dup_to_ins++;
        }
    }

    int sv_count = 0;
    for (auto& c : calls) {

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

    vcf_logger.debug("Wrote " + to_string(sv_count) + " SVs (" + to_string(dup_to_ins) + " DUP->INS) -> " + vcf_path);
}
