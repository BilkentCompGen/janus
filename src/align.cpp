#include "align.h"

#include <parasail.h>

std::vector<VariantCall> alignAndCall(const std::string& n_seq, const std::string& t_seq, uint64_t r_start, const std::string& chrom) {
    std::vector<VariantCall> vars;
    if (n_seq.empty() || t_seq.empty()) return vars;

    auto* m = parasail_matrix_create("ACGTN", 2, -1);
    auto* r = parasail_nw_trace(n_seq.c_str(), n_seq.size(), t_seq.c_str(), t_seq.size(), 3, 1, m);
    auto* c = parasail_result_get_cigar(r, n_seq.c_str(), n_seq.size(), t_seq.c_str(), t_seq.size(), m);

    if (!c || c->len == 0) {
        if (c) parasail_cigar_free(c);
        parasail_result_free(r); parasail_matrix_free(m);
        return vars;
    }

    uint64_t r_pos = r_start;
    size_t n_idx = 0, t_idx = 0;

    for (int i = 0; i < c->len; ++i) {
        char op = parasail_cigar_decode_op(c->seq[i]);
        uint32_t len = parasail_cigar_decode_len(c->seq[i]);

        if (op == 'M' || op == '=' || op == 'X') {
            for (uint32_t j = 0; j < len; ++j) {
                if (n_seq[n_idx + j] != t_seq[t_idx + j]) {
                    vars.push_back({chrom, std::string(1, n_seq[n_idx+j]), std::string(1, t_seq[t_idx+j]), "SNP", r_pos + j, 0});
                }
            }
            r_pos += len; n_idx += len; t_idx += len;
        } else if (op == 'I') { // DEL in tumor
            if (len >= 1 && n_idx + len <= n_seq.size()) {
                uint64_t v_pos = (n_idx > 0) ? r_pos - 1 : r_pos;
                std::string v_ref = (n_idx > 0) ? n_seq.substr(n_idx - 1, 1) : "N";
                vars.push_back({chrom, v_ref, "<DEL>", "DEL", v_pos, -(int64_t)len});
            }
            r_pos += len; n_idx += len;
        } else if (op == 'D') { // INS in tumor
            if (len >= 1 && t_idx + len <= t_seq.size()) {
                uint64_t v_pos = r_pos;
                std::string v_ref = (n_idx < n_seq.size()) ? n_seq.substr(n_idx, 1) : "N";
                vars.push_back({chrom, v_ref, "<INS>", "INS", v_pos, (int64_t)len});
            }
            t_idx += len;
        }
    }

    parasail_cigar_free(c); parasail_result_free(r); parasail_matrix_free(m);

    // Filter artifact SNPs nearby INDELs
    std::vector<VariantCall> filtered_vars;
    for (size_t i = 0; i < vars.size(); ++i) {
        bool is_artifact = false;
        if (vars[i].svtype == "SNP") {
            for (size_t j = 0; j < vars.size(); ++j) {
                if (i == j) continue;
                if (vars[j].svtype == "INS" || vars[j].svtype == "DEL") {
                    long indel_start = vars[j].pos;
                    long indel_end = vars[j].pos;
                    if (vars[j].svtype == "DEL") {
                        indel_end = vars[j].pos + std::abs(vars[j].svlen);
                    }
                    
                    long dist_s = std::abs((long)vars[i].pos - indel_start);
                    long dist_e = std::abs((long)vars[i].pos - indel_end);
                    
                    if (dist_s <= 5 || dist_e <= 5) { // 5bp window from either end
                        is_artifact = true;
                        break;
                    }
                }
            }
        }
        if (!is_artifact) filtered_vars.push_back(vars[i]);
    }
    return filtered_vars;
}
