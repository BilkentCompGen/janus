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
            r_pos += len; n_idx += len; t_idx += len;
        } else if (op == 'I') {
            if (n_idx + len <= n_seq.size()) {
                VariantCall v = {chrom, "", "", "DEL", r_pos, -(int64_t)len};
                if (n_idx > 0) {
                    v.pos = r_pos - 1;
                    v.ref_allele = n_seq.substr(n_idx - 1, len + 1);
                    v.alt_allele = n_seq.substr(n_idx - 1, 1);
                } else {
                    v.ref_allele = n_seq.substr(n_idx, len); v.alt_allele = ".";
                }
                vars.push_back(v);
            }
            r_pos += len; n_idx += len;
        } else if (op == 'D') {
            if (t_idx + len <= t_seq.size()) {
                VariantCall v = {chrom, "", "", "INS", r_pos, (int64_t)len};
                if (n_idx > 0 && n_idx <= n_seq.size()) {
                    v.pos = r_pos - 1;
                    v.ref_allele = n_seq.substr(n_idx - 1, 1);
                    v.alt_allele = n_seq.substr(n_idx - 1, 1) + t_seq.substr(t_idx, len);
                } else {
                    v.ref_allele = "N"; v.alt_allele = "N" + t_seq.substr(t_idx, len);
                }
                vars.push_back(v);
            }
            t_idx += len;
        }
    }

    parasail_cigar_free(c); parasail_result_free(r); parasail_matrix_free(m);
    return vars;
}
