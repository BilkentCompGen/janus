#include <parasail.h>
#include <algorithm>
#include <cstdint>

#include "utils.h"
#include "sv.h"

using namespace std;

static Logger sv_logger("sv");

SVCall::SVCall(const string& normal_seq, const string& tumor_seq, const string& svtype, int64_t start, int64_t length, int64_t ref_pos) : normal_seq(normal_seq), tumor_seq(tumor_seq), svtype(svtype), start(start), length(length), ref_pos(ref_pos) {}

SVRegion detectSV(const string& sequence, const string& query) { // sequence is larger than query
    int i = 0, j = sequence.size() - 1, k = query.size() - 1;
    while (sequence[i] == query[i]) i++;
    while (sequence[j] == query[k]) { j--; k--; }

    if (i > j) return {-1, 0, "no_call"};
    if (i == k && i == j) return {i, 1, "snv"};
    if (k < i) return {i, j + 1 - i, "indel"};
    return {i, j + 1 - i, "multi"};
}

vector<SVCall> callMulti(const string& normal_seq, const string& tumor_seq, const SVRegion& /* region */, int64_t ref_offset) {
    int match = 5, mismatch = -4, gap_open = 50, gap_extend = 1;
    parasail_matrix_t* matrix = parasail_matrix_create("ACGT", match, mismatch);

    // tumor = query (s1), normal = reference (s2)
    // CIGAR: I = bases in tumor not in normal, D = bases in normal not in tumor
    parasail_result_t* result = parasail_nw_trace(tumor_seq.c_str(), (int)tumor_seq.size(), normal_seq.c_str(), (int)normal_seq.size(), gap_open, gap_extend, matrix);
    parasail_cigar_t* cigar = parasail_result_get_cigar(result, tumor_seq.c_str(), (int)tumor_seq.size(), normal_seq.c_str(), (int)normal_seq.size(), matrix);

    // Parse CIGAR into (op, len) pairs
    struct CigarOp { char op; int64_t len; };
    vector<CigarOp> ops;
    for (int i = 0; i < cigar->len; i++) {
        ops.push_back({parasail_cigar_decode_op(cigar->seq[i]), parasail_cigar_decode_len(cigar->seq[i])});
    }

    parasail_cigar_free(cigar);
    parasail_result_free(result);
    parasail_matrix_free(matrix);

    // Group CIGAR operations into variant blocks separated by long match anchors
    // (>= MIN_ANCHOR bp). Short matches between variant ops are absorbed into the
    // block, preventing fragmentation of complex indel regions.
    //
    // Strategy:
    //  - Blocks WITH indels (D/I ops): consolidate into a single call per block.
    //    This fixes fragmentation where nearby DEL+INS produce hundreds of tiny
    //    alternating ops instead of a few clean calls.
    //  - Blocks WITHOUT indels (only X mismatches + short matches): emit each
    //    mismatch position as an individual 1bp SNP call, preserving SNP resolution.
    const int64_t MIN_ANCHOR = 50;

    vector<SVCall> calls;
    int64_t pos = 0, query_pos = 0;

    // Emit a single consolidated SVCall from an indel-containing block
    auto emitIndelBlock = [&](int64_t ref_start, int64_t qry_start, int64_t ref_len, int64_t qry_len) {
        if (ref_len == 0 && qry_len == 0) return;

        string ref_allele = ref_len > 0 ? normal_seq.substr(ref_start, ref_len) : "";
        string alt_allele = qry_len > 0 ? tumor_seq.substr(qry_start, qry_len) : "";

        string svtype;
        int64_t sv_len;

        if (qry_len == 0) {
            svtype = "DEL"; sv_len = ref_len;
        } else if (ref_len == 0) {
            svtype = "INS"; sv_len = qry_len;
        } else if (ref_len > qry_len) {
            svtype = "DEL"; sv_len = ref_len;
        } else if (qry_len > ref_len) {
            svtype = "INS"; sv_len = qry_len;
        } else {
            svtype = "SNP"; sv_len = ref_len;
        }

        calls.emplace_back(ref_allele, alt_allele, svtype, ref_start, sv_len, ref_offset + ref_start);
    };

    // Emit individual 1bp SNP calls from a mismatch-only block (no indels)
    auto emitSnpBlock = [&](size_t start_idx, size_t end_idx, int64_t ref_start, int64_t qry_start) {
        int64_t p = ref_start, qp = qry_start;
        for (size_t j = start_idx; j < end_idx; j++) {
            if (ops[j].op == 'X') {
                for (int64_t k = 0; k < ops[j].len; k++) {
                    calls.emplace_back(
                        string(1, normal_seq[p + k]),
                        string(1, tumor_seq[qp + k]),
                        "SNP", p + k, 1, ref_offset + p + k);
                }
            }
            if (ops[j].op == '=' || ops[j].op == 'M' || ops[j].op == 'X') {
                p += ops[j].len;
                qp += ops[j].len;
            }
        }
    };

    size_t i = 0;
    while (i < ops.size()) {
        char op = ops[i].op;
        int64_t len = ops[i].len;
        bool is_anchor = (op == '=' || op == 'M') && len >= MIN_ANCHOR;

        if (is_anchor) {
            pos += len;
            query_pos += len;
            i++;
            continue;
        }

        // Collect a variant block: all ops until the next long match anchor
        size_t block_start_idx = i;
        int64_t block_ref_start = pos, block_qry_start = query_pos;
        int64_t block_ref_len = 0, block_qry_len = 0;
        bool has_indel = false;

        while (i < ops.size()) {
            op = ops[i].op;
            len = ops[i].len;
            if ((op == '=' || op == 'M') && len >= MIN_ANCHOR) break;

            if (op == 'D' || op == 'I') has_indel = true;

            if (op == '=' || op == 'M' || op == 'X') {
                block_ref_len += len;
                block_qry_len += len;
            } else if (op == 'D') {
                block_ref_len += len;
            } else if (op == 'I') {
                block_qry_len += len;
            }
            i++;
        }

        if (has_indel) {
            emitIndelBlock(block_ref_start, block_qry_start, block_ref_len, block_qry_len);
        } else {
            emitSnpBlock(block_start_idx, i, block_ref_start, block_qry_start);
        }

        pos = block_ref_start + block_ref_len;
        query_pos = block_qry_start + block_qry_len;
    }

    sv_logger.debug("callMulti: " + to_string(ops.size()) + " CIGAR ops -> " + to_string(calls.size()) + " calls (anchor=" + to_string(MIN_ANCHOR) + "bp)");
    return calls;
}

vector<SVCall> callSV(const Bubble& bubble, int64_t ref_offset) {
    const string& normal_seq = bubble.normal_seq;
    const string& tumor_seq = bubble.tumor_seq;

    SVRegion region = (normal_seq.size() >= tumor_seq.size()) ? detectSV(normal_seq, tumor_seq) : detectSV(tumor_seq, normal_seq);
    if (region.start == -1) return {};

    if (region.call == "snv") return { SVCall(string(1, normal_seq[region.start]), string(1, tumor_seq[region.start]), "SNP", region.start, region.length, ref_offset + region.start) };
    if (region.call == "indel") {
        if (normal_seq.size() > tumor_seq.size()) return { SVCall(normal_seq.substr(region.start, region.length), "", "DEL", region.start, region.length, ref_offset + region.start) };
        else return { SVCall("", tumor_seq.substr(region.start, region.length), "INS", region.start, region.length, ref_offset + region.start) };
    }
    if (region.call == "multi") return callMulti(normal_seq, tumor_seq, region, ref_offset);
    return {};
}
