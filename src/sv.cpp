#include <parasail.h>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <unordered_set>

#include "ref.h"
#include "utils.h"
#include "sv.h"

using namespace std;

static Logger sv_logger("sv");

namespace {

bool isNormalOnly(const Node& n) {
    return (n.colors & 0x1) && !(n.colors & 0x2);
}

bool isTumorOnly(const Node& n) {
    return (n.colors & 0x2) && !(n.colors & 0x1);
}

bool isShared(const Node& n) {
    return (n.colors & 0x3) == 0x3;
}

struct EdgeColorCounts { int normal_only = 0, tumor_only = 0, shared = 0; };

EdgeColorCounts countSuccessorColors(const LCPGraph& graph, node_hash h) {
    EdgeColorCounts c;
    auto it = graph.next_edges.find(h);
    if (it == graph.next_edges.end()) return c;
    for (node_hash nh : it->second) {
        auto nit = graph.nodes.find(nh);
        if (nit == graph.nodes.end()) continue;
        if (isNormalOnly(nit->second)) c.normal_only++;
        else if (isTumorOnly(nit->second)) c.tumor_only++;
        else if (isShared(nit->second)) c.shared++;
    }
    return c;
}

EdgeColorCounts countPredecessorColors(const LCPGraph& graph, node_hash h) {
    EdgeColorCounts c;
    auto it = graph.prev_edges.find(h);
    if (it == graph.prev_edges.end()) return c;
    for (node_hash ph : it->second) {
        auto pit = graph.nodes.find(ph);
        if (pit == graph.nodes.end()) continue;
        if (isNormalOnly(pit->second)) c.normal_only++;
        else if (isTumorOnly(pit->second)) c.tumor_only++;
        else if (isShared(pit->second)) c.shared++;
    }
    return c;
}

struct AnchoredCandidate {
    node_hash hash;
    int64_t ref_start, ref_end;
};

vector<SVCall> detectFlippedDUPs(
    const LCPGraph& graph,
    const RefIndex& ref_index,
    const BubbleResult& markers,
    const unordered_set<string>& bubble_pairs)
{
    // Step 1: Filter bubble starts/ends to those with unique ref anchors
    // Step 2: Classify by T=1 S>=1 N=0 pattern
    vector<AnchoredCandidate> dup_starts, dup_ends;

    for (node_hash h : markers.bubble_starts) {
        int64_t rs = -1, re = -1;
        if (!ref_index.anchorSpan(h, rs, re)) continue;
        EdgeColorCounts c = countSuccessorColors(graph, h);
        if (c.tumor_only >= 1 && c.shared >= 1 && c.normal_only == 0)
            dup_starts.push_back({h, rs, re});
    }

    for (node_hash h : markers.bubble_ends) {
        int64_t rs = -1, re = -1;
        if (!ref_index.anchorSpan(h, rs, re)) continue;
        EdgeColorCounts c = countPredecessorColors(graph, h);
        if (c.tumor_only >= 1 && c.shared >= 1 && c.normal_only == 0)
            dup_ends.push_back({h, rs, re});
    }

    // Step 3: Merge starts and ends into one sorted event list, pair E→S in order
    struct Event {
        int64_t ref_start;
        node_hash hash;
        bool is_start; // true = bubble start (near GT DUP end), false = bubble end (near GT DUP start)
    };
    vector<Event> events;
    for (const auto& ds : dup_starts) events.push_back({ds.ref_start, ds.hash, true});
    for (const auto& de : dup_ends) events.push_back({de.ref_start, de.hash, false});
    sort(events.begin(), events.end(), [](const Event& a, const Event& b) { return a.ref_start < b.ref_start; });

    vector<SVCall> calls;
    int64_t pending_end_pos = -1;
    node_hash pending_end_hash = 0;

    for (const auto& ev : events) {
        if (!ev.is_start) {
            // Bubble end (near DUP start) — remember it
            pending_end_pos = ev.ref_start;
            pending_end_hash = ev.hash;
        } else if (pending_end_pos >= 0) {
            // Bubble start (near DUP end) — pair with pending end
            string pair_key = to_string(ev.hash) + ":" + to_string(pending_end_hash);
            if (!bubble_pairs.count(pair_key)) {
                int64_t pos = pending_end_pos;
                int64_t svlen = ev.ref_start - pending_end_pos;
                if (svlen > 0) {
                    calls.emplace_back("", "", "DUP", 0, svlen, pos);
                }
            }
            pending_end_pos = -1;
            pending_end_hash = 0;
        }
    }

    sv_logger.debug("detectFlippedDUPs: candidate_starts=" + to_string(dup_starts.size())
        + " candidate_ends=" + to_string(dup_ends.size())
        + " paired=" + to_string(calls.size()));

    return calls;
}

} // namespace

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

vector<SVCall> extractAndCallSVs(const LCPGraph& graph, const RefIndex& ref_index, const BubbleResult& markers) {
    vector<SVCall> calls;
    vector<Bubble> bubbles = graph.extractBubbles(markers);

    unordered_set<string> bubble_pairs;
    for (const auto& b : bubbles) {
        bubble_pairs.insert(to_string(b.start_node) + ":" + to_string(b.end_node));
    }

    int anchored_bubbles = 0;
    for (const auto& b : bubbles) {
        int64_t ref_offset = ref_index.anchorBubble(b, graph);
        if (ref_offset < 0) continue;
        anchored_bubbles++;
        vector<SVCall> bubble_calls = callSV(b, ref_offset);
        calls.insert(calls.end(), bubble_calls.begin(), bubble_calls.end());
    }

    vector<SVCall> dup_calls = detectFlippedDUPs(graph, ref_index, markers, bubble_pairs);
    calls.insert(calls.end(), dup_calls.begin(), dup_calls.end());

    sv_logger.debug(
        "extractAndCallSVs: bubbles=" + to_string(bubbles.size())
        + " anchored_bubbles=" + to_string(anchored_bubbles)
        + " flipped_dups=" + to_string(dup_calls.size())
        + " total_calls=" + to_string(calls.size())
    );

    return calls;
}
