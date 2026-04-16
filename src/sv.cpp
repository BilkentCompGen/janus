#include <parasail.h>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <fstream>
#include <unordered_set>

#include "diag.h"
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

// Post-process callMulti results: reclassify SNP clusters as INV where
// reverse-complement alignment confirms an inversion.
// Handles both multi-SNP clusters (small INVs) and single large SNP blocks (medium INVs).
static void reclassifySnpAsInv(vector<SVCall>& calls,
    const string& normal_seq, const string& tumor_seq, int64_t ref_offset) {

    constexpr int64_t MAX_GAP = 10;
    constexpr int64_t MIN_SPAN = 10;
    constexpr double INV_THRESHOLD = 0.65;

    // Collect SNP indices sorted by ref_pos
    vector<size_t> snp_idx;
    for (size_t i = 0; i < calls.size(); i++)
        if (calls[i].svtype == "SNP") snp_idx.push_back(i);
    if (snp_idx.empty()) return;

    sort(snp_idx.begin(), snp_idx.end(),
         [&](size_t a, size_t b) { return calls[a].ref_pos < calls[b].ref_pos; });

    // Build clusters (gap <= MAX_GAP between consecutive SNP calls)
    struct Cluster { vector<size_t> indices; int64_t ref_start, ref_end; };
    vector<Cluster> clusters;
    Cluster cur;
    cur.indices.push_back(snp_idx[0]);
    cur.ref_start = calls[snp_idx[0]].ref_pos;
    cur.ref_end = calls[snp_idx[0]].ref_pos + calls[snp_idx[0]].length;

    for (size_t i = 1; i < snp_idx.size(); i++) {
        auto& c = calls[snp_idx[i]];
        if (c.ref_pos - cur.ref_end <= MAX_GAP) {
            cur.indices.push_back(snp_idx[i]);
            cur.ref_end = max(cur.ref_end, c.ref_pos + c.length);
        } else {
            if (cur.ref_end - cur.ref_start >= MIN_SPAN) clusters.push_back(cur);
            cur = {};
            cur.indices.push_back(snp_idx[i]);
            cur.ref_start = c.ref_pos;
            cur.ref_end = c.ref_pos + c.length;
        }
    }
    if (cur.ref_end - cur.ref_start >= MIN_SPAN) clusters.push_back(cur);

    // Verify each cluster with revcomp alignment
    unordered_set<size_t> remove_set;
    vector<SVCall> inv_calls;

    int sw_match = 2, sw_mismatch = -4, sw_gap_open = 5, sw_gap_extend = 1;
    parasail_matrix_t* matrix = parasail_matrix_create("ACGT", sw_match, sw_mismatch);

    for (auto& cl : clusters) {
        int64_t local_start = cl.ref_start - ref_offset;
        int64_t local_end = cl.ref_end - ref_offset;
        if (local_start < 0 || local_end > (int64_t)normal_seq.size()
            || local_end > (int64_t)tumor_seq.size()) continue;

        string norm_region = normal_seq.substr(local_start, local_end - local_start);
        string tum_region = tumor_seq.substr(local_start, local_end - local_start);
        if (norm_region.empty() || tum_region.empty()) continue;

        string rc = reverse_complement(norm_region);
        parasail_result_t* result = parasail_sw(
            tum_region.c_str(), (int)tum_region.size(),
            rc.c_str(), (int)rc.size(),
            sw_gap_open, sw_gap_extend, matrix);

        int perfect = (int)min(tum_region.size(), rc.size()) * sw_match;
        double score = (perfect > 0) ? (double)result->score / perfect : 0.0;
        parasail_result_free(result);

        if (score >= INV_THRESHOLD) {
            for (size_t idx : cl.indices) remove_set.insert(idx);
            int64_t span = cl.ref_end - cl.ref_start;
            inv_calls.emplace_back(norm_region, tum_region, "INV",
                                   local_start, span, cl.ref_start);
            sv_logger.debug("SNP cluster → INV: pos=" + to_string(cl.ref_start)
                + " len=" + to_string(span) + " sw=" + to_string(score)
                + " snps=" + to_string(cl.indices.size()));
        }
    }

    parasail_matrix_free(matrix);

    if (!remove_set.empty()) {
        vector<SVCall> filtered;
        for (size_t i = 0; i < calls.size(); i++)
            if (!remove_set.count(i)) filtered.push_back(calls[i]);
        filtered.insert(filtered.end(), inv_calls.begin(), inv_calls.end());
        calls = filtered;
    }
}

vector<SVCall> callSV(const Bubble& bubble, int64_t ref_offset) {
    const string& normal_seq = bubble.normal_seq;
    const string& tumor_seq = bubble.tumor_seq;

    constexpr int64_t MAX_ALIGN_SIZE = 5000;

    SVRegion region = (normal_seq.size() >= tumor_seq.size()) ? detectSV(normal_seq, tumor_seq) : detectSV(tumor_seq, normal_seq);
    if (region.start == -1) return {};

    if (region.call == "snv") return { SVCall(string(1, normal_seq[region.start]), string(1, tumor_seq[region.start]), "SNP", region.start, region.length, ref_offset + region.start) };
    if (region.call == "indel") {
        if (normal_seq.size() > tumor_seq.size()) return { SVCall(normal_seq.substr(region.start, region.length), "", "DEL", region.start, region.length, ref_offset + region.start) };
        else return { SVCall("", tumor_seq.substr(region.start, region.length), "INS", region.start, region.length, ref_offset + region.start) };
    }
    if (region.call == "multi") {
        // Skip parasail alignment for huge bubbles — NW on >5kb uses O(N*M) memory
        if ((int64_t)normal_seq.size() <= MAX_ALIGN_SIZE && (int64_t)tumor_seq.size() <= MAX_ALIGN_SIZE) {
            auto results = callMulti(normal_seq, tumor_seq, region, ref_offset);
            reclassifySnpAsInv(results, normal_seq, tumor_seq, ref_offset);
            return results;
        }
        return {};
    }
    return {};
}

vector<SVCall> extractAndCallSVs(const LCPGraph& graph, const RefIndex& ref_index, const BubbleResult& markers,
    const vector<SVCall>& large_inv_regions) {
    vector<SVCall> calls;
    vector<Bubble> bubbles = graph.extractBubbles(markers);

    unordered_set<string> bubble_pairs;
    for (const auto& b : bubbles) {
        bubble_pairs.insert(to_string(b.start_node) + ":" + to_string(b.end_node));
    }

    int anchored_bubbles = 0;
    int skipped_in_inv = 0;
    int unanchored_bubbles = 0;
    auto& diag = DiagLog::instance();

    for (const auto& b : bubbles) {
        int64_t ref_offset = ref_index.anchorBubble(b, graph);

        if (ref_offset < 0) {
            unanchored_bubbles++;
            if (diag.enabled()) {
                DiagBubble db;
                db.start_hash = b.start_node; db.end_hash = b.end_node;
                db.anchor_pos = -1; db.anchor_method = "failed";
                db.normal_seq_len = b.normal_seq.size();
                db.tumor_seq_len = b.tumor_seq.size();
                db.skip_reason = "anchor_failed";
                db.calls_emitted = 0;
                diag.addBubble(db);
            }
            continue;
        }
        anchored_bubbles++;

        // Skip bubbles inside core-index detected large INV regions
        // Padding accounts for core boundary offset between tumor and reference coordinates
        constexpr int64_t INV_SKIP_PAD = 500;
        bool inside_large_inv = false;
        for (const auto& li : large_inv_regions) {
            if (ref_offset >= li.ref_pos - INV_SKIP_PAD && ref_offset <= li.ref_pos + li.length + INV_SKIP_PAD) {
                inside_large_inv = true;
                break;
            }
        }
        // Also skip bubbles where both sequences are large — likely undetected INV
        // or merged variant that would blow up parasail O(N*M)
        constexpr int64_t MAX_BUBBLE_SEQ = 5000;
        bool too_large = false;
        if (!inside_large_inv
            && (int64_t)b.normal_seq.size() > MAX_BUBBLE_SEQ
            && (int64_t)b.tumor_seq.size() > MAX_BUBBLE_SEQ) {
            inside_large_inv = true;
            too_large = true;
        }

        if (inside_large_inv) skipped_in_inv++;

        vector<SVCall> bubble_calls = callSV(b, ref_offset);

        if (diag.enabled()) {
            DiagBubble db;
            db.start_hash = b.start_node; db.end_hash = b.end_node;
            db.anchor_pos = ref_offset; db.anchor_method = "ok";
            db.normal_seq_len = b.normal_seq.size();
            db.tumor_seq_len = b.tumor_seq.size();
            db.skip_reason = inside_large_inv ? (too_large ? "too_large_>5kb" : "inside_large_inv") : "";
            db.calls_emitted = (int)bubble_calls.size();
            string types;
            for (const auto& c : bubble_calls) {
                if (!types.empty()) types += ",";
                types += c.svtype;
            }
            db.call_types = types;
            diag.addBubble(db);

            for (const auto& c : bubble_calls) {
                diag.addCall({c.svtype, c.ref_pos, c.length, "bubble",
                    "nlen=" + to_string(b.normal_seq.size()) + " tlen=" + to_string(b.tumor_seq.size())});
            }
        }

        calls.insert(calls.end(), bubble_calls.begin(), bubble_calls.end());
    }

    vector<SVCall> dup_calls = detectFlippedDUPs(graph, ref_index, markers, bubble_pairs);
    calls.insert(calls.end(), dup_calls.begin(), dup_calls.end());

    sv_logger.debug(
        "extractAndCallSVs: bubbles=" + to_string(bubbles.size())
        + " anchored=" + to_string(anchored_bubbles)
        + " unanchored=" + to_string(unanchored_bubbles)
        + " skipped_in_inv=" + to_string(skipped_in_inv)
        + " flipped_dups=" + to_string(dup_calls.size())
        + " total_calls=" + to_string(calls.size())
    );

    if (diag.enabled()) {
        diag.addStat("total_bubbles", to_string(bubbles.size()));
        diag.addStat("anchored_bubbles", to_string(anchored_bubbles));
        diag.addStat("unanchored_bubbles", to_string(unanchored_bubbles));
        diag.addStat("skipped_in_inv", to_string(skipped_in_inv));
        diag.addStat("flipped_dups", to_string(dup_calls.size()));
        diag.addStat("total_calls", to_string(calls.size()));

        for (const auto& dc : dup_calls) {
            diag.addCall({dc.svtype, dc.ref_pos, dc.length, "flipped_dup", "."});
        }
    }

    return calls;
}

vector<SVCall> detectLargeINVs(const string& tumor_seq, int depth,
    const unordered_set<uint64_t>& ref_hashes,
    const unordered_set<uint64_t>& rcref_hashes,
    const RefIndex& ref_index) {

    // Extract tumor LCP cores
    struct lps lcp;
    init_lps(&lcp, tumor_seq.c_str(), tumor_seq.size());
    lps_deepen(&lcp, depth);

    // Classify each tumor core
    // 0=neither, 1=ref_only, 2=rcref_only, 3=both
    struct ClassifiedCore { uint64_t start, end; uint8_t cls; node_hash hash; };
    vector<ClassifiedCore> cores(lcp.size);
    for (int i = 0; i < lcp.size; ++i) {
        string region(tumor_seq.c_str() + lcp.cores[i].start, lcp.cores[i].end - lcp.cores[i].start);
        uint64_t h = hash_string(region, 0);
        bool in_ref = ref_hashes.count(h) > 0;
        bool in_rcref = rcref_hashes.count(h) > 0;
        uint8_t cls = 0;
        if (in_ref && in_rcref) cls = 3;
        else if (in_ref) cls = 1;
        else if (in_rcref) cls = 2;
        cores[i] = {lcp.cores[i].start, lcp.cores[i].end, cls, h};
    }
    free_lps(&lcp);

    // Find runs of rcref-only cores (gap-tolerant)
    constexpr int MAX_GAP = 3;
    constexpr int MIN_RCREF_CORES = 5;
    constexpr int64_t MIN_INV_BP = 500;

    // Collect raw runs of rcref-only cores
    struct InvRun { size_t core_start, core_end; int rcref_count; };
    vector<InvRun> runs;

    size_t i = 0;
    while (i < cores.size()) {
        if (cores[i].cls != 2) { i++; continue; }

        size_t run_start = i;
        int rcref_count = 0;
        int gap = 0;
        size_t run_end = i;

        while (i < cores.size()) {
            if (cores[i].cls == 2) {
                rcref_count++;
                gap = 0;
                run_end = i;
                i++;
            } else {
                gap++;
                if (gap > MAX_GAP) break;
                i++;
            }
        }

        if (rcref_count >= MIN_RCREF_CORES) {
            runs.push_back({run_start, run_end, rcref_count});
        }
    }

    // Merge adjacent runs only if gap cores are NOT ref-only (i.e. the gap is
    // ambiguous "both" or "neither" cores, not clearly forward-ref sequence)
    // If gap has any ref_only cores, it's a real break between separate inversions.
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t r = 0; r + 1 < runs.size(); r++) {
            size_t gap_start = runs[r].core_end + 1;
            size_t gap_end = runs[r+1].core_start;
            if (gap_start >= gap_end) {
                // Overlapping or adjacent runs — always merge
                runs[r].core_end = runs[r+1].core_end;
                runs[r].rcref_count += runs[r+1].rcref_count;
                runs.erase(runs.begin() + r + 1);
                merged = true;
                break;
            }

            // Count gap core classifications
            int gap_ref_only = 0, gap_total = 0;
            for (size_t g = gap_start; g < gap_end; g++) {
                gap_total++;
                if (cores[g].cls == 1) gap_ref_only++;
            }

            // Merge if no ref-only cores in the gap (all are both/neither/rcref)
            if (gap_ref_only == 0) {
                runs[r].core_end = runs[r+1].core_end;
                runs[r].rcref_count += runs[r+1].rcref_count;
                runs.erase(runs.begin() + r + 1);
                merged = true;
                break;
            }
        }
    }

    // Convert runs to calls, anchoring to reference coordinates
    // Find flanking ref-only cores to map tumor positions → ref positions
    vector<SVCall> calls;
    for (const auto& run : runs) {
        int64_t tumor_span = (int64_t)cores[run.core_end].end - (int64_t)cores[run.core_start].start;
        if (tumor_span < MIN_INV_BP) continue;

        // Look left for nearest ref-only core with unique ref anchor
        int64_t ref_start = -1;
        for (int k = (int)run.core_start - 1; k >= 0; k--) {
            if (cores[k].cls == 1 || cores[k].cls == 3) {
                ref_start = ref_index.anchor(cores[k].hash);
                if (ref_start >= 0) break;
            }
        }

        // Look right for nearest ref-only core with unique ref anchor
        int64_t ref_end = -1;
        for (size_t k = run.core_end + 1; k < cores.size(); k++) {
            if (cores[k].cls == 1 || cores[k].cls == 3) {
                ref_end = ref_index.anchor(cores[k].hash);
                if (ref_end >= 0) break;
            }
        }

        if (ref_start < 0 && ref_end < 0) {
            sv_logger.debug("Large INV skipped (no ref anchor): tumor_pos=" + to_string(cores[run.core_start].start)
                + " tumor_span=" + to_string(tumor_span));
            continue;
        }

        SVCall call;
        call.svtype = "INV";
        if (ref_start >= 0 && ref_end >= 0) {
            call.ref_pos = ref_start;
            call.length = ref_end - ref_start;
        } else if (ref_start >= 0) {
            call.ref_pos = ref_start;
            call.length = tumor_span;
        } else {
            call.ref_pos = ref_end - tumor_span;
            call.length = tumor_span;
        }

        if (call.length < MIN_INV_BP) continue;

        sv_logger.debug("Large INV detected via core-index: ref_pos=" + to_string(call.ref_pos)
            + " len=" + to_string(call.length) + " rcref_cores=" + to_string(run.rcref_count)
            + " (tumor_pos=" + to_string(cores[run.core_start].start) + ")");
        calls.push_back(call);
    }

    sv_logger.debug("detectLargeINVs: " + to_string(cores.size()) + " tumor cores, "
        + to_string(runs.size()) + " runs after merge, "
        + to_string(calls.size()) + " large INVs detected");

    if (DiagLog::instance().enabled()) {
        // Core classification counts
        int cls0 = 0, cls1 = 0, cls2 = 0, cls3 = 0;
        for (const auto& c : cores) {
            if (c.cls == 0) cls0++;
            else if (c.cls == 1) cls1++;
            else if (c.cls == 2) cls2++;
            else cls3++;
        }
        DiagLog::instance().addStat("large_inv_tumor_cores", to_string(cores.size()));
        DiagLog::instance().addStat("large_inv_cls0_novel", to_string(cls0));
        DiagLog::instance().addStat("large_inv_cls1_ref_only", to_string(cls1));
        DiagLog::instance().addStat("large_inv_cls2_rcref_only", to_string(cls2));
        DiagLog::instance().addStat("large_inv_cls3_both", to_string(cls3));
        DiagLog::instance().addStat("large_inv_runs", to_string(runs.size()));
        DiagLog::instance().addStat("large_inv_calls", to_string(calls.size()));

        for (const auto& c : calls) {
            DiagLog::instance().addCall({c.svtype, c.ref_pos, c.length, "large_inv",
                "rcref_core_pipeline"});
        }
    }

    return calls;
}

// ---------------------------------------------------------------------------
// Large DEL detection via ref-core absence in tumor
// ---------------------------------------------------------------------------
vector<SVCall> detectLargeDELs(const RefIndex& ref_index,
    const unordered_set<uint64_t>& tumor_hashes,
    const vector<SVCall>& large_inv_regions) {

    constexpr int MAX_GAP = 3;
    constexpr int MIN_ABSENT_CORES = 5;
    constexpr int64_t MIN_DEL_BP = 500;
    constexpr int64_t INV_OVERLAP_PAD = 500;

    // Classify each reference core:
    // 0 = ambiguous (non-unique hash), 1 = present in tumor, 2 = absent
    const auto& ref_cores = ref_index.cores;
    vector<uint8_t> cls(ref_cores.size());
    int n_present = 0, n_absent = 0, n_ambig = 0;

    for (size_t i = 0; i < ref_cores.size(); i++) {
        auto it = ref_index.index.find(ref_cores[i].hash);
        if (it == ref_index.index.end() || it->second.size() != 1) {
            cls[i] = 0; n_ambig++;
        } else if (tumor_hashes.count(ref_cores[i].hash)) {
            cls[i] = 1; n_present++;
        } else {
            cls[i] = 2; n_absent++;
        }
    }

    // Find gap-tolerant runs of absent cores
    // Key: present (cls=1) cores are HARD boundaries — always split
    //      ambiguous (cls=0) are soft gaps — allow up to MAX_GAP consecutive
    struct DelRun { size_t core_start, core_end; int absent_count; };
    vector<DelRun> runs;

    size_t i = 0;
    while (i < ref_cores.size()) {
        if (cls[i] != 2) { i++; continue; }

        size_t run_start = i;
        int absent_count = 0;
        int soft_gap = 0;
        size_t run_end = i;

        while (i < ref_cores.size()) {
            if (cls[i] == 2) {
                absent_count++;
                soft_gap = 0;
                run_end = i;
                i++;
            } else if (cls[i] == 1) {
                break; // hard boundary
            } else {
                soft_gap++;
                if (soft_gap > MAX_GAP) break;
                i++;
            }
        }

        if (absent_count >= MIN_ABSENT_CORES) {
            runs.push_back({run_start, run_end, absent_count});
        }
    }

    // Filter out runs overlapping large INV regions
    vector<DelRun> filtered_runs;
    for (const auto& run : runs) {
        int64_t run_start_bp = (int64_t)ref_cores[run.core_start].start;
        int64_t run_end_bp = (int64_t)ref_cores[run.core_end].end;
        bool overlaps_inv = false;
        for (const auto& inv : large_inv_regions) {
            if (run_start_bp < inv.ref_pos + inv.length + INV_OVERLAP_PAD &&
                run_end_bp > inv.ref_pos - INV_OVERLAP_PAD) {
                overlaps_inv = true;
                break;
            }
        }
        if (!overlaps_inv) filtered_runs.push_back(run);
    }

    // Convert to SVCalls — tighten boundaries using flanking present cores
    // to avoid overestimating deletion size (absent core run extends past
    // actual breakpoints due to boundary-spanning cores whose hash changed)
    vector<SVCall> calls;
    for (const auto& run : filtered_runs) {
        // Default to outer boundaries of absent core run
        int64_t ref_start = (int64_t)ref_cores[run.core_start].start;
        int64_t ref_end = (int64_t)ref_cores[run.core_end].end;

        // Tighten left: use end of last present core before the run
        for (int k = (int)run.core_start - 1; k >= 0; k--) {
            if (cls[k] == 1) {
                ref_start = (int64_t)ref_cores[k].end;
                break;
            }
        }
        // Tighten right: use start of first present core after the run
        for (size_t k = run.core_end + 1; k < ref_cores.size(); k++) {
            if (cls[k] == 1) {
                ref_end = (int64_t)ref_cores[k].start;
                break;
            }
        }

        int64_t span = ref_end - ref_start;
        if (span < MIN_DEL_BP) continue;

        SVCall call;
        call.svtype = "DEL";
        call.ref_pos = ref_start;
        call.length = span;

        sv_logger.debug("Large DEL detected via core-index: ref_pos=" + to_string(call.ref_pos)
            + " len=" + to_string(call.length) + " absent_cores=" + to_string(run.absent_count));
        calls.push_back(call);
    }

    sv_logger.debug("detectLargeDELs: " + to_string(ref_cores.size()) + " ref cores, "
        + to_string(n_present) + " present, " + to_string(n_absent) + " absent, "
        + to_string(n_ambig) + " ambiguous, "
        + to_string(runs.size()) + " raw runs, "
        + to_string(filtered_runs.size()) + " after INV filter, "
        + to_string(calls.size()) + " large DELs detected");

    if (DiagLog::instance().enabled()) {
        DiagLog::instance().addStat("large_del_ref_cores", to_string(ref_cores.size()));
        DiagLog::instance().addStat("large_del_present", to_string(n_present));
        DiagLog::instance().addStat("large_del_absent", to_string(n_absent));
        DiagLog::instance().addStat("large_del_ambiguous", to_string(n_ambig));
        DiagLog::instance().addStat("large_del_raw_runs", to_string(runs.size()));
        DiagLog::instance().addStat("large_del_after_inv_filter", to_string(filtered_runs.size()));
        DiagLog::instance().addStat("large_del_calls", to_string(calls.size()));
        for (const auto& c : calls) {
            DiagLog::instance().addCall({c.svtype, c.ref_pos, c.length, "large_del",
                "absent_core_pipeline"});
        }
    }

    return calls;
}

// ---------------------------------------------------------------------------
// Large INS detection via novel tumor core runs
// ---------------------------------------------------------------------------
vector<SVCall> detectLargeINS(const string& tumor_seq, int depth,
    const unordered_set<uint64_t>& ref_hashes,
    const unordered_set<uint64_t>& rcref_hashes,
    const RefIndex& ref_index) {

    constexpr int MAX_GAP = 3;
    constexpr int MIN_NOVEL_CORES = 5;
    constexpr int64_t MIN_INS_BP = 500;

    // Extract tumor LCP cores
    struct lps lcp;
    init_lps(&lcp, tumor_seq.c_str(), tumor_seq.size());
    lps_deepen(&lcp, depth);

    // Classify each tumor core:
    // 0 = novel, 1 = ref-only, 2 = rcref-only, 3 = both
    struct ClassifiedCore { uint64_t start, end; uint8_t cls; node_hash hash; };
    vector<ClassifiedCore> cores(lcp.size);
    for (int ci = 0; ci < lcp.size; ++ci) {
        string region(tumor_seq.c_str() + lcp.cores[ci].start, lcp.cores[ci].end - lcp.cores[ci].start);
        uint64_t h = hash_string(region, 0);
        bool in_ref = ref_hashes.count(h) > 0;
        bool in_rcref = rcref_hashes.count(h) > 0;
        uint8_t c = 0;
        if (in_ref && in_rcref) c = 3;
        else if (in_ref) c = 1;
        else if (in_rcref) c = 2;
        cores[ci] = {(uint64_t)lcp.cores[ci].start, (uint64_t)lcp.cores[ci].end, c, h};
    }
    free_lps(&lcp);

    // Find gap-tolerant runs of novel cores (cls=0)
    // ref-matching (cls=1) cores are hard boundaries
    struct InsRun { size_t core_start, core_end; int novel_count; };
    vector<InsRun> runs;

    size_t i = 0;
    while (i < cores.size()) {
        if (cores[i].cls != 0) { i++; continue; }

        size_t run_start = i;
        int novel_count = 0;
        int gap = 0;
        size_t run_end = i;

        while (i < cores.size()) {
            if (cores[i].cls == 0) {
                novel_count++;
                gap = 0;
                run_end = i;
                i++;
            } else if (cores[i].cls == 1) {
                break; // ref-matching = hard boundary
            } else {
                gap++;
                if (gap > MAX_GAP) break;
                i++;
            }
        }

        if (novel_count >= MIN_NOVEL_CORES) {
            runs.push_back({run_start, run_end, novel_count});
        }
    }

    // Convert runs to calls — tighten boundaries using flanking ref-matching
    // cores to avoid overestimating insertion size
    vector<SVCall> calls;
    for (const auto& run : runs) {
        // Default to outer boundaries of novel core run
        int64_t tumor_start = (int64_t)cores[run.core_start].start;
        int64_t tumor_end = (int64_t)cores[run.core_end].end;

        // Tighten left: use end of last ref-matching core before the run
        int64_t ref_pos = -1;
        for (int k = (int)run.core_start - 1; k >= 0; k--) {
            if (cores[k].cls == 1 || cores[k].cls == 3) {
                tumor_start = (int64_t)cores[k].end;
                if (ref_pos < 0) ref_pos = ref_index.anchor(cores[k].hash);
                break;
            }
        }
        // Tighten right: use start of first ref-matching core after the run
        for (size_t k = run.core_end + 1; k < cores.size(); k++) {
            if (cores[k].cls == 1 || cores[k].cls == 3) {
                tumor_end = (int64_t)cores[k].start;
                if (ref_pos < 0) ref_pos = ref_index.anchor(cores[k].hash);
                break;
            }
        }

        int64_t tumor_span = tumor_end - tumor_start;
        if (tumor_span < MIN_INS_BP) continue;

        if (ref_pos < 0) {
            sv_logger.debug("Large INS skipped (no ref anchor): tumor_pos=" + to_string(cores[run.core_start].start)
                + " tumor_span=" + to_string(tumor_span));
            continue;
        }

        SVCall call;
        call.svtype = "INS";
        call.ref_pos = ref_pos;
        call.length = tumor_span;

        sv_logger.debug("Large INS detected via core-index: ref_pos=" + to_string(call.ref_pos)
            + " len=" + to_string(call.length) + " novel_cores=" + to_string(run.novel_count)
            + " (tumor_pos=" + to_string(cores[run.core_start].start) + ")");
        calls.push_back(call);
    }

    sv_logger.debug("detectLargeINS: " + to_string(cores.size()) + " tumor cores, "
        + to_string(runs.size()) + " novel runs, "
        + to_string(calls.size()) + " large INSs detected");

    if (DiagLog::instance().enabled()) {
        DiagLog::instance().addStat("large_ins_novel_runs", to_string(runs.size()));
        DiagLog::instance().addStat("large_ins_calls", to_string(calls.size()));
        for (const auto& c : calls) {
            DiagLog::instance().addCall({c.svtype, c.ref_pos, c.length, "large_ins",
                "novel_core_pipeline"});
        }
    }

    return calls;
}
