#include "debug.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <queue>
#include <sstream>
#include <unordered_set>

static Logger debug_logger("debug");

namespace debug {

uint32_t getDebugColor(uint32_t colors) {
    if (colors & (1U << 3)) return COL_BUBBLE_START;
    if (colors & (1U << 4)) return COL_BUBBLE_END;
    if (colors & (1U << 5)) return COL_DEL_MARKED;
    if (colors & (1U << 6)) return COL_INS_MARKED;
    if (colors & (1U << 7)) return COL_SNP_MARKED;
    uint32_t own = colors & 0x07;
    switch (own) {
        case 0x07: return COL_SHARED;
        case 0x06: return COL_NORMAL_TUMOR;
        case 0x05: return COL_REF_TUMOR;
        case 0x04: return COL_TUMOR_ONLY;
        case 0x03: return COL_REF_NORMAL;
        case 0x02: return COL_NORMAL_ONLY;
        case 0x01: return COL_REF_ONLY;
        default: return COL_DEFAULT;
    }
}

std::map<uint32_t, uint32_t> buildColorMap() {
    std::map<uint32_t, uint32_t> m;
    for (uint32_t i = 0; i < 256; ++i) m[i] = getDebugColor(i);
    return m;
}

ColorStats collectColorStats(const LCPGraph& graph) {
    ColorStats s; s.total_nodes = graph.nodes.size();
    for (const auto& [h, n] : graph.nodes) {
        uint32_t c = n.colors;
        if (c & (1U << 3)) { s.bubble_starts++; continue; }
        if (c & (1U << 4)) { s.bubble_ends++; continue; }
        if (c & (1U << 5)) { s.del_marked++; continue; }
        if (c & (1U << 6)) { s.ins_marked++; continue; }
        if (c & (1U << 7)) { s.snp_marked++; continue; }
        uint32_t own = c & 0x07;
        switch (own) {
            case 0x07: s.shared++; break;
            case 0x06: s.normal_tumor++; break;
            case 0x05: s.ref_tumor++; break;
            case 0x04: s.tumor_only++; break;
            case 0x03: s.ref_normal++; break;
            case 0x02: s.normal_only++; break;
            case 0x01: s.ref_only++; break;
            default: s.no_color++; break;
        }
    }
    return s;
}

void logOutput(const std::string& msg, const std::string& report_file) {
    if (!report_file.empty()) {
        std::ofstream f(report_file, std::ios::app);
        if (f) f << msg << std::endl;
    } else {
        debug_logger.debug(msg);
    }
}

void printColorStats(const ColorStats& s, const std::string& report_file) {
    logOutput("=== NODE COLOR STATISTICS ===", report_file);
    logOutput("Total: " + std::to_string(s.total_nodes), report_file);
    logOutput("Bubble Start: " + std::to_string(s.bubble_starts) + ", End: " + std::to_string(s.bubble_ends), report_file);
    logOutput("DEL: " + std::to_string(s.del_marked) + ", INS: " + std::to_string(s.ins_marked) + ", SNP: " + std::to_string(s.snp_marked), report_file);
    logOutput("Tumor-only: " + std::to_string(s.tumor_only) + ", Normal-only: " + std::to_string(s.normal_only) + ", Shared: " + std::to_string(s.shared), report_file);
    logOutput("Ref+Normal: " + std::to_string(s.ref_normal) + ", Ref+Tumor: " + std::to_string(s.ref_tumor) + ", Normal+Tumor: " + std::to_string(s.normal_tumor), report_file);
    logOutput("Ref-only: " + std::to_string(s.ref_only) + ", No color: " + std::to_string(s.no_color), report_file);
}

BubbleStats analyzeBubbles(const std::vector<Bubble>& bubbles, const std::string& normal_seq, const std::string& tumor_seq) {
    BubbleStats s; s.total = bubbles.size();
    for (const auto& b : bubbles) {
        if (b.ref_s == UINT64_MAX || b.norm_s == UINT64_MAX || b.tum_s == UINT64_MAX) { s.invalid_coords++; continue; }
        if (b.norm_s >= b.norm_e || b.norm_e > normal_seq.size() || b.tum_s >= b.tum_e || b.tum_e > tumor_seq.size()) { s.out_of_bounds++; continue; }
        s.valid++;
        size_t n_len = b.norm_e - b.norm_s, t_len = b.tum_e - b.tum_s, max_len = std::max(n_len, t_len);
        if (max_len <= 100) s.size_0_100++;
        else if (max_len <= 1000) s.size_100_1k++;
        else if (max_len <= 10000) s.size_1k_10k++;
        else if (max_len <= 50000) s.size_10k_50k++;
        else s.size_50k_plus++;
        if (max_len > 50000) s.large_bubbles++; else s.aligned++;
    }
    return s;
}

void printBubbleStats(const BubbleStats& s, const std::string& report_file) {
    logOutput("=== BUBBLE STATISTICS ===", report_file);
    logOutput("Total: " + std::to_string(s.total) + ", Valid: " + std::to_string(s.valid) + ", Invalid: " + std::to_string(s.invalid_coords) + ", OOB: " + std::to_string(s.out_of_bounds), report_file);
    logOutput("Size dist: 0-100=" + std::to_string(s.size_0_100) + ", 100-1K=" + std::to_string(s.size_100_1k) + ", 1K-10K=" + std::to_string(s.size_1k_10k) + ", 10K-50K=" + std::to_string(s.size_10k_50k) + ", 50K+=" + std::to_string(s.size_50k_plus), report_file);
    logOutput("Aligned: " + std::to_string(s.aligned) + ", Large (size-based): " + std::to_string(s.large_bubbles), report_file);
}

std::vector<Variant> parseVCFForAnalysis(const std::string& path, const std::string& chrom_filter, const std::vector<std::string>& types) {
    std::vector<Variant> vars;
    std::ifstream file(path); if (!file.is_open()) { debug_logger.debug("Warning: Could not open VCF: " + path); return vars; }
    std::unordered_set<std::string> type_set(types.begin(), types.end());
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string chrom, id, ref, alt, qual, filter, info; uint64_t pos;
        if (!(iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) continue;
        if (!chrom_filter.empty() && chrom != chrom_filter) continue;
        std::string svtype = "SNP"; uint64_t end = pos; int64_t svlen = 0;
        std::istringstream info_stream(info); std::string item;
        while (std::getline(info_stream, item, ';')) {
            size_t eq = item.find('=');
            if (eq != std::string::npos) {
                std::string key = item.substr(0, eq), val = item.substr(eq + 1);
                if (key == "SVTYPE") svtype = val;
                else if (key == "END") end = std::stoull(val);
                else if (key == "SVLEN") svlen = std::stoll(val);
            }
        }
        if (type_set.find(svtype) == type_set.end()) continue;
        vars.push_back({chrom, svtype, pos, end, svlen});
    }
    return vars;
}

AnalysisResult analyzeVariants(const std::vector<Variant>& truth, const std::vector<Variant>& pred, const std::string& type, int64_t pos_tolerance, double min_overlap) {
    AnalysisResult r;
    std::vector<Variant> t_list, p_list;
    for (const auto& v : truth) if (v.type == type) t_list.push_back(v);
    for (const auto& v : pred) if (v.type == type) p_list.push_back(v);
    r.truth_count = t_list.size(); r.pred_count = p_list.size();
    std::vector<bool> truth_matched(t_list.size(), false), pred_matched(p_list.size(), false);
    for (size_t pi = 0; pi < p_list.size(); ++pi) {
        const auto& p = p_list[pi]; int best_ti = -1; double best_score = 0;
        for (size_t ti = 0; ti < t_list.size(); ++ti) {
            if (truth_matched[ti]) continue;
            const auto& t = t_list[ti]; if (p.chrom != t.chrom) continue;
            double score = 0;
            if (type == "SNP") { if (p.pos == t.pos) score = 1.0; }
            else if (type == "INS") { int64_t dist = std::abs((int64_t)p.pos - (int64_t)t.pos); if (dist <= pos_tolerance) score = 1.0 - (double)dist / pos_tolerance; }
            else { int64_t ov_start = std::max(p.pos, t.pos), ov_end = std::min(p.end, t.end); if (ov_start < ov_end) { int64_t ov_len = ov_end - ov_start; score = std::min((double)ov_len / t.length(), (double)ov_len / p.length()); } }
            if (score >= min_overlap && score > best_score) { best_score = score; best_ti = ti; }
        }
        if (best_ti >= 0) { 
            truth_matched[best_ti] = true; 
            pred_matched[pi] = true; 
            r.tp++; 
            r.matches.push_back({t_list[best_ti], p_list[pi], best_score});
        }
    }
    for (size_t ti = 0; ti < t_list.size(); ++ti) if (!truth_matched[ti]) { r.missed.push_back(t_list[ti]); r.fn++; }
    for (size_t pi = 0; pi < p_list.size(); ++pi) if (!pred_matched[pi]) { r.false_pos.push_back(p_list[pi]); r.fp++; }
    if (r.tp + r.fp > 0) r.precision = (double)r.tp / (r.tp + r.fp);
    if (r.tp + r.fn > 0) r.recall = (double)r.tp / (r.tp + r.fn);
    if (r.precision + r.recall > 0) r.f1 = 2 * r.precision * r.recall / (r.precision + r.recall);
    return r;
}

void printAnalysisResult(const AnalysisResult& r, const std::string& type, const std::string& report_file) {
    logOutput("=== " + type + " ANALYSIS ===", report_file);
    logOutput("Truth: " + std::to_string(r.truth_count) + ", Pred: " + std::to_string(r.pred_count), report_file);
    logOutput("TP: " + std::to_string(r.tp) + ", FP: " + std::to_string(r.fp) + ", FN: " + std::to_string(r.fn), report_file);
    char buf[256]; std::snprintf(buf, sizeof(buf), "Precision: %.4f, Recall: %.4f, F1: %.4f", r.precision, r.recall, r.f1);
    logOutput(buf, report_file);
    
    if (!r.matches.empty()) {
        logOutput("Matches details:", report_file);
        for (const auto& m : r.matches) {
            std::snprintf(buf, sizeof(buf), "GT: %s:%lu (Len: %ld) | Pred: %s:%lu (Len: %ld) | Score: %.2f", 
                m.truth.chrom.c_str(), m.truth.pos, m.truth.length(),
                m.pred.chrom.c_str(), m.pred.pos, m.pred.length(), m.score);
            logOutput(buf, report_file);
        }
    }

    if (!r.missed.empty()) {
        logOutput("Missed (False Negatives - GT not called):", report_file);
        for (const auto& v : r.missed) {
             std::snprintf(buf, sizeof(buf), "  GT: %s:%lu (Len: %ld)", v.chrom.c_str(), v.pos, v.length());
             logOutput(buf, report_file);
        }
    }

    if (!r.false_pos.empty()) {
        logOutput("Extra (False Positives - Pred not in GT):", report_file);
        for (const auto& v : r.false_pos) {
             std::snprintf(buf, sizeof(buf), "  Pred: %s:%lu (Len: %ld)", v.chrom.c_str(), v.pos, v.length());
             logOutput(buf, report_file);
        }
    }
}

void runFullAnalysis(const std::string& truth_vcf, const std::string& pred_vcf, const std::string& chrom, const std::string& report_file) {
    logOutput("=== VCF VARIANT ANALYSIS ===", report_file);
    logOutput("Truth: " + truth_vcf + ", Pred: " + pred_vcf, report_file);
    auto truth = parseVCFForAnalysis(truth_vcf, chrom), pred = parseVCFForAnalysis(pred_vcf, chrom);
    logOutput("Loaded " + std::to_string(truth.size()) + " truth, " + std::to_string(pred.size()) + " predicted", report_file);
    for (const std::string& type : {"DEL", "INS", "SNP"}) {
        auto r = analyzeVariants(truth, pred, type);
        if (r.truth_count > 0 || r.pred_count > 0) printAnalysisResult(r, type, report_file);
    }
}

std::vector<TaggedPosition> parseTagFile(const std::string& path) {
    std::vector<TaggedPosition> tags;
    std::ifstream file(path); if (!file.is_open()) { debug_logger.debug("Warning: Could not open tag file: " + path); return tags; }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string chrom, label = "";
        uint64_t pos;
        size_t colon = line.find(':');
        if (colon != std::string::npos) {
            chrom = line.substr(0, colon);
            std::string rest = line.substr(colon + 1);
            size_t space = rest.find_first_of(" \t");
            if (space != std::string::npos) { pos = std::stoull(rest.substr(0, space)); label = rest.substr(space + 1); }
            else pos = std::stoull(rest);
        } else {
            if (!(iss >> chrom >> pos)) continue;
            iss >> label;
        }
        tags.push_back({chrom, pos, label});
        debug_logger.debug("Tag: " + chrom + ":" + std::to_string(pos) + " " + label);
    }
    debug_logger.debug("Loaded " + std::to_string(tags.size()) + " tagged positions");
    return tags;
}

std::unordered_set<node_hash> findNodesNearPositions(const LCPGraph& graph, const std::vector<Bubble>& bubbles, const std::vector<TaggedPosition>& tags, int64_t tolerance) {
    std::unordered_set<node_hash> result;
    for (const auto& tag : tags) {
        for (const auto& b : bubbles) {
            if (b.ref_s == UINT64_MAX) continue;
            int64_t dist = std::abs((int64_t)b.ref_s - (int64_t)tag.pos);
            if (dist <= tolerance) {
                result.insert(b.start_node); result.insert(b.end_node);
                // BFS to find all nodes in bubble
                std::queue<node_hash> q; std::unordered_set<node_hash> visited;
                auto start_it = graph.next_edges.find(b.start_node);
                if (start_it != graph.next_edges.end()) for (node_hash nh : start_it->second) { q.push(nh); visited.insert(nh); }
                while (!q.empty()) {
                    node_hash cur = q.front(); q.pop();
                    if (cur == b.end_node) continue;
                    result.insert(cur);
                    auto next_it = graph.next_edges.find(cur);
                    if (next_it != graph.next_edges.end()) {
                        for (node_hash nh : next_it->second) {
                            if (!visited.count(nh)) { q.push(nh); visited.insert(nh); }
                        }
                    }
                }
                debug_logger.debug("Tagged bubble at " + tag.chrom + ":" + std::to_string(tag.pos) + " -> " + std::to_string(result.size()) + " nodes");
            }
        }
    }
    return result;
}

std::map<node_hash, uint32_t> buildTaggedColorMap(const LCPGraph& graph, const std::unordered_set<node_hash>& tagged_nodes) {
    std::map<node_hash, uint32_t> node_colors;
    constexpr uint32_t BLACK = 0x000000;
    for (const auto& [hash, node] : graph.nodes) {
        if (tagged_nodes.count(hash)) node_colors[hash] = getDebugColor(node.colors);
        else node_colors[hash] = BLACK;
    }
    debug_logger.debug("Tagged color map: " + std::to_string(tagged_nodes.size()) + " colored, " + std::to_string(graph.nodes.size() - tagged_nodes.size()) + " black");
    return node_colors;
}

}

