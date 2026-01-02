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
    
    r.truth_count = t_list.size();
    r.pred_count = p_list.size();

    // 1. Calculate all possible pairwise scores
    struct Candidate {
        size_t t_idx, p_idx;
        double score;
    };
    std::vector<Candidate> candidates;

    for (size_t i = 0; i < t_list.size(); ++i) {
        const auto& t = t_list[i];
        for (size_t j = 0; j < p_list.size(); ++j) {
            const auto& p = p_list[j];
            if (t.chrom != p.chrom) continue;

            double score = 0;
            if (type == "SNP") {
                if (p.pos == t.pos) score = 1.0;
            } else if (type == "INS") {
                // Strict Proportional Tolerance (GIAB/Wittyer standard)
                // 25% of size, min 10bp, capped at 500bp
                int64_t min_len = std::min(t.length(), p.length());
                int64_t prop_tolerance = std::clamp((int64_t)(min_len * 0.25), 10L, 500L);
                
                int64_t dist = std::abs((int64_t)p.pos - (int64_t)t.pos);
                if (dist <= prop_tolerance) {
                    // Size similarity score only (distance check is binary filter)
                    double size_ratio = (double)std::min(t.length(), p.length()) / std::max(t.length(), p.length());
                    score = size_ratio; 
                }
            } else { // DEL
                int64_t ov_start = std::max(p.pos, t.pos);
                int64_t ov_end = std::min(p.end, t.end);
                if (ov_start < ov_end) {
                    double overlap_len = (double)(ov_end - ov_start);
                    // Reciprocal overlap score
                    double score_t = overlap_len / t.length();
                    double score_p = overlap_len / p.length();
                    score = std::min(score_t, score_p);
                }
            }

            if (score > 0) candidates.push_back({i, j, score});
        }
    }

    // 2. Sort candidates by Score (Desc) -> Truth Length (Desc)
    std::sort(candidates.begin(), candidates.end(), [&](const Candidate& a, const Candidate& b) {
        if (std::abs(a.score - b.score) > 1e-6) return a.score > b.score;
        return t_list[a.t_idx].length() > t_list[b.t_idx].length();
    });

    // 3. Greedy Matching "Best-First"
    std::vector<bool> t_matched(t_list.size(), false);
    std::vector<bool> p_matched(p_list.size(), false);

    for (const auto& c : candidates) {
        if (c.score < min_overlap) break; // Optimization: sorted list
        if (t_matched[c.t_idx] || p_matched[c.p_idx]) continue;

        t_matched[c.t_idx] = true;
        p_matched[c.p_idx] = true;
        r.tp_count++;
        r.true_positives.push_back({t_list[c.t_idx], p_list[c.p_idx], c.score});
    }

    // 4. Merge Detection: One Prediction covering Multiple GTs (DEL only)
    if (type == "DEL") {
        for (size_t j = 0; j < p_list.size(); ++j) {
            if (p_matched[j]) continue;
            const auto& p = p_list[j];
            
            // Find all unmatched GTs that this prediction overlaps
            std::vector<std::pair<size_t, Variant>> covered;
            for (size_t i = 0; i < t_list.size(); ++i) {
                if (t_matched[i]) continue;
                const auto& t = t_list[i];
                if (p.chrom != t.chrom) continue;
                
                // Check overlap with 100bp tolerance
                int64_t ov = std::max(0L, (int64_t)std::min(p.end + 100, t.end) - (int64_t)std::max(p.pos, t.pos));
                if (ov > 0) covered.push_back({i, t});
            }
            
            // Need at least 2 GTs for a merge
            if (covered.size() >= 2) {
                std::sort(covered.begin(), covered.end(), [](const auto& a, const auto& b) { return a.second.pos < b.second.pos; });
                
                // Check clustering: gaps <= 500bp between consecutive GTs
                bool is_clustered = true;
                for (size_t k = 1; k < covered.size(); ++k) {
                    if ((int64_t)covered[k].second.pos - (int64_t)covered[k-1].second.end > 500) {
                        is_clustered = false;
                        break;
                    }
                }
                
                if (is_clustered) {
                    // Calculate fill ratio
                    int64_t total_gt_len = 0;
                    std::vector<Variant> covered_truths;
                    for (const auto& [idx, v] : covered) {
                        total_gt_len += v.length();
                        covered_truths.push_back(v);
                    }
                    double fill_ratio = (double)total_gt_len / p.length();
                    
                    // Require >= 75% fill
                    if (fill_ratio >= 0.75) {
                        p_matched[j] = true;
                        for (const auto& [idx, v] : covered) t_matched[idx] = true;
                        
                        r.merge_matches.push_back({p, covered_truths, fill_ratio});
                        r.tp_count += covered.size(); // Each GT counts as a TP
                    }
                }
            }
        }
    }

    // 5. Categorize Remaining Variants
    // False Positives: Check if they overlap ANY truth (Support/Fragmentation) or are pure Artifacts
    for (size_t j = 0; j < p_list.size(); ++j) {
        if (p_matched[j]) continue;
        const auto& p = p_list[j];
        
        // Find best overlapping truth (even if already matched)
        int best_t_idx = -1;
        double best_ov = 0.0;

        for (size_t i = 0; i < t_list.size(); ++i) {
            const auto& t = t_list[i];
            if (p.chrom != t.chrom) continue;
            
            double score = 0;
            if (type == "SNP") { if (p.pos == t.pos) score = 1.0; }
            else if (type == "INS") {
                 int64_t dist = std::abs((int64_t)p.pos - (int64_t)t.pos);
                 if (dist <= 1000) score = 1.0; // Loose association
            } else {
                int64_t ov = std::max(0L, (int64_t)std::min(p.end, t.end) - (int64_t)std::max(p.pos, t.pos));
                if (ov > 0) score = (double)ov / p.length(); // How much of P is explained by T?
            }

            if (score > best_ov) { best_ov = score; best_t_idx = i; }
        }

        if (best_t_idx != -1 && best_ov > 0.01) {
            r.fragmented_matches.push_back({p, t_list[best_t_idx], best_ov});
        } else {
            r.artifacts.push_back(p);
        }
        r.fp_count++;
    }

    // False Negatives: Check if they are partially covered (PartialMiss) or completely Missed
    for (size_t i = 0; i < t_list.size(); ++i) {
        if (t_matched[i]) continue;
        const auto& t = t_list[i];

        std::vector<Variant> supporters;
        double covered_len = 0;
        
        for (size_t j = 0; j < p_list.size(); ++j) {
            const auto& p = p_list[j];
            if (p.chrom != t.chrom) continue;
            
            if (type == "DEL") {
                int64_t ov = std::max(0L, (int64_t)std::min(p.end, t.end) - (int64_t)std::max(p.pos, t.pos));
                if (ov > 0) {
                    supporters.push_back(p);
                    covered_len += ov; 
                }
            } else if (type == "INS") {
                 int64_t dist = std::abs((int64_t)p.pos - (int64_t)t.pos);
                 if (dist <= 1000) supporters.push_back(p);
            } else { // SNP
                if (p.pos == t.pos) supporters.push_back(p);
            }
        }

        if (!supporters.empty()) {
            double cov_ratio = (type == "DEL") ? std::min(1.0, covered_len / t.length()) : 1.0;
            r.partial_matches.push_back({t, supporters, cov_ratio});
        } else {
            r.missed.push_back(t);
        }
        r.fn_count++;
    }

    if (r.tp_count + r.fp_count > 0) r.precision = (double)r.tp_count / (r.tp_count + r.fp_count);
    if (r.tp_count + r.fn_count > 0) r.recall = (double)r.tp_count / (r.tp_count + r.fn_count);
    if (r.precision + r.recall > 0) r.f1 = 2 * r.precision * r.recall / (r.precision + r.recall);
    
    return r;
}

void printAnalysisResult(AnalysisResult& r, const std::string& type, const std::string& report_file) {
    auto log = [&](const std::string& msg) {
        if (!report_file.empty()) { std::ofstream f(report_file, std::ios::app); if (f) f << msg << std::endl; }
        else debug_logger.debug(msg);
    };

    struct ColWidths { int t=7, p=7, s=7, l=7; } w;
    char buf[1024];

    // SORTING HELPERS
    auto sort_matches = [](auto& list) {
        std::sort(list.begin(), list.end(), [](const auto& a, const auto& b) {
            if (std::abs(a.score - b.score) > 1e-4) return a.score > b.score;
            return a.truth.length() > b.truth.length();
        });
    };
    auto sort_artifacts = [](auto& list) {
        std::sort(list.begin(), list.end(), [](const auto& a, const auto& b) { return a.length() > b.length(); });
    };
    
    // 1. TRUE POSITIVES (Matches)
    log("\n### " + type + " - TRUE POSITIVES (" + std::to_string(r.true_positives.size()) + ")");
    if (!r.true_positives.empty()) {
        sort_matches(r.true_positives);
        log("Score | GT Position        | GT Len | Pred Position      | Pred Len | Deviation");
        log("------|--------------------|--------|--------------------|----------|----------");
        for (const auto& m : r.true_positives) {
            std::snprintf(buf, sizeof(buf), "%.2f  | %-18s | %-6ld | %-18s | %-8ld | dPOS=%ld dSIZE=%ld",
                m.score,
                (m.truth.chrom + ":" + std::to_string(m.truth.pos)).c_str(), m.truth.length(),
                (m.pred.chrom + ":" + std::to_string(m.pred.pos)).c_str(), m.pred.length(),
                (int64_t)m.truth.pos - (int64_t)m.pred.pos, m.truth.length() - m.pred.length());
            log(buf);
        }
    }

    // 1b. MERGE MATCHES (One Prediction -> Multiple GTs)
    if (!r.merge_matches.empty()) {
        log("\n### " + type + " - MERGE MATCHES (" + std::to_string(r.merge_matches.size()) + " predictions covering multiple GTs)");
        log("These are single predictions that correctly capture multiple clustered ground truth events.");
        for (const auto& mm : r.merge_matches) {
            std::snprintf(buf, sizeof(buf), "\n[PRED] %s:%lu (Len: %ld) - Fill Ratio: %.1f%% - Covers %zu GTs:", 
                mm.pred.chrom.c_str(), mm.pred.pos, mm.pred.length(), mm.fill_ratio * 100.0, mm.covered_truths.size());
            log(buf);
            
            int64_t total_gt_len = 0;
            for (const auto& t : mm.covered_truths) {
                total_gt_len += t.length();
                std::snprintf(buf, sizeof(buf), "     -> GT %s:%lu (Len: %ld)", t.chrom.c_str(), t.pos, t.length());
                log(buf);
            }
            std::snprintf(buf, sizeof(buf), "     Total GT Length: %ld / Pred Length: %ld", total_gt_len, mm.pred.length());
            log(buf);
        }
    }

    // 2. PARTIAL MATCHES (Size Mismatches & Fragmentation)
    size_t frag_count = r.partial_matches.size();
    if (frag_count > 0) {
        std::vector<AnalysisResult::PartialMiss> single_event_mismatches;
        std::vector<AnalysisResult::PartialMiss> fragmented_events;

        for (const auto& pm : r.partial_matches) {
            if (pm.supporting_preds.size() == 1) single_event_mismatches.push_back(pm);
            else fragmented_events.push_back(pm);
        }

        auto print_partial_block = [&](const std::string& title, const std::string& desc, std::vector<AnalysisResult::PartialMiss>& list) {
            if (list.empty()) return;
            log("\n### " + type + " - " + title + " (" + std::to_string(list.size()) + ")");
            log(desc);
            
            std::sort(list.begin(), list.end(), [](const auto& a, const auto& b) {
                return a.truth.length() > b.truth.length();
            });

            for (const auto& pm : list) {
                std::snprintf(buf, sizeof(buf), "\n[GT] %s:%lu (Len: %ld) - Coverage: %.1f%%", 
                    pm.truth.chrom.c_str(), pm.truth.pos, pm.truth.length(), pm.coverage * 100.0);
                log(buf);
                
                auto sorted_preds = pm.supporting_preds;
                std::sort(sorted_preds.begin(), sorted_preds.end(), [](const auto& a, const auto& b) { return a.pos < b.pos; });
                
                std::string pred_line = "     Covered by: ";
                for (size_t i = 0; i < sorted_preds.size(); ++i) {
                    if (i > 0) pred_line += ", ";
                    pred_line += std::to_string(sorted_preds[i].pos) + " (Len: " + std::to_string(sorted_preds[i].length()) + ")";
                    if (i >= 4) { pred_line += " ... (" + std::to_string(sorted_preds.size() - 5) + " more)"; break; }
                }
                log(pred_line);
            }
        };

        print_partial_block("SIZE / POSITION MISMATCHES", 
            "Single predictions that cover the Truth event but failed the strict overlap threshold (Reciprocal Overlap < 0.5).", 
            single_event_mismatches);
            
        print_partial_block("FRAGMENTED MATCHES", 
            "Truth events covered by MULTIPLE split predictions (Shattered calls).", 
            fragmented_events);
    }

    // 3. ARTIFACTS (Pure False Positives)
    if (!r.artifacts.empty()) {
        log("\n### " + type + " - UNMATCHED ARTIFACTS / NOISE (" + std::to_string(r.artifacts.size()) + ")");
        log("Predictions with NO overlap to any Ground Truth.");
        sort_artifacts(r.artifacts);
        
        // Compact printing for short artifacts
        std::string small_artifacts = "";
        int small_count = 0;
        
        log("Listed by Size (Descending):");
        for (const auto& a : r.artifacts) {
            if (a.length() < 20) {
                if (small_count++ > 0) small_artifacts += ", ";
                small_artifacts += std::to_string(a.pos) + "(L" + std::to_string(a.length()) + ")";
                if (small_count % 10 == 0) small_artifacts += "\n     ";
                continue;
            }
            std::snprintf(buf, sizeof(buf), " - %s:%lu (Len: %ld)", a.chrom.c_str(), a.pos, a.length());
            log(buf);
        }
        if (small_count > 0) {
            log(" + " + std::to_string(small_count) + " small artifacts (<20bp):");
            log("     " + small_artifacts);
        }
    }

    // 4. MISSED (Pure False Negatives)
    if (!r.missed.empty()) {
        log("\n### " + type + " - COMPLETELY MISSED (" + std::to_string(r.missed.size()) + ")");
        log("Ground Truth events with NO predictive coverage.");
        sort_artifacts(r.missed);
        for (const auto& m : r.missed) {
            std::snprintf(buf, sizeof(buf), " - %s:%lu (Len: %ld)", m.chrom.c_str(), m.pos, m.length());
            log(buf);
        }
    }
    log("\n" + std::string(60, '=') + "\n");
}

void runFullAnalysis(const std::string& truth_vcf, const std::string& pred_vcf, const std::string& chrom, const std::string& report_file) {
    auto log = [&](const std::string& msg) {
        if (!report_file.empty()) { std::ofstream f(report_file, std::ios::app); if (f) f << msg << std::endl; }
        else debug_logger.debug(msg);
    };

    log("=== VCF VARIANT ANALYSIS ===");
    log("Truth: " + truth_vcf + ", Pred: " + pred_vcf);
    auto truth = parseVCFForAnalysis(truth_vcf, chrom), pred = parseVCFForAnalysis(pred_vcf, chrom);
    log("Loaded " + std::to_string(truth.size()) + " truth, " + std::to_string(pred.size()) + " predicted");

    // FILTER ARTIFACTS: Remove predicted SNPs within 5bp of any predicted Indel
    std::vector<Variant> indels;
    for (const auto& v : pred) if (v.type != "SNP") indels.push_back(v);
    
    if (!indels.empty()) {
        std::vector<Variant> filtered_pred;
        size_t kept_snps = 0, removed_snps = 0;
        
        for (const auto& v : pred) {
            if (v.type == "SNP") {
                bool is_artifact = false;
                for (const auto& i : indels) {
                    if (v.chrom == i.chrom && std::abs((int64_t)v.pos - (int64_t)i.pos) <= 5) {
                        is_artifact = true;
                        break;
                    }
                }
                if (!is_artifact) { filtered_pred.push_back(v); kept_snps++; }
                else removed_snps++;
            } else {
                filtered_pred.push_back(v);
            }
        }
        pred = filtered_pred;
        log("Filtered " + std::to_string(removed_snps) + " SNPs close to Indels (<5bp)");
    }
    log("");

    std::vector<std::string> types = {"DEL", "INS", "SNP"};
    std::map<std::string, AnalysisResult> results;

    // 1. Collect results
    for (const std::string& type : types) {
        results[type] = analyzeVariants(truth, pred, type);
    }

    // 2. Summary Table (All Variants)
    log("### Summary Statistics (All Variants)");
    log("TYPE | Truth | Pred | TP (Exact) | FP (Total) | FN (Total) | Precision | Recall | F1");
    log("-----|-------|------|------------|------------|------------|-----------|--------|----");
    for (const std::string& type : types) {
        const auto& r = results[type];
        char buf[512];
        std::snprintf(buf, sizeof(buf), "%-4s | %5zu | %4zu | %10zu | %10zu | %10zu | %9.4f | %6.4f | %.4f",
            type.c_str(), r.truth_count, r.pred_count, r.tp_count, r.fp_count, r.fn_count, r.precision, r.recall, r.f1);
        log(buf);
    }
    log("");

    // 3. Detailed analysis
    for (const std::string& type : types) {
        if (results[type].truth_count > 0 || results[type].pred_count > 0) {
            printAnalysisResult(results[type], type, report_file);
        }
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

