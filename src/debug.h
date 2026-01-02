#ifndef DEBUG_H
#define DEBUG_H

#include <cstdint>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#include "graph.h"

namespace debug {

// Color constants (RGB hex)
constexpr uint32_t COL_BUBBLE_START = 0x32CD32, COL_BUBBLE_END = 0xFF4500;
constexpr uint32_t COL_DEL_MARKED = 0xFFD700, COL_INS_MARKED = 0x00FFFF, COL_SNP_MARKED = 0xFF69B4;
constexpr uint32_t COL_TUMOR_ONLY = 0xDC143C, COL_NORMAL_ONLY = 0x4169E1, COL_SHARED = 0x90EE90;
constexpr uint32_t COL_REF_NORMAL = 0x87CEEB, COL_REF_TUMOR = 0xFF6347, COL_NORMAL_TUMOR = 0x9370DB;
constexpr uint32_t COL_REF_ONLY = 0xC0C0C0, COL_DEFAULT = 0x404040;

uint32_t getDebugColor(uint32_t colors);
std::map<uint32_t, uint32_t> buildColorMap();

struct ColorStats {
    size_t total_nodes = 0, bubble_starts = 0, bubble_ends = 0;
    size_t del_marked = 0, ins_marked = 0, snp_marked = 0;
    size_t tumor_only = 0, normal_only = 0, shared = 0;
    size_t ref_normal = 0, ref_tumor = 0, normal_tumor = 0, ref_only = 0, no_color = 0;
};

ColorStats collectColorStats(const LCPGraph& graph);
void printColorStats(const ColorStats& stats, const std::string& report_file = "");

struct BubbleStats {
    size_t total = 0, valid = 0, invalid_coords = 0, out_of_bounds = 0, large_bubbles = 0, aligned = 0;
    size_t size_0_100 = 0, size_100_1k = 0, size_1k_10k = 0, size_10k_50k = 0, size_50k_plus = 0;
};

BubbleStats analyzeBubbles(const std::vector<Bubble>& bubbles, const std::string& normal_seq, const std::string& tumor_seq);
void printBubbleStats(const BubbleStats& stats, const std::string& report_file = "");

struct Variant {
    std::string chrom, type;
    uint64_t pos, end;
    int64_t svlen;
    int64_t length() const { return type == "SNP" ? 1 : (svlen != 0 ? std::abs(svlen) : (int64_t)(end - pos)); }
};

struct AnalysisResult {
    struct Match {
        Variant truth, pred;
        double score;
    };
    // False Positive that actually supports a Truth variant (but wasn't the best match)
    struct FragmentedSupport {
        Variant pred;
        Variant target_truth;
        double overlap;
    };
    // False Negative that was partially covered but not well enough to be a Match
    struct PartialMiss {
        Variant truth;
        std::vector<Variant> supporting_preds;
        double coverage;
    };
    // Merge Match: One prediction covers multiple GT events
    struct MergeMatch {
        Variant pred;
        std::vector<Variant> covered_truths;
        double fill_ratio;  // sum(GT lengths) / pred length
    };

    size_t truth_count = 0, pred_count = 0;
    size_t tp_count = 0, fp_count = 0, fn_count = 0;
    double precision = 0.0, recall = 0.0, f1 = 0.0;

    std::vector<Match> true_positives;
    std::vector<MergeMatch> merge_matches;              // One pred -> Many GTs
    std::vector<FragmentedSupport> fragmented_matches;  // FP but useful
    std::vector<PartialMiss> partial_matches;           // FN but touched
    std::vector<Variant> artifacts;                     // Pure FP (no overlap)
    std::vector<Variant> missed;                        // Pure FN (no coverage)
};

std::vector<Variant> parseVCFForAnalysis(const std::string& path, const std::string& chrom_filter = "", const std::vector<std::string>& types = {"DEL", "INS", "SNP"});
AnalysisResult analyzeVariants(const std::vector<Variant>& truth, const std::vector<Variant>& pred, const std::string& type, int64_t pos_tolerance = 50, double min_overlap = 0.5);
void printAnalysisResult(const AnalysisResult& result, const std::string& type, const std::string& report_file = "");
void runFullAnalysis(const std::string& truth_vcf, const std::string& pred_vcf, const std::string& chrom = "", const std::string& report_file = "");

// Tagged variant highlighting - only show specific regions, everything else black
struct TaggedPosition { std::string chrom; uint64_t pos; std::string label; };
std::vector<TaggedPosition> parseTagFile(const std::string& path);
std::unordered_set<node_hash> findNodesNearPositions(const LCPGraph& graph, const std::vector<Bubble>& bubbles, const std::vector<TaggedPosition>& tags, int64_t tolerance = 100);
std::map<node_hash, uint32_t> buildTaggedColorMap(const LCPGraph& graph, const std::unordered_set<node_hash>& tagged_nodes);

}

#endif
