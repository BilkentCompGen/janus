// svp_clean — Remove SNPs and small indels from a FASTA to isolate structural variants.
//
// Two-stage approach for real-data SV calling:
//   1. svp_clean: detect SNPs/small-indels via sequence anchoring, patch them out
//   2. svp:       run the normal SV pipeline on the cleaned tumor
//
// Algorithm:
//   Divide both ref and tumor into fixed-size windows. Hash each window. Windows
//   with identical hashes are exact matches — anchor points. Between anchors,
//   compare the gap sequences. Gaps differing by <= max_indel bp are SNPs or
//   small indels — replace tumor gap with ref gap. Larger gaps are potential SVs
//   and are preserved unchanged.
//
// Why fixed windows instead of LCP cores:
//   LCP core boundaries shift when nearby SNPs change the local parsing context.
//   Fixed windows have stable boundaries and reliably anchor shared regions.

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"
#include "src/utils.h"

using namespace std;

static Logger logger("svp_clean");

// ---------------------------------------------------------------------------
// Fast hash for a window of sequence
// ---------------------------------------------------------------------------
static uint64_t window_hash(const char* seq, size_t len) {
    // FNV-1a 64-bit
    uint64_t h = 14695981039346656037ULL;
    for (size_t i = 0; i < len; ++i) {
        h ^= (uint64_t)(unsigned char)seq[i];
        h *= 1099511628211ULL;
    }
    return h;
}

// ---------------------------------------------------------------------------
// Anchor: a matched window between ref and tumor
// ---------------------------------------------------------------------------
struct Anchor {
    uint32_t ref_pos;
    uint32_t tum_pos;
    uint32_t length;  // = window size
};

// ---------------------------------------------------------------------------
// Find anchors between ref and tumor using fixed-size window hashing
// ---------------------------------------------------------------------------
static vector<Anchor> find_anchors(const string& ref_seq, const string& tum_seq,
                                    int window, int step) {
    // Build sparse hash map from ref: hash every S-th window
    // Only keep hashes that appear exactly once (unique anchors)
    unordered_map<uint64_t, uint32_t> ref_windows;
    unordered_map<uint64_t, uint8_t> ref_count;  // 0=unseen, 1=unique, 2+=dup

    for (uint32_t i = 0; i + window <= ref_seq.size(); i += step) {
        // Skip windows with N's
        bool has_n = false;
        for (int k = 0; k < window; ++k) {
            if (ref_seq[i + k] == 'N' || ref_seq[i + k] == 'n') { has_n = true; break; }
        }
        if (has_n) continue;

        uint64_t h = window_hash(ref_seq.c_str() + i, window);
        auto& cnt = ref_count[h];
        if (cnt == 0) {
            ref_windows[h] = i;
            cnt = 1;
        } else if (cnt == 1) {
            ref_windows.erase(h);  // no longer unique
            cnt = 2;
        }
    }

    logger.debug("  ref unique windows: " + to_string(ref_windows.size()));

    // Scan tumor at EVERY position to find matches (handles SV-caused shifts)
    // Use step=1 for dense scanning — essential when SVs shift positions
    vector<Anchor> raw_anchors;
    for (uint32_t j = 0; j + window <= tum_seq.size(); ++j) {
        // Quick N check on first char
        if (tum_seq[j] == 'N' || tum_seq[j] == 'n') continue;

        uint64_t h = window_hash(tum_seq.c_str() + j, window);
        auto it = ref_windows.find(h);
        if (it != ref_windows.end()) {
            // Verify exact match
            if (ref_seq.compare(it->second, window, tum_seq, j, window) == 0) {
                raw_anchors.push_back({it->second, j, (uint32_t)window});
                j += window - 1;  // skip past this match
            }
        }
    }

    logger.debug("  raw matches: " + to_string(raw_anchors.size()));

    // Sort by tumor position
    sort(raw_anchors.begin(), raw_anchors.end(),
         [](const Anchor& a, const Anchor& b) { return a.tum_pos < b.tum_pos; });

    // Find the dominant offset (ref_pos - tum_pos) using histogram binning.
    // The two assemblies may have different starts (N-padding, telomeres),
    // so the offset varies across the chromosome. Use local offset tracking.
    //
    // Strategy: sweep through sorted matches. Build collinear segments greedily.
    // For each match, check if it extends the current segment (ref increases too).
    // If ref goes backwards, start a new segment. Keep all segments.

    // First: find the longest increasing subsequence of ref_pos (patience sorting)
    // This gives the optimal collinear chain through all structural rearrangements.
    vector<Anchor> chained;
    if (raw_anchors.empty()) return chained;

    // Patience-sort-based LIS on ref_pos values
    // For 500k elements this is O(n log n)
    size_t n = raw_anchors.size();
    vector<uint32_t> tails;     // smallest ref_pos ending an increasing subseq of length i+1
    vector<size_t> tail_idx;    // index in raw_anchors for that element
    vector<int64_t> parent(n, -1);
    vector<size_t> indices;     // for backtracking

    for (size_t i = 0; i < n; ++i) {
        uint32_t rp = raw_anchors[i].ref_pos;
        auto it = lower_bound(tails.begin(), tails.end(), rp);
        size_t pos = it - tails.begin();
        if (it == tails.end()) {
            tails.push_back(rp);
            tail_idx.push_back(i);
        } else {
            *it = rp;
            tail_idx[pos] = i;
        }
        parent[i] = (pos > 0) ? (int64_t)tail_idx[pos - 1] : -1;
    }

    // Backtrack to find the actual LIS
    vector<size_t> lis_indices;
    int64_t idx = tail_idx.back();
    while (idx >= 0) {
        lis_indices.push_back(idx);
        idx = parent[idx];
    }
    reverse(lis_indices.begin(), lis_indices.end());

    for (size_t i : lis_indices)
        chained.push_back(raw_anchors[i]);

    logger.debug("  LIS chain: " + to_string(chained.size()) + " anchors");
    return chained;
}

// ---------------------------------------------------------------------------
// Extend anchors: grow each anchor match in both directions as far as sequences agree
// ---------------------------------------------------------------------------
static void extend_anchors(vector<Anchor>& anchors,
                            const string& ref_seq, const string& tum_seq) {
    for (size_t i = 0; i < anchors.size(); ++i) {
        auto& a = anchors[i];

        // Limits from neighbors
        uint32_t ref_lo = (i > 0) ? anchors[i-1].ref_pos + anchors[i-1].length : 0;
        uint32_t tum_lo = (i > 0) ? anchors[i-1].tum_pos + anchors[i-1].length : 0;
        uint32_t ref_hi = (i + 1 < anchors.size()) ? anchors[i+1].ref_pos : ref_seq.size();
        uint32_t tum_hi = (i + 1 < anchors.size()) ? anchors[i+1].tum_pos : tum_seq.size();

        // Extend left
        while (a.ref_pos > ref_lo && a.tum_pos > tum_lo &&
               ref_seq[a.ref_pos - 1] == tum_seq[a.tum_pos - 1]) {
            a.ref_pos--;
            a.tum_pos--;
            a.length++;
        }

        // Extend right
        uint32_t ref_end = a.ref_pos + a.length;
        uint32_t tum_end = a.tum_pos + a.length;
        while (ref_end < ref_hi && tum_end < tum_hi &&
               ref_seq[ref_end] == tum_seq[tum_end]) {
            ref_end++;
            tum_end++;
            a.length++;
        }
    }
}

// ---------------------------------------------------------------------------
// Build cleaned sequence from anchors + gap patching
// ---------------------------------------------------------------------------
struct CleanStats {
    uint64_t anchors;
    uint64_t anchor_bases;
    uint64_t gaps_total;
    uint64_t gaps_identical;
    uint64_t gaps_patched;
    uint64_t gaps_preserved;
    uint64_t snps_fixed;
    uint64_t indels_fixed;
};

static string build_cleaned(const string& ref_seq, const string& tum_seq,
                              const vector<Anchor>& anchors,
                              int max_indel, CleanStats& stats) {
    string cleaned;
    cleaned.reserve(tum_seq.size());

    stats = {};
    stats.anchors = anchors.size();

    uint32_t tum_pos = 0;

    for (size_t i = 0; i < anchors.size(); ++i) {
        const auto& a = anchors[i];

        // --- Gap before this anchor ---
        uint32_t tum_gap_end = a.tum_pos;
        if (tum_gap_end > tum_pos) {
            uint32_t ref_gap_start = (i > 0) ?
                anchors[i-1].ref_pos + anchors[i-1].length : 0;
            uint32_t ref_gap_end = a.ref_pos;

            string ref_gap = ref_seq.substr(ref_gap_start, ref_gap_end - ref_gap_start);
            string tum_gap = tum_seq.substr(tum_pos, tum_gap_end - tum_pos);

            stats.gaps_total++;

            if (ref_gap == tum_gap) {
                cleaned += tum_gap;
                stats.gaps_identical++;
            } else {
                int size_diff = abs((int)ref_gap.size() - (int)tum_gap.size());
                if (size_diff <= max_indel) {
                    cleaned += ref_gap;
                    stats.gaps_patched++;
                    if (ref_gap.size() == tum_gap.size())
                        stats.snps_fixed++;
                    else
                        stats.indels_fixed++;
                } else {
                    cleaned += tum_gap;
                    stats.gaps_preserved++;
                }
            }
        }

        // --- Anchor region: copy from tumor (it's identical to ref) ---
        cleaned += tum_seq.substr(a.tum_pos, a.length);
        stats.anchor_bases += a.length;
        tum_pos = a.tum_pos + a.length;
    }

    // --- Trailing gap ---
    if (tum_pos < tum_seq.size()) {
        uint32_t ref_trail_start = anchors.empty() ? 0 :
            anchors.back().ref_pos + anchors.back().length;
        string ref_trail = (ref_trail_start < ref_seq.size()) ?
            ref_seq.substr(ref_trail_start) : "";
        string tum_trail = tum_seq.substr(tum_pos);

        stats.gaps_total++;
        if (ref_trail == tum_trail) {
            cleaned += tum_trail;
            stats.gaps_identical++;
        } else {
            int size_diff = abs((int)ref_trail.size() - (int)tum_trail.size());
            if (size_diff <= max_indel) {
                cleaned += ref_trail;
                stats.gaps_patched++;
            } else {
                cleaned += tum_trail;
                stats.gaps_preserved++;
            }
        }
    }

    return cleaned;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    LCP_INIT();

    ArgParser parser("Remove SNPs and small indels from tumor FASTA using window-hash anchoring.");
    parser.add("ref",       "--ref",       "-r", "Reference FASTA",                      "string", true);
    parser.add("tumor",     "--tumor",     "-t", "Tumor FASTA to clean",                 "string", true);
    parser.add("chr",       "--chr",       "-c", "Chromosome to process (0 = all 1-22)", "int",    false);
    parser.add("depth",     "--depth",     "-d", "LCP depth (unused, for compat)",       "int",    false);
    parser.add("window",    "--window",    "-w", "Anchor window size (default 100)",     "int",    false);
    parser.add("max_indel", "--max-indel", "-m", "Max indel size to clean (default 10)", "int",    false);
    parser.add("output",    "--output",    "-o", "Output cleaned FASTA path",            "string", true);

    auto args = parser.parse(argc, argv);
    string ref_file    = any_cast<string>(args["ref"]);
    string tumor_file  = any_cast<string>(args["tumor"]);
    string output_file = any_cast<string>(args["output"]);
    int chr            = args.count("chr") ? any_cast<int>(args["chr"]) : 0;
    int window         = args.count("window") ? any_cast<int>(args["window"]) : 100;
    int max_indel      = args.count("max_indel") ? any_cast<int>(args["max_indel"]) : 10;
    int step           = window;  // non-overlapping windows

    bool all_chr = (chr == 0);

    ofstream out(output_file);
    if (!out.is_open()) {
        cerr << "ERROR: cannot open output file: " << output_file << endl;
        return 1;
    }

    vector<string> chr_list;
    if (all_chr) {
        for (int c = 1; c <= 22; ++c) chr_list.push_back(to_string(c));
    } else {
        chr_list.push_back(to_string(chr));
    }

    uint64_t total_patched = 0, total_preserved = 0, total_snps = 0, total_indels = 0;

    for (const string& chr_str : chr_list) {
        // --- Read reference chromosome ---
        string ref_seq;
        {
            FastaStreamer streamer(ref_file, false);
            for (const FastaPair& pair : streamer) {
                string chrom = pair.first;
                size_t p = chrom.find("_hap");
                if (p != string::npos) chrom = chrom.substr(0, p);
                if (chrom == chr_str) { ref_seq = pair.second; break; }
            }
        }
        if (ref_seq.empty()) {
            logger.debug("Skipping chr" + chr_str + ": not in reference");
            continue;
        }

        // --- Read tumor chromosome ---
        string tum_seq;
        {
            FastaStreamer streamer(tumor_file, false);
            for (const FastaPair& pair : streamer) {
                string chrom = pair.first;
                string bare = chrom;
                size_t hap_pos = bare.find("_hap");
                if (hap_pos != string::npos) bare = bare.substr(0, hap_pos);
                for (const char* suffix : {"_PATERNAL", "_MATERNAL"}) {
                    size_t p = bare.find(suffix);
                    if (p != string::npos) { bare = bare.substr(0, p); break; }
                }
                if (bare.substr(0, 3) == "chr") bare = bare.substr(3);
                if (bare == chr_str) { tum_seq = pair.second; break; }
            }
        }
        if (tum_seq.empty()) {
            logger.debug("Skipping chr" + chr_str + ": not in tumor");
            continue;
        }

        logger.debug("chr" + chr_str + ": ref=" + to_string(ref_seq.size()) +
                     "bp tumor=" + to_string(tum_seq.size()) + "bp");

        // --- Find anchors ---
        auto anchors = find_anchors(ref_seq, tum_seq, window, step);
        logger.debug("chr" + chr_str + ": raw anchors=" + to_string(anchors.size()));

        // --- Extend anchors ---
        extend_anchors(anchors, ref_seq, tum_seq);

        uint64_t anchor_bp = 0;
        for (const auto& a : anchors) anchor_bp += a.length;
        logger.debug("chr" + chr_str + ": after extend, anchor_bp=" + to_string(anchor_bp) +
                     " (" + to_string(anchor_bp * 100 / tum_seq.size()) + "% of tumor)");

        // --- Build cleaned sequence ---
        CleanStats stats;
        string cleaned = build_cleaned(ref_seq, tum_seq, anchors, max_indel, stats);

        logger.debug("chr" + chr_str + ": gaps=" + to_string(stats.gaps_total) +
                     " identical=" + to_string(stats.gaps_identical) +
                     " patched=" + to_string(stats.gaps_patched) +
                     " (snps=" + to_string(stats.snps_fixed) +
                     " indels=" + to_string(stats.indels_fixed) + ")" +
                     " preserved=" + to_string(stats.gaps_preserved));

        total_patched   += stats.gaps_patched;
        total_preserved += stats.gaps_preserved;
        total_snps      += stats.snps_fixed;
        total_indels    += stats.indels_fixed;

        // --- Write cleaned chromosome ---
        out << ">" << chr_str << "\n";
        for (size_t i = 0; i < cleaned.size(); i += 80)
            out << cleaned.substr(i, min((size_t)80, cleaned.size() - i)) << "\n";
    }

    out.close();

    logger.debug("TOTAL: patched=" + to_string(total_patched) +
                 " (snps=" + to_string(total_snps) +
                 " indels=" + to_string(total_indels) + ")" +
                 " preserved=" + to_string(total_preserved));
    logger.debug("Cleaned FASTA: " + output_file);

    return 0;
}
