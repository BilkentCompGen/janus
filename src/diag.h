#ifndef DIAG_H
#define DIAG_H

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

// Lightweight diagnostic logger for SVP pipeline.
// When enabled (via --debug), captures per-bubble and per-call details
// to a TSV file for post-hoc analysis.

struct DiagBubble {
    uint64_t start_hash, end_hash;
    int64_t anchor_pos;            // -1 if anchoring failed
    std::string anchor_method;     // "direct", "backward", "forward", "failed"
    size_t normal_seq_len, tumor_seq_len;
    std::string skip_reason;       // empty if not skipped
    int calls_emitted;
    std::string call_types;        // e.g. "DEL,SNP,SNP" or "INV"
};

struct DiagCall {
    std::string svtype;
    int64_t ref_pos, length;
    std::string source;            // "bubble", "large_inv", "flipped_dup", "callMulti", "reclassify_inv"
    std::string detail;            // extra info (alignment score, etc.)
};

class DiagLog {
public:
    static DiagLog& instance();

    void enable(const std::string& output_prefix);
    bool enabled() const { return enabled_; }

    void addBubble(const DiagBubble& b);
    void addCall(const DiagCall& c);
    void addStat(const std::string& key, const std::string& value);

    void flush();

private:
    DiagLog() = default;
    bool enabled_ = false;
    std::string prefix_;
    std::vector<DiagBubble> bubbles_;
    std::vector<DiagCall> calls_;
    std::vector<std::pair<std::string, std::string>> stats_;
};

#endif // DIAG_H
