#include "diag.h"
#include <iostream>

DiagLog& DiagLog::instance() {
    static DiagLog log;
    return log;
}

void DiagLog::enable(const std::string& output_prefix) {
    enabled_ = true;
    prefix_ = output_prefix;
}

void DiagLog::addBubble(const DiagBubble& b) {
    if (!enabled_) return;
    bubbles_.push_back(b);
}

void DiagLog::addCall(const DiagCall& c) {
    if (!enabled_) return;
    calls_.push_back(c);
}

void DiagLog::addStat(const std::string& key, const std::string& value) {
    if (!enabled_) return;
    stats_.emplace_back(key, value);
}

void DiagLog::flush() {
    if (!enabled_) return;

    // Write bubble diagnostics
    {
        std::string path = prefix_ + "_bubbles.tsv";
        std::ofstream out(path);
        out << "start_hash\tend_hash\tanchor_pos\tanchor_method\tnormal_len\ttumor_len\tskip_reason\tcalls_emitted\tcall_types\n";
        for (const auto& b : bubbles_) {
            out << b.start_hash << "\t" << b.end_hash << "\t"
                << b.anchor_pos << "\t" << b.anchor_method << "\t"
                << b.normal_seq_len << "\t" << b.tumor_seq_len << "\t"
                << (b.skip_reason.empty() ? "." : b.skip_reason) << "\t"
                << b.calls_emitted << "\t"
                << (b.call_types.empty() ? "." : b.call_types) << "\n";
        }
        std::cerr << "[diag] Wrote " << bubbles_.size() << " bubble records -> " << path << "\n";
    }

    // Write call diagnostics
    {
        std::string path = prefix_ + "_calls.tsv";
        std::ofstream out(path);
        out << "svtype\tref_pos\tlength\tsource\tdetail\n";
        for (const auto& c : calls_) {
            out << c.svtype << "\t" << c.ref_pos << "\t" << c.length << "\t"
                << c.source << "\t"
                << (c.detail.empty() ? "." : c.detail) << "\n";
        }
        std::cerr << "[diag] Wrote " << calls_.size() << " call records -> " << path << "\n";
    }

    // Write summary stats
    {
        std::string path = prefix_ + "_stats.tsv";
        std::ofstream out(path);
        out << "key\tvalue\n";
        for (const auto& [k, v] : stats_) {
            out << k << "\t" << v << "\n";
        }
        std::cerr << "[diag] Wrote " << stats_.size() << " stats -> " << path << "\n";
    }
}
