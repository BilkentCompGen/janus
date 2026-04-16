#include <cstdint>
#include <string>
#include <unordered_set>

#include "lib/lcptools/encoding.h"
#include "lib/lcptools/lps.h"

#include "pipeline.h"
#include "utils.h"

using namespace std;

static Logger pipeline_logger("pipeline");

void processGenome(const string& file, int color, uint32_t& sid, string& seq,
                   LCPGraph& graph, const string& chr_str, int depth) {
    FastaStreamer streamer(file, false);
    for (const FastaPair& pair : streamer) {
        string chrom = pair.first;
        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
        if (chrom != chr_str) continue;

        sid = graph.seq_names.size();
        graph.seq_names.push_back(pair.first);
        seq = pair.second;

        struct lps lcp_string;
        init_lps(&lcp_string, seq.c_str(), seq.size());
        lps_deepen(&lcp_string, depth);
        graph.addString(&lcp_string, color, seq.c_str(), sid);
        free_lps(&lcp_string);
    }
}

unordered_set<uint64_t> buildRefHashes(const string& ref_file, const string& chr_str,
                                        int depth, RefIndex& ref_index,
                                        string& ref_sequence) {
    unordered_set<uint64_t> ref_hashes;
    FastaStreamer ref_streamer(ref_file, false);
    for (const FastaPair& pair : ref_streamer) {
        string chrom = pair.first;
        size_t hap_pos = chrom.find("_hap");
        if (hap_pos != string::npos) chrom = chrom.substr(0, hap_pos);
        if (chrom != chr_str) continue;

        ref_sequence = pair.second;
        struct lps lcp_string;
        init_lps(&lcp_string, ref_sequence.c_str(), ref_sequence.size());
        lps_deepen(&lcp_string, depth);
        ref_index.addString(&lcp_string, ref_sequence.c_str());
        for (int i = 0; i < lcp_string.size; ++i) {
            string region(ref_sequence.c_str() + lcp_string.cores[i].start,
                          lcp_string.cores[i].end - lcp_string.cores[i].start);
            ref_hashes.insert(hash_string(region, 0));
        }
        free_lps(&lcp_string);
    }
    return ref_hashes;
}

unordered_set<uint64_t> buildRcrefHashes(const string& ref_sequence, int depth) {
    string rcref = reverse_complement(ref_sequence);
    struct lps rcref_lcp;
    init_lps(&rcref_lcp, rcref.c_str(), rcref.size());
    lps_deepen(&rcref_lcp, depth);

    unordered_set<uint64_t> rcref_hashes;
    for (int i = 0; i < rcref_lcp.size; ++i) {
        string region(rcref.c_str() + rcref_lcp.cores[i].start,
                      rcref_lcp.cores[i].end - rcref_lcp.cores[i].start);
        rcref_hashes.insert(hash_string(region, 0));
    }
    free_lps(&rcref_lcp);

    pipeline_logger.debug("Built rcref hash set: " + to_string(rcref_hashes.size()) + " unique hashes");
    return rcref_hashes;
}

unordered_set<uint64_t> buildTumorHashes(const string& tumor_sequence, int depth) {
    struct lps tumor_lcp;
    init_lps(&tumor_lcp, tumor_sequence.c_str(), tumor_sequence.size());
    lps_deepen(&tumor_lcp, depth);

    unordered_set<uint64_t> tumor_hashes;
    for (int i = 0; i < tumor_lcp.size; ++i) {
        string region(tumor_sequence.c_str() + tumor_lcp.cores[i].start,
                      tumor_lcp.cores[i].end - tumor_lcp.cores[i].start);
        tumor_hashes.insert(hash_string(region, 0));
    }
    free_lps(&tumor_lcp);

    pipeline_logger.debug("Built tumor hash set: " + to_string(tumor_hashes.size()) + " unique hashes");
    return tumor_hashes;
}
