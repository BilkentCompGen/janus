#include <algorithm>
#include <string>
#include <unordered_set>
#include <vector>

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include "graph.h"
#include "ref.h"

std::vector<DeletionRegion> parseVCFDeletions(const std::string& vcf_path) {
    std::vector<DeletionRegion> deletions;

    htsFile* fp = bcf_open(vcf_path.c_str(), "r");
    if (!fp) return deletions;

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        bcf_close(fp);
        return deletions;
    }

    bcf1_t* rec = bcf_init();

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_INFO);

        char* svtype = nullptr;
        int n_svtype = 0;
        if (bcf_get_info_string(hdr, rec, "SVTYPE", &svtype, &n_svtype) < 0) continue;

        if (std::string(svtype) != "DEL") {
            free(svtype);
            continue;
        }
        free(svtype);

        int64_t* end_pos_ptr = nullptr;
        int n_end = 0;
        if (bcf_get_info_int64(hdr, rec, "END", &end_pos_ptr, &n_end) < 0 || n_end == 0) continue;
        
        int64_t end_pos = end_pos_ptr[0];
        free(end_pos_ptr);

        uint64_t start_pos = static_cast<uint64_t>(rec->pos) + 1;
        if (end_pos <= static_cast<int64_t>(start_pos)) continue;

        std::string chrom = bcf_hdr_id2name(hdr, rec->rid);
        if (chrom.rfind("chr", 0) == 0) chrom = chrom.substr(3);

        deletions.push_back({chrom, start_pos, static_cast<uint64_t>(end_pos)});
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return deletions;
}

// TODO: Implement findDeletionNodeHashes when node positions are tracked
// This function requires node position tracking which was removed.
// Will be re-implemented when VCF functionality is fully integrated.
std::unordered_set<node_hash> findDeletionNodeHashes(const std::vector<DeletionRegion>& /* deletions */, const LCPGraph& /* graph */) {
    std::unordered_set<node_hash> deletion_hashes;
    
    // Implementation pending - requires node position tracking
    // for (const auto& [hash, node] : graph.nodes) {
    //     auto positions = graph.getPositions(hash, 0);
    //     ...
    // }
    
    return deletion_hashes;
}

