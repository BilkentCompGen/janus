#include "vcf.h"
#include "graph.h"

#include <htslib/vcf.h>
#include <htslib/hts.h>

#include <unordered_set>
#include <vector>
#include <string>
#include <algorithm>


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

std::unordered_set<node_hash> findDeletionNodeHashes(const std::vector<DeletionRegion>& deletions, const LCPGraph& graph) {
    std::unordered_set<node_hash> deletion_hashes;
    
    std::unordered_map<std::string, std::vector<DeletionRegion>> deletions_by_chrom;
    for (const auto& del : deletions) {
        deletions_by_chrom[del.chrom].push_back(del);
    }

    for (const auto& [hash, node] : graph.nodes) {
        auto positions = graph.getPositions(hash, 0);
        
        for (const auto& pos : positions) {
            std::string core_chrom = pos.chrom;
            if (core_chrom.rfind("chr", 0) == 0) core_chrom = core_chrom.substr(3);
            
            auto chrom_it = deletions_by_chrom.find(core_chrom);
            if (chrom_it == deletions_by_chrom.end()) continue;
            
            for (const auto& del : chrom_it->second) {
                uint64_t overlap_start = std::max(pos.start, del.start);
                uint64_t overlap_end = std::min(pos.end, del.end);
                
                if (overlap_end > overlap_start) {
                    deletion_hashes.insert(hash);
                    goto next_node;
                }
            }
        }
        next_node:;
    }

    return deletion_hashes;
}