#include <algorithm>
#include <queue>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include "align.h"
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

std::unordered_set<node_hash> findDeletionNodeHashes(
    const std::vector<DeletionRegion>& deletions, 
    const LCPGraph& graph, 
    uint32_t ref_seq_id
) {
    std::unordered_set<node_hash> deletion_hashes;
    
    for (const auto& [hash, node] : graph.nodes) {
        for (const auto& occ : node.occurrences) {
            if (occ.seq_id == ref_seq_id) {
                for (const auto& del : deletions) {
                    if (occ.start < del.end && occ.end > del.start) {
                        deletion_hashes.insert(hash);
                        break;
                    }
                }
            }
            if (deletion_hashes.count(hash)) break;
        }
    }
    
    return deletion_hashes;
}

std::vector<Region> mergeOverlappingRegions(std::vector<Region>& regions) {
    if (regions.empty()) return {};
    
    // Sort by start position
    std::sort(regions.begin(), regions.end(), [](const Region& a, const Region& b) {
        return a.start < b.start;
    });
    
    std::vector<Region> merged;
    merged.push_back(regions[0]);
    
    for (size_t i = 1; i < regions.size(); ++i) {
        Region& last = merged.back();
        const Region& curr = regions[i];
        
        // Check for overlap (including adjacent)
        if (curr.start <= last.end) {
            // Merge: extend end if necessary
            last.end = std::max(last.end, curr.end);
        } else {
            // No overlap, add new region
            merged.push_back(curr);
        }
    }
    
    return merged;
}

std::pair<uint64_t, uint64_t> findAnchorWindow(
    node_hash target, 
    const LCPGraph& graph,
    uint32_t ref_seq_id
) {
    uint64_t pre_end = 0;
    uint64_t post_start = UINT64_MAX;
    
    // BFS backward to find pre-anchor
    {
        std::queue<node_hash> q;
        std::unordered_set<node_hash> visited;
        q.push(target);
        visited.insert(target);
        
        while (!q.empty()) {
            node_hash curr = q.front();
            q.pop();
            
            auto node_it = graph.nodes.find(curr);
            if (node_it != graph.nodes.end()) {
                // Check if this node has reference occurrence
                for (const auto& occ : node_it->second.occurrences) {
                    if (occ.seq_id == ref_seq_id) {
                        pre_end = std::max(pre_end, occ.end);
                        goto found_pre;
                    }
                }
            }
            
            // Continue BFS backward
            auto prev_it = graph.prev_edges.find(curr);
            if (prev_it != graph.prev_edges.end()) {
                for (node_hash prev : prev_it->second) {
                    if (visited.find(prev) == visited.end()) {
                        visited.insert(prev);
                        q.push(prev);
                    }
                }
            }
        }
        found_pre:;
    }
    
    // BFS forward to find post-anchor
    {
        std::queue<node_hash> q;
        std::unordered_set<node_hash> visited;
        q.push(target);
        visited.insert(target);
        
        while (!q.empty()) {
            node_hash curr = q.front();
            q.pop();
            
            auto node_it = graph.nodes.find(curr);
            if (node_it != graph.nodes.end()) {
                for (const auto& occ : node_it->second.occurrences) {
                    if (occ.seq_id == ref_seq_id) {
                        post_start = std::min(post_start, occ.start);
                        goto found_post;
                    }
                }
            }
            
            auto next_it = graph.next_edges.find(curr);
            if (next_it != graph.next_edges.end()) {
                for (node_hash next : next_it->second) {
                    if (visited.find(next) == visited.end()) {
                        visited.insert(next);
                        q.push(next);
                    }
                }
            }
        }
        found_post:;
    }
    
    return {pre_end, post_start};
}

std::vector<VariantRegion> collectVariantRegions(
    const LCPGraph& graph, 
    uint32_t ref_seq_id,
    const std::string& chrom
) {
    std::vector<Region> del_regions;  // Normal-only (DEL in tumor)
    
    for (const auto& [hash, node] : graph.nodes) {
        uint32_t colors = node.colors;
        bool has_normal = (colors & 0x2);  // Color 1
        bool has_tumor = (colors & 0x4);   // Color 2
        
        if (has_normal && !has_tumor) {
            // Normal-only node -> DEL in tumor
            // Get reference coordinates directly
            for (const auto& occ : node.occurrences) {
                if (occ.seq_id == ref_seq_id) {
                    del_regions.push_back({occ.start, occ.end});
                    break;
                }
            }
        }
        // INS detection disabled for now
    }
    
    // Merge overlapping regions
    std::vector<Region> merged_del = mergeOverlappingRegions(del_regions);
    
    // Convert to VariantRegion
    std::vector<VariantRegion> variants;
    
    for (const auto& r : merged_del) {
        variants.push_back({chrom, r.start, r.end, "DEL"});
    }
    
    // Sort by position
    std::sort(variants.begin(), variants.end(), [](const VariantRegion& a, const VariantRegion& b) {
        return a.ref_start < b.ref_start;
    });
    
    return variants;
}

void exportVCF(
    const std::string& vcf_path,
    const std::vector<VariantRegion>& regions
) {
    htsFile* fp = bcf_open(vcf_path.c_str(), "w");
    if (!fp) return;
    
    bcf_hdr_t* hdr = bcf_hdr_init("w");
    
    // Add VCF version
    bcf_hdr_append(hdr, "##fileformat=VCFv4.2");
    
    // Add contig (assuming chromosome 19 for now)
    bcf_hdr_append(hdr, "##contig=<ID=19,length=58617616>");
    
    // Add INFO fields
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">");
    
    // No sample columns - sites-only VCF
    bcf_hdr_add_sample(hdr, nullptr);
    
    // Write header
    if (bcf_hdr_write(fp, hdr) < 0) {
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return;
    }
    
    // Write each variant
    bcf1_t* rec = bcf_init();
    
    for (const auto& var : regions) {
        bcf_clear(rec);
        
        // Set chromosome
        rec->rid = bcf_hdr_name2id(hdr, var.chrom.c_str());
        if (rec->rid < 0) continue;
        
        // Set position (0-based in htslib)
        rec->pos = static_cast<hts_pos_t>(var.ref_start - 1);
        
        // Set REF and ALT
        std::string alleles = "N,<" + var.svtype + ">";
        bcf_update_alleles_str(hdr, rec, alleles.c_str());
        
        // Set QUAL to missing (.)
        bcf_float_set_missing(rec->qual);
        
        // Set FILTER to PASS
        int32_t pass = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
        bcf_update_filter(hdr, rec, &pass, 1);
        
        // Set INFO fields
        bcf_update_info_string(hdr, rec, "SVTYPE", var.svtype.c_str());
        
        int32_t end_pos = static_cast<int32_t>(var.ref_end);
        bcf_update_info_int32(hdr, rec, "END", &end_pos, 1);
        
        int32_t svlen = static_cast<int32_t>(var.ref_end - var.ref_start);
        if (var.svtype == "DEL") svlen = -svlen;
        bcf_update_info_int32(hdr, rec, "SVLEN", &svlen, 1);
        
        // Write record
        if (bcf_write(fp, hdr, rec) < 0) continue;
    }
    
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}

void exportAlignedVCF(
    const std::string& vcf_path,
    const std::vector<VariantCall>& variants
) {
    htsFile* fp = bcf_open(vcf_path.c_str(), "w");
    if (!fp) return;
    
    bcf_hdr_t* hdr = bcf_hdr_init("w");
    
    // Add VCF version
    bcf_hdr_append(hdr, "##fileformat=VCFv4.2");
    
    // Add contig
    bcf_hdr_append(hdr, "##contig=<ID=19,length=58617616>");
    
    // Add INFO fields
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variant\">");
    
    // No sample columns
    bcf_hdr_add_sample(hdr, nullptr);
    
    if (bcf_hdr_write(fp, hdr) < 0) {
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return;
    }
    
    bcf1_t* rec = bcf_init();
    
    for (const auto& var : variants) {
        bcf_clear(rec);
        
        // Set chromosome
        rec->rid = bcf_hdr_name2id(hdr, var.chrom.c_str());
        if (rec->rid < 0) continue;
        
        // Set position (0-based in htslib)
        rec->pos = static_cast<hts_pos_t>(var.pos - 1);
        
        // Set REF and ALT alleles
        std::string alleles = var.ref_allele + "," + var.alt_allele;
        bcf_update_alleles_str(hdr, rec, alleles.c_str());
        
        // Set QUAL to missing
        bcf_float_set_missing(rec->qual);
        
        // Set FILTER to PASS
        int32_t pass = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
        bcf_update_filter(hdr, rec, &pass, 1);
        
        // Set INFO fields
        bcf_update_info_string(hdr, rec, "SVTYPE", var.svtype.c_str());
        
        int32_t end_pos = static_cast<int32_t>(var.pos + std::abs(var.svlen));
        bcf_update_info_int32(hdr, rec, "END", &end_pos, 1);
        
        int32_t svlen = static_cast<int32_t>(var.svlen);
        bcf_update_info_int32(hdr, rec, "SVLEN", &svlen, 1);
        
        if (bcf_write(fp, hdr, rec) < 0) continue;
    }
    
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
}
