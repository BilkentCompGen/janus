#include "align.h"

#include <parasail.h>
#include <algorithm>
#include <cstdint>

// ============================================================================
// Helper Functions
// ============================================================================

VarType mapType(const std::string& t) {
    if (t == "DEL") return VarType::DEL;
    if (t == "INS") return VarType::INS;
    if (t == "SNP") return VarType::SNP;
    return VarType::COMPLEX;
}

std::string unmapType(VarType t) {
    switch(t) {
        case VarType::DEL: return "DEL";
        case VarType::INS: return "INS";
        case VarType::SNP: return "SNP";
        default: return "COMPLEX";
    }
}

ChatterMetrics detectChatter(const std::vector<MergeVariant>& vars, 
                             int64_t start, int64_t end,
                             double density_threshold) {
    ChatterMetrics m = {0, 0, 0.0, false};
    for (const auto& v : vars) {
        if (v.pos >= start && v.pos < end) m.variant_count++;
    }
    m.span = static_cast<int32_t>(end - start);
    m.density = m.span > 0 ? static_cast<double>(m.variant_count) / m.span : 0;
    m.is_chatter = m.density > density_threshold && m.variant_count > 3;
    return m;
}

static std::vector<std::vector<size_t>> findClusters(const std::vector<MergeVariant>& vars, 
                                                      int32_t proximity_window = 50,
                                                      int32_t max_cluster_span = 500) {
    std::vector<std::vector<size_t>> clusters;
    if (vars.empty()) return clusters;
    
    std::vector<size_t> current_cluster = {0};
    int64_t cluster_start = vars[0].pos;
    
    for (size_t i = 1; i < vars.size(); ++i) {
        int64_t gap = vars[i].pos - vars[current_cluster.back()].end();
        int64_t span = vars[i].end() - cluster_start;
        bool within_proximity = gap <= proximity_window;
        bool within_span = span <= max_cluster_span;
        bool same_chrom = vars[i].chrom == vars[current_cluster[0]].chrom;
        
        if (within_proximity && within_span && same_chrom) {
            current_cluster.push_back(i);
        } else {
            clusters.push_back(std::move(current_cluster));
            current_cluster = {i};
            cluster_start = vars[i].pos;
        }
    }
    clusters.push_back(std::move(current_cluster));
    return clusters;
}

static MergeVariant mergeCluster(const std::vector<MergeVariant>& vars, 
                                  const std::vector<size_t>& indices) {
    if (indices.size() == 1) return vars[indices[0]];
    
    int32_t total_del = 0, total_ins = 0;
    int64_t min_pos = INT64_MAX;
    std::string chrom;
    
    for (size_t idx : indices) {
        const MergeVariant& v = vars[idx];
        chrom = v.chrom;
        min_pos = std::min(min_pos, v.pos);
        if (v.type == VarType::DEL) total_del += v.length;
        else if (v.type == VarType::INS) total_ins += v.length;
    }
    
    MergeVariant merged;
    merged.chrom = chrom;
    merged.pos = min_pos;
    merged.merge_count = static_cast<int32_t>(indices.size());
    
    if (total_del > 0 && total_ins > 0) {
        merged.type = VarType::COMPLEX;
        merged.length = std::abs(total_del - total_ins);
        merged.ref_seq = "N";
        merged.alt_seq = "<CMP>";
    } else if (total_del > 0) {
        merged.type = VarType::DEL;
        merged.length = total_del;
        merged.ref_seq = "N";
        merged.alt_seq = "<DEL>";
    } else if (total_ins > 0) {
        merged.type = VarType::INS;
        merged.length = total_ins;
        merged.ref_seq = "N";
        merged.alt_seq = "<INS>";
    } else {
        merged.type = VarType::SNP;
        merged.length = 1;
        merged.ref_seq = "N";
        merged.alt_seq = "N";
    }
    return merged;
}

std::vector<MergeVariant> mergeIndels(std::vector<MergeVariant>& variants) {
    if (variants.size() < 2) return variants;
    
    std::sort(variants.begin(), variants.end(), 
        [](const MergeVariant& a, const MergeVariant& b) { return a.pos < b.pos; });
    
    auto clusters = findClusters(variants);
    std::vector<MergeVariant> result;
    result.reserve(clusters.size());
    
    for (const auto& cluster : clusters) {
        if (cluster.size() == 1) {
            result.push_back(variants[cluster[0]]);
        } else {
            result.push_back(mergeCluster(variants, cluster));
        }
    }
    return result;
}

std::vector<VariantCall> alignAndCall(const std::string& n_seq, const std::string& t_seq, uint64_t r_start, const std::string& chrom) {
    std::vector<VariantCall> vars;
    if (n_seq.empty() || t_seq.empty()) return vars;

    auto* m = parasail_matrix_create("ACGTN", 2, -1);
    
    // Gap Penalty Tuning: High Open (25), Low Extend (2)
    int gap_open = 25;
    int gap_extend = 2;
    
    size_t min_len = std::min(n_seq.size(), t_seq.size());
    if (min_len > 5000) { gap_open = 50; }

    auto* r = parasail_nw_trace(n_seq.c_str(), n_seq.size(), t_seq.c_str(), t_seq.size(), gap_open, gap_extend, m);
    auto* c = parasail_result_get_cigar(r, n_seq.c_str(), n_seq.size(), t_seq.c_str(), t_seq.size(), m);

    if (!c || c->len == 0) {
        if (c) parasail_cigar_free(c);
        parasail_result_free(r); parasail_matrix_free(m);
        return vars;
    }

    uint64_t r_pos = r_start;
    size_t n_idx = 0, t_idx = 0;
    
    std::vector<VariantCall> raw_vars;

    for (int i = 0; i < c->len; ++i) {
        char op = parasail_cigar_decode_op(c->seq[i]);
        uint32_t len = parasail_cigar_decode_len(c->seq[i]);

        if (op == 'M' || op == '=' || op == 'X') {
            for (uint32_t j = 0; j < len; ++j) {
                if (n_seq[n_idx + j] != t_seq[t_idx + j]) {
                    raw_vars.push_back({chrom, std::string(1, n_seq[n_idx+j]), std::string(1, t_seq[t_idx+j]), "SNP", r_pos + j, 0});
                }
            }
            r_pos += len; n_idx += len; t_idx += len;
        } else if (op == 'I') { // DEL (bases in ref not in query, assuming query=n_seq/ref_seq for align purposes? Wait. alignAndCall(N, T). N is first. I means Insertion in N relative to T. So T has DELETION.)
            // Logic confirmed: I -> DEL in T.
            if (len >= 1 && n_idx + len <= n_seq.size()) {
                uint64_t v_pos = (n_idx > 0) ? r_pos - 1 : r_pos;
                std::string v_ref = "N";
                if(n_idx > 0) v_ref = n_seq.substr(n_idx - 1, 1);
                
                raw_vars.push_back({chrom, v_ref, "<DEL>", "DEL", v_pos, -(int64_t)len});
            }
            r_pos += len; n_idx += len;
        } else if (op == 'D') { // INS (bases in T not in N)
            if (len >= 1 && t_idx + len <= t_seq.size()) {
                uint64_t v_pos = r_pos;
                std::string v_ref = (n_idx < n_seq.size()) ? n_seq.substr(n_idx, 1) : "N";
                raw_vars.push_back({chrom, v_ref, "<INS>", "INS", v_pos, (int64_t)len});
            }
            t_idx += len;
        }
    }

    parasail_cigar_free(c); parasail_result_free(r); parasail_matrix_free(m);

    // 1. Detect Chatter
    std::vector<MergeVariant> merge_vars;
    for(const auto& v : raw_vars) {
        MergeVariant mv;
        mv.chrom = v.chrom;
        mv.pos = (int64_t)v.pos;
        mv.length = (int32_t)std::abs(v.svlen);
        mv.type = mapType(v.svtype);
        mv.ref_seq = v.ref_allele;
        mv.alt_seq = v.alt_allele;
        merge_vars.push_back(mv);
    }
    
    int64_t region_end = r_start + std::max(n_seq.size(), t_seq.size());
    ChatterMetrics metrics = detectChatter(merge_vars, (int64_t)r_start, region_end);
    
    if (metrics.is_chatter) {
        // Fallback: Size Based Call
        int64_t size_diff = (int64_t)n_seq.size() - (int64_t)t_seq.size();
        if (std::abs(size_diff) > 10) {
            VariantCall fallback;
            fallback.chrom = chrom;
            fallback.pos = r_start;
            fallback.ref_allele = "N";
            
            if (size_diff > 0) { // Deletion
                fallback.alt_allele = "<DEL>";
                fallback.svtype = "DEL";
                fallback.svlen = -size_diff;
            } else { // Insertion
                fallback.alt_allele = "<INS>";
                fallback.svtype = "INS";
                fallback.svlen = -size_diff; // Positive
            }
            return {fallback};
        }
    }

    // 2. Merge Indels
    auto merged_output = mergeIndels(merge_vars);
    
    std::vector<VariantCall> final_vars;
    for(const auto& v : merged_output) {
        VariantCall vc;
        vc.chrom = v.chrom;
        vc.pos = (uint64_t)v.pos;
        vc.svtype = unmapType(v.type);
        vc.ref_allele = v.ref_seq;
        vc.alt_allele = v.alt_seq;
        
        // Fix up svlen sign
        vc.svlen = (int64_t)v.length;
        if (vc.svtype == "DEL") vc.svlen = -vc.svlen;
        
        final_vars.push_back(vc);
    }
    
    // 3. Artifact Filter (Post-Merge)
    std::vector<VariantCall> filtered_vars;
    for (const auto& v : final_vars) {
        bool is_artifact = false;
        if (v.svtype == "SNP") {
             for (const auto& other : final_vars) {
                 if (other.svtype == "INS" || other.svtype == "DEL") {
                     int64_t start = other.pos;
                     int64_t end = start + std::abs(other.svlen);
                     if ((int64_t)v.pos >= start - 20 && (int64_t)v.pos <= end + 20) {
                         is_artifact = true;
                         break;
                     }
                 }
             }
        }
        if (!is_artifact) filtered_vars.push_back(v);
    }

    return filtered_vars;
}
