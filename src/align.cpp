#include "align.h"
#include "utils.h"

#include <parasail.h>

static Logger align_logger("align");

std::vector<VariantCall> alignAndCall(
    const std::string& normal_seq,
    const std::string& tumor_seq,
    uint64_t ref_start,
    const std::string& chrom
) {
    std::vector<VariantCall> variants;
    
    if (normal_seq.empty() || tumor_seq.empty()) {
        return variants;
    }
    
    // Create substitution matrix (match=2, mismatch=-1)
    parasail_matrix_t* matrix = parasail_matrix_create("ACGTN", 2, -1);
    
    // Global alignment with traceback
    // s1 = normal (reference), s2 = tumor (query)
    parasail_result_t* result = parasail_nw_trace(
        normal_seq.c_str(), (int)normal_seq.size(),
        tumor_seq.c_str(), (int)tumor_seq.size(),
        3,  // gap open penalty
        1,  // gap extend penalty
        matrix
    );
    
    if (!result) {
        parasail_matrix_free(matrix);
        return variants;
    }
    
    // Get CIGAR string
    parasail_cigar_t* cigar = parasail_result_get_cigar(
        result,
        normal_seq.c_str(), (int)normal_seq.size(),
        tumor_seq.c_str(), (int)tumor_seq.size(),
        matrix
    );
    
    if (!cigar || cigar->len == 0) {
        parasail_result_free(result);
        parasail_matrix_free(matrix);
        return variants;
    }
    
    // Parse CIGAR to find insertions and deletions
    // CIGAR operations: M=match/mismatch, I=insertion, D=deletion
    uint64_t ref_pos = ref_start;  // Position in reference/normal (1-based)
    size_t normal_idx = 0;         // Index into normal_seq (0-based)
    size_t tumor_idx = 0;          // Index into tumor_seq (0-based)
    
    for (int i = 0; i < cigar->len; ++i) {
        uint32_t cigar_int = cigar->seq[i];
        char op = parasail_cigar_decode_op(cigar_int);
        uint32_t len = parasail_cigar_decode_len(cigar_int);
        
        switch (op) {
            case 'M':
            case '=':
            case 'X':
                // Match or mismatch - advance both positions
                ref_pos += len;
                normal_idx += len;
                tumor_idx += len;
                break;
                
            case 'I':
                // In parasail CIGAR with s1=normal, s2=tumor:
                // I = bases in s1 (normal) not in s2 (tumor) = DELETION in tumor
                // This is the opposite of standard SAM CIGAR convention!
                if (len >= 1 && normal_idx + len <= normal_seq.size()) {
                    VariantCall var;
                    var.chrom = chrom;
                    
                    // For VCF format, include anchor base before deletion
                    if (normal_idx > 0) {
                        var.pos = ref_pos - 1;
                        var.ref_allele = normal_seq.substr(normal_idx - 1, len + 1);
                        var.alt_allele = normal_seq.substr(normal_idx - 1, 1);
                    } else {
                        var.pos = ref_pos;
                        var.ref_allele = normal_seq.substr(normal_idx, len);
                        var.alt_allele = ".";
                    }
                    
                    var.svtype = "DEL";
                    var.svlen = -(int64_t)len;
                    variants.push_back(var);
                }
                ref_pos += len;
                normal_idx += len;
                break;
                
            case 'D':
                // In parasail CIGAR with s1=normal, s2=tumor:
                // D = bases in s2 (tumor) not in s1 (normal) = INSERTION in tumor
                if (len >= 1 && tumor_idx + len <= tumor_seq.size()) {
                    VariantCall var;
                    var.chrom = chrom;
                    
                    if (normal_idx > 0 && normal_idx <= normal_seq.size()) {
                        var.pos = ref_pos - 1;
                        var.ref_allele = normal_seq.substr(normal_idx - 1, 1);
                        var.alt_allele = normal_seq.substr(normal_idx - 1, 1) + tumor_seq.substr(tumor_idx, len);
                    } else {
                        var.pos = ref_pos;
                        var.ref_allele = "N";
                        var.alt_allele = "N" + tumor_seq.substr(tumor_idx, len);
                    }
                    
                    var.svtype = "INS";
                    var.svlen = (int64_t)len;
                    variants.push_back(var);
                }
                tumor_idx += len;
                break;
                
            default:
                // S (soft clip), H (hard clip), N (skip) - shouldn't happen in global alignment
                break;
        }
    }
    
    // Cleanup
    parasail_cigar_free(cigar);
    parasail_result_free(result);
    parasail_matrix_free(matrix);
    
    align_logger.debug("Found " + std::to_string(variants.size()) + " variants from alignment");
    
    return variants;
}
