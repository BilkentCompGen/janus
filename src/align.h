#ifndef ALIGN_H
#define ALIGN_H

#include <cstdint>
#include <string>
#include <vector>

struct VariantCall {
    std::string chrom, ref_allele, alt_allele, svtype;
    uint64_t pos;
    int64_t svlen;
};

// Variant type for merging logic
enum class VarType : uint8_t { DEL = 0, INS = 1, SNP = 2, COMPLEX = 3 };

struct MergeVariant {
    std::string chrom;
    int64_t pos;
    int32_t length;
    VarType type;
    std::string ref_seq;
    std::string alt_seq;
    int32_t merge_count = 1;
    int64_t end() const { return pos + (type == VarType::DEL ? length : 1); }
};

struct ChatterMetrics {
    int32_t variant_count;
    int32_t span;
    double density;
    bool is_chatter;
};

// Main alignment function
std::vector<VariantCall> alignAndCall(const std::string& n_seq, const std::string& t_seq, uint64_t r_start, const std::string& chrom);

// Helper functions
ChatterMetrics detectChatter(const std::vector<MergeVariant>& vars, int64_t start, int64_t end, double density_threshold = 0.1);
std::vector<MergeVariant> mergeIndels(std::vector<MergeVariant>& variants);

#endif
