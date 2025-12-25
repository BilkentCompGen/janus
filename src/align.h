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

std::vector<VariantCall> alignAndCall(const std::string& n_seq, const std::string& t_seq, uint64_t r_start, const std::string& chrom);

#endif
