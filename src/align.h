#ifndef ALIGN_H
#define ALIGN_H

#include <cstdint>
#include <string>
#include <vector>

struct VariantCall {
    std::string chrom;
    uint64_t pos;           // 1-based VCF position
    std::string ref_allele;
    std::string alt_allele;
    std::string svtype;     // DEL or INS
    int64_t svlen;          // Negative for DEL, positive for INS
};

// Align normal and tumor sequences and extract variants
// normal_seq: the reference/normal sequence (longer in case of deletion)
// tumor_seq: the tumor sequence (shorter in case of deletion)
// ref_start: reference coordinate of the bubble start (1-based)
// chrom: chromosome name
std::vector<VariantCall> alignAndCall(
    const std::string& normal_seq,
    const std::string& tumor_seq,
    uint64_t ref_start,
    const std::string& chrom
);

#endif // ALIGN_H
