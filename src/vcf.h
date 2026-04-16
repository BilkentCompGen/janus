#ifndef VCF_H
#define VCF_H

#include <string>
#include <vector>

#include "sv.h"

void writeVCF(const std::string& vcf_path,
              std::vector<SVCall>& calls,
              const std::string& ref_sequence,
              const std::string& chr_str);

#endif
