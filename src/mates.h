#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct mate { uint32_t L, R, f, site; uint64_t best[3]; std::string seg; };

void find_mates(std::vector<mate>& w, const std::vector<uint64_t>& g, const char* reads, int t);
