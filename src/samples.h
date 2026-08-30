#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct Sample { std::string name, path; uint8_t normal; };

std::vector<Sample> samples_read(const std::string& tsv);

void samples_lock(const std::vector<Sample>& s, const std::string& tsv);
