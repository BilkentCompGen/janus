#pragma once
#include "samples.h"
#include "snape.h"
#include <cstdint>
#include <string>
#include <vector>

void run_graph_reads(const std::vector<Sample>& samples, const std::string& dir,
                     int level, int margin, int min_mult, int threads, uint64_t max_bases = 0,
                     int index_n = 8, int pad_cores = 50);
