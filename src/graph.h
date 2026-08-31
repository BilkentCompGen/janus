#pragma once
#include "samples.h"
#include "snape.h"
#include <cstdint>
#include <string>
#include <vector>

void run_graph_reads(const std::vector<Sample>& samples, const std::string& dir,
                     int level, int margin, int min_mult, int threads, int index_n = 0, int pad_cores = 50);
