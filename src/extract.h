#pragma once
#include <cstdint>
#include <vector>
#include "core.h"
#include "lps.h"

static const int LEV = 5; // LCP parsing level whose cores are the graph's nodes

struct scr { rdcore r; lcp_soa o[2]; uint64_t* mem; uint32_t *lab, *lsz; std::vector<uint32_t> map; std::vector<char> hs; };

void cores(scr* s, rdcore* r, const uint8_t* rec, int lev);
