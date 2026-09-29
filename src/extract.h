#pragma once
#include <cstdint>
#include <vector>
#include "core.h"
#include "io.h"
#include "lps.h"

struct scr { rdcore r; lcp_soa o[2]; bits bw; uint64_t nc, *mem; uint32_t *lab, *lsz; std::vector<uint32_t> map; std::vector<char> hs; };

void cores(scr* s, const uint8_t* rec, int lev);
int cmd_extract(const char* in, const char* out, int t, int l);
