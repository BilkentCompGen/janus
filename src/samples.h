#pragma once
#include <cstdint>
#include <string>
#include <vector>

// The sample sheet: one normal and N tumour lines, each with a NAME.
//
// A `.svb` read carries no name of its own -- graph.cpp synthesises one from the read's index within
// its file -- so two tumour files both call their reads r0, r1, ... . `snape call` counts support by
// distinct contig name, and a contig is named after the read it was cut from, so without a per-line
// name those two lines' support silently merges into one. The sheet is the only place that name can
// come from, which is why it is a file and not a repeated flag.
//
//   name  role    path
//   N01   normal  hifi/normal.svb
//   T01   tumor   hifi/t01.svb
//   T02   tumor   hifi/t02.svb
//
// Index 0 is ALWAYS the normal; 1..N are the tumour lines in sheet order. That order is the count
// plane order in graph.bin and the sample column order in the VCF, so it is fixed on the first run
// and verified against <dir>/samples.tsv on every resume -- a reordered sheet would silently
// reinterpret an existing graph.
struct Sample { std::string name, path; uint8_t normal; };

std::vector<Sample> samples_read(const std::string& tsv);

// Write the sheet into the work dir on the first run; on a resume, verify the given sheet is the one
// the existing artifacts were built from and abort if it is not.
void samples_lock(const std::vector<Sample>& s, const std::string& tsv);
