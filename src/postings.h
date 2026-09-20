#pragma once
#include "snape.h"
#include <cstdint>
#include <string>
#include <vector>

struct Post {
    uint32_t read, ord, st, clen, dl, dr;
    uint8_t samp;
};

struct Index {
    uint32_t pad = 0, N = 0;
    uint64_t K = 0;
    std::vector<uint32_t> off;
    std::vector<uint64_t> w;
    void load(const std::string& path);
    uint32_t count(uint32_t node) const;
    Post get(uint32_t node, uint32_t j) const;
};

struct SvbFile {
    FILE* f = nullptr;
    std::vector<uint64_t> off;
    uint64_t nreads = 0;
    void open(const std::string& svb);
    void close();
    bool read(uint64_t no, std::string& seq);
};

struct SigLookup {
    static constexpr int BITS = 24;
    std::vector<uint32_t> off;
    const uint64_t* v = nullptr;
    void build(const std::vector<uint64_t>& sig);
    uint32_t find(uint64_t s) const;
};

// Opens a .svb, validates its 32-byte header ("SVB1" magic + version), and
// leaves the file positioned right after it, ready for svb_read_seq(). The
// same 32-byte layout (magic[4] + ver[4] + nreads[8] + total[8] + reserved[8])
// is written by `snape pack` and walked by `snape svbidx` -- keep all three
// in sync if it ever changes.
FILE* open_svb(const std::string& path);
bool svb_read_seq(FILE* f, std::string& seq);
std::vector<uint64_t> graph_sigs(const std::string& graph_bin, uint64_t& K);
bool core_bytes(const Index& ix, std::vector<SvbFile>& svb, uint32_t node,
                std::string& out, std::string& scratch);
