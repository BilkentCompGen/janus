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
    uint32_t count(uint32_t node) const { return off[node + 1] - off[node]; }
    Post get(uint32_t node, uint32_t j) const;
};

struct SvbFile {
    FILE* f = nullptr;
    std::vector<uint64_t> off;
    uint64_t nreads = 0, maxlen = 0;
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

std::vector<uint64_t> graph_sigs(const std::string& graph_bin, uint64_t& K);

bool core_bytes(const Index& ix, std::vector<SvbFile>& svb, uint32_t node,
                std::string& out, std::string& scratch);
