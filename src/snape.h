#pragma once
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

extern "C" void extract(const char* str, int len, uint64_t offset, int target_level,
                             int* n_out, const uint64_t** starts_out,
                             const uint64_t** ends_out, const uint64_t** sigs_out);

extern "C" void extract2(const char* str, int len, int target_level,
                              int* nf, const uint64_t** stf, const uint64_t** enf, const uint64_t** sgf,
                              int* nr, const uint64_t** str_, const uint64_t** enr, const uint64_t** sgr);
uint64_t core_sig(const char* key, int len);

using cov_t = uint8_t;
constexpr cov_t COV_MAX = 255;

struct SigCount { uint64_t sig; uint32_t slot; };

struct CEdge { uint32_t src, dst; };

constexpr uint32_t SEQ_RAW = 0x80000000u;
inline uint64_t seq_bytes(uint32_t len) {
    const uint64_t n = len & ~SEQ_RAW;
    return (len & SEQ_RAW) ? n : (n + 3) / 4;
}

struct Reads {
    FILE* f = nullptr;
    bool svb = false;
    uint64_t idx = 0;
    std::string tag;

    std::vector<char> buf;
    size_t pos = 0, end = 0;
    std::string carry;
};

Reads reads_open(const std::string& path, const std::string& sample);
bool reads_next(Reads& r, std::string& seq, std::string* name);
void reads_close(Reads& r);

FILE* xopen(const std::string& p, const char* m);
void  wr(FILE* f, const void* p, size_t n);
void  rd(FILE* f, void* p, size_t n);
void  atomic_close(FILE* f, const std::string& tmp, const std::string& dst);
